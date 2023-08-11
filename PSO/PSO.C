//ROOT headers
#include <TFile.h>
#include <TROOT.h>
#include <TTree.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TF1.h>


//Cpp headers
#include <iostream>
#include <fstream> //For log txt
#include <vector>
#include <random>
#include <cmath>
#include <chrono> //Timer
using namespace std;

// CUSTOM STREAM BUFFER CLASS FOR OUTPUT LOG FILE
class DualStreamBuffer : public std::streambuf
{
public:
  // Constructor
  DualStreamBuffer(std::streambuf *primaryBuffer, std::streambuf *secondaryBuffer)
      : primaryBuffer(primaryBuffer), secondaryBuffer(secondaryBuffer) {}


  int_type overflow(int_type c) override
  {
    if (c != EOF)
    {
      if (primaryBuffer != nullptr)
        primaryBuffer->sputc(c);

      if (secondaryBuffer != nullptr)
        secondaryBuffer->sputc(c);
    }
    return c;
  }



private:
  std::streambuf *primaryBuffer;   // Primary stream buffer (std::cout)
  std::streambuf *secondaryBuffer; // Secondary stream buffer (log file)
};
// END OF LOG FILE CLASS


struct particle
{
  double position[4];
  double pbest[4];
  double velocity[4];
  double signal[2];
  double significance;
};


void print_particle(const particle& p, double *gbest, double gbest_significance, double *max_significance, double *max_signal, double *max_signal_er, int i, int j)
{
  cout << "    --------------------------------------------------------------------------" << endl;
  cout << "                                 Particle (" << i << ", " << j << "):" << endl << endl; 
  cout << "    Position:                        (" << p.position[0] << ", " << p.position[1] << ", " << p.position[2] << ", " << p.position[3] << ")" << endl  << endl;
  cout << "    Best Position:                   (" << p.pbest[0] << ", " << p.pbest[1] << ", " << p.pbest[2] << ", " << p.pbest[3] << ")" << endl  << endl;
  cout << "    Velocity:                        (" << p.velocity[0] << ", " << p.velocity[1] << ", " << p.velocity[2] << ", " << p.velocity[3] << ")" << endl  << endl;
  cout << "    Current Signal:                  " << p.signal[0] << " +- " << p.signal[1] << endl << endl;
  cout << "    Max Signal:                      " << max_signal[j] << " +- " << max_signal_er[j] << endl << endl;
  cout << "    Current Significance:            " << p.significance << endl << endl;
  cout << "    Max Significance:                " << max_significance[j] << endl << endl;   
  cout << "    Global Best:                     (" << gbest[0] << ", " << gbest[1] << ", " << gbest[2] << ", " << gbest[3] << ")" << endl  << endl;
  cout << "    Global Max Significance:         " << gbest_significance << endl << endl;   
  cout << "    --------------------------------------------------------------------------" << endl << endl;
}


//Apply boundaries on particles 
void apply_bounds(vector<vector<particle>> &swarm, Float_t dLepR_bounds[2], Float_t dMetZPhi_bounds[2], Float_t met_tst_bounds[2], Float_t MetOHT_bounds[2], int i, int j, int k)
{
  Float_t bounds[4][2] = {
      {dLepR_bounds[0], dLepR_bounds[1]},
      {dMetZPhi_bounds[0], dMetZPhi_bounds[1]},
      {met_tst_bounds[0], met_tst_bounds[1]},
      {MetOHT_bounds[0], MetOHT_bounds[1]}};

  if (swarm[i][j].position[k] < bounds[k][0])
  {
    swarm[i][j].position[k] = bounds[k][0];
  }
  else if (swarm[i][j].position[k] > bounds[k][1])
  {
    swarm[i][j].position[k] = bounds[k][1];
  }
}

//Update the search space particles before entering next iteration
void update_particle(vector<vector<particle>> &swarm, double *gbest, int i, int iterations, int n_particles)
{
  uniform_real_distribution<> uni_dist(0., 1.);
  random_device rd;
  mt19937 gen(rd());

  float_t w_min = 0.6;
  float_t w_max = 0.9;
  float_t w = w_max - (w_max - w_min) * i / iterations;
  // float_t w = 0.7;
  float_t c1 = 1.47;   
  float_t c2 = 1.47;
  float_t r1 = uni_dist(gen);
  float_t r2 = uni_dist(gen);
  //Clerc p.40 - (w, c1, c2) = (0.7, 1.47, 1.47) or (0.6, 1.62, 1.62)

  for (int j = 0; j < n_particles; j++)
  {
    for (int k = 0; k < 4; k++)
    {
      swarm[i][j].velocity[k] = w * swarm[i - 1][j].velocity[k] + c1 * r1 * (swarm[i - 1][j].pbest[k] - swarm[i - 1][j].position[k]) + c2 * r2 * (gbest[k] - swarm[i - 1][j].position[k]);
      swarm[i][j].position[k] = swarm[i - 1][j].position[k] + swarm[i][j].velocity[k];


      //Boundaries
      Float_t dLepR_bounds[2] = {1., 3.};
      Float_t dMetZPhi_bounds[2] = {2., 4.};
      Float_t met_tst_bounds[2] = {70, 130};
      Float_t MetOHT_bounds[2] = {0.2, 1.4};

      apply_bounds(swarm, dLepR_bounds, dMetZPhi_bounds, met_tst_bounds, MetOHT_bounds, i, j, k);

    }
  }
  return;
}


vector<Float_t> event_counter(particle &particle, TTree *tree)
{

  TH1::SetDefaultSumw2(kTRUE);

  Int_t nentries = (Int_t)tree->GetEntries();

  Double_t M2Lep = 0.;
  Double_t met_tst = 0.;
  Double_t met_signif = 0.;
  Double_t dMetZPhi = 0.;
  // Float_t frac_pT = 0.;
  Double_t MetOHT = 0.;
  Double_t dPhiJ100met = 0.;
  Double_t dLepR = 0.;
  Double_t n_bjets = 0.;
  Double_t leading_pT_lepton = 0;
  Double_t subleading_pT_lepton = 0;
  Double_t Z_pT = 0;
  Double_t n_jets = 0.;
  Double_t detajj = 0;
  Double_t mjj = 0;
  Double_t leading_jet_pt = 0;
  Double_t second_jet_pt = 0;
  Double_t event_3CR = 0.;
  Double_t event_type = 0.;
  Double_t weight;
  vector<Float_t> events;
  events.clear();

  vector<TString> branches = {
      "M2Lep",
      "met_tst",
      "dMetZPhi",
      "MetOHT",
      "dLepR",
      "M2Lep",
      // "leading_pT_lepton",
      // "subleading_pT_lepton",
      // "Z_pT",
      // "n_jets",
      // "n_bjets",
      // "detajj",
      // "mjj",
      // "leading_jet_pt",
      // "second_jet_pt",
      "event_3CR",
      "event_type",
      "global_weight"};


  tree->SetBranchStatus("*", 0);

  for (const auto& branch : branches) {
    tree->SetBranchStatus(branch, 1);
  }

  tree->SetBranchAddress("M2Lep", &M2Lep);
  tree->SetBranchAddress("met_tst", &met_tst);
  tree->SetBranchAddress("met_signif", &met_signif);
  tree->SetBranchAddress("dMetZPhi", &dMetZPhi);
  tree->SetBranchAddress("MetOHT", &MetOHT);
  tree->SetBranchAddress("dLepR", &dLepR);
  tree->SetBranchAddress("M2Lep", &M2Lep);
  tree->SetBranchAddress("global_weight", &weight);

  Double_t signal = 0.;
  Double_t signaler = 0.;

  // Loop over events
  for (Int_t i = 0; i < nentries; i++)
  {
    tree->GetEntry(i);

    if (event_3CR==0 && (event_type==0 || event_type==1) && met_tst>90 && dLepR<1.8 && n_bjets==0 && dMetZPhi>2.6 &&
    dLepR < particle.position[0] && dMetZPhi > particle.position[1] && met_tst > particle.position[2] && MetOHT > particle.position[3])
    {
     
    // if (event_3CR==0 && (event_type==0 || event_type==1) && dLepR < particle.position[0] && dMetZPhi > particle.position[1] && met_tst > particle.position[2] && MetOHT > particle.position[3])
    // {
      signal = signal + weight;
      signaler = signaler + weight * weight;
    }
  }

  cout << "    ENTRIES = " << tree->GetEntries() << endl << endl;
  cout << "          N = " << signal << " +- " << sqrt(signaler) << endl << endl;

  events.push_back(signal);
  events.push_back(sqrt(signaler));

  return events;
}


//Main
void PSO()
{
  // Timer start
  auto start = std::chrono::high_resolution_clock::now();

  // Output log file
  ofstream logFile("./PSO.txt");
  DualStreamBuffer dualBuffer(std::cout.rdbuf(), logFile.rdbuf());
  std::streambuf *oldBuffer = std::cout.rdbuf(&dualBuffer);

  gROOT->SetBatch(kTRUE); //Disable plot popups
  

  //Load root files
  string filepath = "../../data/SAMPLES/SR/";

  TFile *file_data = new TFile((string(filepath) + "DATA.root").c_str());
  TTree *tree_data = file_data->Get<TTree>("tree");

  // Singal = SignalEWK + SignalQCD
  TFile *file_llvv = new TFile((string(filepath) + "llvv.root").c_str());
  TTree *tree_llvv = file_llvv->Get<TTree>("tree");
  TFile *file_llvvjj = new TFile((string(filepath) + "llvvjj.root").c_str());
  TTree *tree_llvvjj = file_llvvjj->Get<TTree>("tree");

  // WZ
  TFile *file_WZ = new TFile((string(filepath) + "WZ.root").c_str());
  TTree *tree_WZ = file_WZ->Get<TTree>("tree");

  // Zjets = Z_jets_ee + Z_jets_mumu
  TFile *file_Z_jets_ee = new TFile((string(filepath) + "Z_jets_ee.root").c_str());
  TTree *tree_Z_jets_ee = file_Z_jets_ee->Get<TTree>("tree");
  TFile *file_Z_jets_mumu = new TFile((string(filepath) + "Z_jets_mumu.root").c_str());
  TTree *tree_Z_jets_mumu = file_Z_jets_mumu->Get<TTree>("tree");

  // top = top + ttbarV_ttbar_VV + Wt
  TFile *file_top = new TFile((string(filepath) + "top.root").c_str());
  TTree *tree_top = file_top->Get<TTree>("tree");
  TFile *file_ttbarV_ttbarVV = new TFile((string(filepath) + "ttbarV_ttbarVV.root").c_str());
  TTree *tree_ttbarV_ttbarVV = file_ttbarV_ttbarVV->Get<TTree>("tree");
  TFile *file_Wt = new TFile((string(filepath) + "Wt.root").c_str());
  TTree *tree_Wt = file_Wt->Get<TTree>("tree");

  // WW
  TFile *file_WW = new TFile((string(filepath) + "WW.root").c_str());
  TTree *tree_WW = file_WW->Get<TTree>("tree");

  // Othr
  TFile *file_llll = new TFile((string(filepath) + "llll.root").c_str());
  TTree *tree_llll = file_llll->Get<TTree>("tree");
  TFile *file_llqq = new TFile((string(filepath) + "llqq.root").c_str());
  TTree *tree_llqq = file_llqq->Get<TTree>("tree");
  TFile *file_VVV = new TFile((string(filepath) + "VVV.root").c_str());
  TTree *tree_VVV = file_VVV->Get<TTree>("tree");
  TFile *file_W_jets = new TFile((string(filepath) + "W_jets.root").c_str());
  TTree *tree_W_jets = file_W_jets->Get<TTree>("tree");
  TFile *file_Ztt = new TFile((string(filepath) + "Ztt.root").c_str());
  TTree *tree_Ztt = file_Ztt->Get<TTree>("tree");
  TFile *file_WZ_jj = new TFile((string(filepath) + "WZ_jj.root").c_str());
  TTree *tree_WZ_jj = file_WZ_jj->Get<TTree>("tree");
  TFile *file_lllljj = new TFile((string(filepath) + "lllljj.root").c_str());
  TTree *tree_lllljj = file_lllljj->Get<TTree>("tree");
  TFile *file_llvvjj_WW = new TFile((string(filepath) + "llvvjj_WW.root").c_str());
  TTree *tree_llvvjj_WW = file_llvvjj_WW->Get<TTree>("tree");


  //PSO ALGORITHM

  int iterations = 30;
  int n_particles = 3;

  double max_signal[n_particles];
  double max_signal_er[n_particles];
  double max_significance[n_particles];
  double gbest[4];
  double gbest_significance = - 1;
  vector<double> gbest_vector;

  for (int index = 0; index < 4; index++) {gbest[index] = -1;}
  for (int index = 0; index < n_particles; index++) {max_significance[index] = -1;}
  for (int index = 0; index < n_particles; index++) {max_signal[index] = -1; max_signal_er[index] = -1;}
  

  // Initialization
  uniform_real_distribution<> uni_dLepR(1., 3. );
  uniform_real_distribution<> uni_dMetZPhi(2., 4.);
  uniform_real_distribution<> uni_met_tst(70, 130);
  uniform_real_distribution<> uni_MetOHT(0.2, 1.4);

  uniform_real_distribution<> uni_dLepR2(-2., 4. );
  uniform_real_distribution<> uni_dMetZPhi2(-2., 6.);
  uniform_real_distribution<> uni_met_tst2(-50, 200);
  uniform_real_distribution<> uni_MetOHT2(-1.2, 1.8);

  

  random_device rd;
  mt19937 gen(rd());

  

  
  vector<vector<particle>> swarm(iterations, vector<particle>(n_particles)); //Define the swarm

  for (int i = 0; i < iterations; i++)
  {
    if (i != 0)
    {
      for (int j = 0; j < n_particles; j++)
      {
        particle &particle = swarm[i][j]; // Define the particle
        update_particle(swarm, gbest, i, iterations, n_particles);
      }
    }
    for (int j = 0; j < n_particles; j++)
    {
  
      particle &particle = swarm[i][j]; // Define the particle
  
      // cout << "   SWARM:  "  << swarm[i] << endl << endl;
      cout << "  ____________________________________________________" << endl << endl;
      cout << "              ITERATION:  " << i << "   PARTICLE:  " << j << "   " << endl;
      cout << "  ____________________________________________________" << endl << endl << endl;
      
      cout << "    ================== DATA ==================    " << endl << endl;
      cout << "    DATA:";
      vector<Float_t> n_data = event_counter(particle, tree_data);

      while (n_data[0] < 1) // Re-initialize until the phase space produces non-zero events
      {
        // Initialize position
        particle.position[0] = uni_dLepR(gen);
        particle.position[1] = uni_dMetZPhi(gen);
        particle.position[2] = uni_met_tst(gen);
        particle.position[3] = uni_MetOHT(gen);

        // Initialize pbest
        particle.pbest[0] = particle.position[0];
        particle.pbest[1] = particle.position[1];
        particle.pbest[2] = particle.position[2];
        particle.pbest[3] = particle.position[3];

        // Initialize velocity
        particle.velocity[0] = uni_dLepR2(gen);
        particle.velocity[1] = uni_dMetZPhi2(gen);
        particle.velocity[2] = uni_met_tst2(gen);
        particle.velocity[3] = uni_MetOHT2(gen);
        cout << "          No entries - Re-Initialize..." << endl << endl;
        cout << "    DATA:";
        n_data = event_counter(particle, tree_data);
      }

      cout << "    ================== SIGNAL ==================    " << endl << endl;
      cout << "    llvv:";
      vector<Float_t> n_llvv = event_counter(particle, tree_llvv);
      cout << "    llvvjj:";
      vector<Float_t> n_llvvjj = event_counter(particle, tree_llvvjj);
  
      cout << "    ================== WZ ==================    " << endl << endl;
      cout << "    WZ:";
      vector<Float_t> n_WZ = event_counter(particle, tree_WZ);
  
      cout << "    ================== Zjets ==================    " << endl << endl;
      cout << "    Z_jets_ee:";
      vector<Float_t> n_Zjets_ee = event_counter(particle, tree_Z_jets_ee);
      cout << "    Z_jets_mumu:";
      vector<Float_t> n_Zjets_mumu = event_counter(particle, tree_Z_jets_mumu);
  
      cout << "    ================== top ==================    " << endl << endl;
      cout << "    Top:";
      vector<Float_t> n_top = event_counter(particle, tree_top);
      cout << "    ttbarV_ttbarVV:";
      vector<Float_t> n_ttbarV_ttbarVV = event_counter(particle, tree_ttbarV_ttbarVV);
      cout << "    Wt:";
      vector<Float_t> n_Wt = event_counter(particle, tree_Wt);
  
      cout << "    ================== WW ==================    " << endl << endl;
      cout << "    WW:";
      vector<Float_t> n_WW = event_counter(particle, tree_WW);
  
      cout << "    ================== Othr ==================    " << endl << endl;
      cout << "    llll:";
      vector<Float_t> n_llll = event_counter(particle, tree_llll);
      cout << "    llqq:";
      vector<Float_t> n_llqq = event_counter(particle, tree_llqq);
      cout << "    VVV:";
      vector<Float_t> n_VVV = event_counter(particle, tree_VVV);
      cout << "    W_jets:";
      vector<Float_t> n_Wjets = event_counter(particle, tree_W_jets);
      cout << "    Ztt:";
      vector<Float_t> n_Ztt = event_counter(particle, tree_Ztt);
      cout << "    WZ_jj:";
      vector<Float_t> n_WZjj = event_counter(particle, tree_WZ_jj);
      cout << "    lllljj:";
      vector<Float_t> n_lllljj = event_counter(particle, tree_lllljj);
      cout << "    llvvjj_WW:";
      vector<Float_t> n_llvvjj_WW = event_counter(particle, tree_llvvjj_WW);
  
  
      Float_t events_data = n_data[0];
      Float_t events_data_er = n_data[1];
      Float_t events_signal = n_llvv[0] + n_llvvjj[0];
      Float_t events_signal_er = sqrt(pow(n_llvv[1], 2) + pow(n_llvvjj[1], 2));
      Float_t events_bkg;
      Float_t events_bkg_er;
      Float_t events_WZ = n_WZ[0];
      Float_t events_WZ_er = n_WZ[1];
      Float_t events_Zjets = n_Zjets_ee[0] + n_Zjets_mumu[0];
      Float_t events_Zjets_er = sqrt(pow(n_Zjets_ee[1], 2) + pow(n_Zjets_mumu[1], 2));
      Float_t events_top = n_top[0] + n_ttbarV_ttbarVV[0] + n_Wt[0];
      Float_t events_top_er = sqrt(pow(n_top[1],2) + pow(n_ttbarV_ttbarVV[1], 2) + pow(n_Wt[1], 2));
      Float_t events_WW = n_WW[0];
      Float_t events_WW_er = n_WW[1];
      Float_t events_othr = n_llll[0] + n_llqq[0] + n_VVV[0] + n_Wjets[0] + n_Ztt[0] + n_WZjj[0] + n_lllljj[0] + n_llvvjj_WW[0];
      Float_t events_othr_er = sqrt(pow(n_llll[1], 2) + pow(n_llqq[1], 2) + pow(n_VVV[1], 2) + pow(n_Wjets[1], 2) + pow(n_Ztt[1], 2) + pow(n_WZjj[1], 2) + pow(n_lllljj[1], 2) + pow(n_llvvjj_WW[1], 2));
  
      events_bkg = events_WZ + events_top + events_WW + events_Zjets + events_othr;
      events_bkg_er = sqrt(pow(events_WZ_er, 2) + pow(events_top_er, 2) + pow(events_WW_er, 2) + pow(events_Zjets, 2) + pow(events_othr_er, 2));
  
  
      Float_t S = events_signal;
      Float_t B = events_bkg;
      particle.significance = sqrt(2 * ((S + B) * log(1 + (S / B)) - S));
      particle.signal[0] = events_signal;
      particle.signal[1] = events_signal_er;


      
      if (particle.significance > max_significance[j])
      { 
        max_significance[j] = particle.significance;
        max_signal[j] = particle.signal[0];
        max_signal_er[j] = particle.signal[1];
  
        for (int index = 0; index < 4; ++index)
        {
          particle.pbest[index] = particle.position[index];
        }

        if (particle.significance > gbest_significance)
        {
          gbest_significance = particle.significance;

          for (int index = 0; index < 4; ++index)
          {
            gbest[index] = particle.pbest[index];
          }
        }
      }
      else
      {
        for (int index = 0; index < 4; ++index)
        {
          swarm[i][j].pbest[index] = swarm[i - 1][j].pbest[index];
        }
      }
      print_particle(particle, gbest, gbest_significance, max_significance, max_signal, max_signal_er, i, j);
    }
    gbest_vector.push_back(gbest_significance);

    // Graph significance vs iterations
    TCanvas *c = new TCanvas("c", "Significance vs Iterations", 800, 600);
    TGraph *graph = new TGraph();
    graph->SetTitle("Significance vs Iterations");
    graph->GetXaxis()->SetTitle("Iterations");
    graph->GetYaxis()->SetTitle("Significance");
    graph->GetXaxis()->CenterTitle();
    graph->GetYaxis()->CenterTitle();
    graph->GetXaxis()->SetTitleOffset(1.2);
    graph->GetYaxis()->SetTitleOffset(1.3);
    graph->GetXaxis()->SetLabelSize(0.04);
    graph->GetYaxis()->SetLabelSize(0.04);
    graph->GetXaxis()->SetRangeUser(0, gbest_vector.size() + 1);

    // Use a single distinctive marker style
    graph->SetMarkerStyle(29);
    graph->SetMarkerSize(2);
    graph->SetLineColor(kBlack);

    for (int i = 0; i < gbest_vector.size(); ++i)
    {
      graph->SetPoint(i, i + 1, gbest_vector[i]);
    }

    graph->Draw("ALP");

    double current_best = gbest_vector[gbest_vector.size() - 1]; 
    double x1 = 0.12, y1 = 0.75;                      
    double x2 = 0.37, y2 = 0.9;                         

    // Adjust legend size to fit the parameter value
    TLegend *legend = new TLegend(x1, y1, x2, y2);
    legend->SetMargin(0.15); 
    legend->SetTextSize(0.03);
    legend->AddEntry((TObject *)0, Form("Current Best: %.3f", current_best), "");
    legend->Draw();

    c->SetGrid();
    c->Draw();
    c->SaveAs("significance_vs_iterations.png");
  }

  cout << "   --------------------------------------------------------------------------" << endl << endl;
  cout << "   Max significance is:  "  << gbest_significance << endl;
  cout << "   Best position:        (" << gbest[0] << ", " << gbest[1] << ", " << gbest[2] << ", " << gbest[3] << ") " << endl;
  cout << "   --------------------------------------------------------------------------" << endl << endl;

  

  // Timer stop
  auto end = chrono::high_resolution_clock::now();
  chrono::duration<float> duration = end - start;
  cout << endl << "   Script executed in: " << int(duration.count() / 60.0) << " minutes" << " and " << int((duration.count() / 60.0 - int(duration.count() / 60.0))*60) << " s" <<  endl << endl;


  // For the log file
  std::cout.rdbuf(oldBuffer); // Restore the original stream buffer for std::cout
  logFile.close(); // Close the log file

  return;
}
