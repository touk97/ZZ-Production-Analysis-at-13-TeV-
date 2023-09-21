//Particle Swarm optimization algorithm for background estimation - 3 variables optimization

//ROOT headers
#include <TFile.h>
#include <TROOT.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TF1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TBrowser.h>


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
  std::streambuf *primaryBuffer;  
  std::streambuf *secondaryBuffer; 
};
// END OF LOG FILE CLASS


struct particle
{
  Float_t position[3];
  Float_t pbest[3];
  Float_t velocity[3];
  Float_t signal[2];
  Float_t significance;
};


void print_particle(const particle& p, Float_t *gbest, Int_t is_new_gbest, Int_t is_new_pbest, Float_t gbest_significance, Float_t *max_significance, Float_t *max_signal, Float_t *max_signal_er, int i, int j)
{
  cout << "    ------------------------------------------------------------------------------" << endl;
  cout << "                                   Particle (" << i + 1 << ", " << j + 1 << "):" << endl << endl; 
  cout << "    Position:                            (" << p.position[0] << ", " << p.position[1] << ", " << p.position[2] << ")" << endl  << endl;
  cout << "    Best Position:                       (" << p.pbest[0] << ", " << p.pbest[1] << ", " << p.pbest[2] << ")" << endl  << endl;
  cout << "    Best Position has been updated:      "  << is_new_pbest << endl  << endl;
  cout << "    Velocity:                            (" << p.velocity[0] << ", " << p.velocity[1] << ", " << p.velocity[2] <<  ")" << endl  << endl;
  cout << "    Current Signal:                      "  << p.signal[0] << " +- " << p.signal[1] << endl << endl;
  cout << "    Max Signal:                          "  << max_signal[j] << " +- " << max_signal_er[j] << endl << endl;
  cout << "    Current Significance:                "  << p.significance << endl << endl;
  cout << "    Max Significance:                    "  << max_significance[j] << endl << endl;   
  cout << "    Global Best:                         (" << gbest[0] << ", " << gbest[1] << ", " << gbest[2] << ")" << endl  << endl;
  cout << "    Global Best has been updated:        "  << is_new_gbest << endl  << endl;
  cout << "    Global Max Significance:             "  << gbest_significance << endl << endl;   
  cout << "    ----------------------------------------------------------------------------" << endl << endl;
}


//Apply boundaries on particles 
void apply_bounds(vector<vector<particle>> &swarm, Float_t bounds[3][2], int i, int j, int k)
{
  uniform_real_distribution<> uni_dist(0, 1);
  random_device rd;
  mt19937 gen(rd());

  if (swarm[i][j].position[k] <= bounds[k][0])
  {
    swarm[i][j].position[k] = bounds[k][0];
    swarm[i][j].velocity[k] = -swarm[i][j].velocity[k] * uni_dist(gen);
  }
  else if (swarm[i][j].position[k] > bounds[k][1])
  {
    swarm[i][j].position[k] = bounds[k][1];
    swarm[i][j].velocity[k] = -swarm[i][j].velocity[k] * uni_dist(gen);
  }
}

//Update the search space particles before entering next iteration
void update_swarm(vector<vector<particle>> &swarm, Float_t bounds[3][2], Float_t *gbest, int i, int iterations, int n_particles)
{
  uniform_real_distribution<> uni_dist(0., 1.);
  random_device rd;
  mt19937 gen(rd());

  //Big w values result in many entries in boundaries

  // float_t w_min = 0.4;
  // float_t w_max = 0.7;
  // float_t w = w_max - (w_max - w_min) * i /  iterations;
  float_t w = 0.8;
  float_t cmax = 1.47;   
  float_t r1 = uni_dist(gen);
  float_t r2 = uni_dist(gen);
  //Clerc p.40 - (w, cmax) = (0.7, 1.47) or (0.6, 1.62)

  for (int j = 0; j < n_particles; j++)
  {
    for (int k = 0; k < 4; k++)
    {
      swarm[i+1][j].velocity[k] = w * swarm[i][j].velocity[k] + cmax * r1 * (swarm[i][j].pbest[k] - swarm[i][j].position[k]) + cmax * r2 * (gbest[k] - swarm[i][j].position[k]);
      swarm[i+1][j].position[k] = swarm[i][j].position[k] + swarm[i+1][j].velocity[k];

      apply_bounds(swarm, bounds, i+1, j, k);
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
      // "M2Lep",
      "met_tst",
      "dMetZPhi",
      "MetOHT",
      // "dLepR",
      // "leading_pT_lepton",
      // "subleading_pT_lepton",
      // "Z_pT",
      // "n_jets",
      // "n_bjets",
      // "detajj",
      // "mjj",
      // "leading_jet_pt",
      // "second_jet_pt",
      // "event_3CR",
      // "event_type",
      "global_weight"};


  tree->SetBranchStatus("*", 0);

  for (const auto& branch : branches) {
    tree->SetBranchStatus(branch, 1);
  }

  // tree->SetBranchAddress("M2Lep", &M2Lep);
  tree->SetBranchAddress("met_tst", &met_tst);
  tree->SetBranchAddress("met_signif", &met_signif);
  tree->SetBranchAddress("dMetZPhi", &dMetZPhi);
  tree->SetBranchAddress("MetOHT", &MetOHT);
  tree->SetBranchAddress("global_weight", &weight);

  Double_t signal = 0.;
  Double_t signaler = 0.;

  // Loop over events
  for (Int_t i = 0; i < nentries; i++)
  {
    tree->GetEntry(i);
    
    //Any additional criteria is already accounted for during the sample files generation
    //FID.VOLUME: 80 < M2Lep < 100, met_tst > 70, dLepR < 1.8, dMetZPhi > 2.2
    if (dMetZPhi > particle.position[0] && met_tst > particle.position[1] && MetOHT > particle.position[2])
    {
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
void PSO_3VAR()
{

  gROOT->SetBatch(kTRUE); // Disable plot popups

  // Timer start
  auto start = std::chrono::high_resolution_clock::now();

  // Output log file
  ofstream logFile("./3var/PSO_3VAR.txt");
  DualStreamBuffer dualBuffer(std::cout.rdbuf(), logFile.rdbuf());
  std::streambuf *oldBuffer = std::cout.rdbuf(&dualBuffer);

  //------------PSO ALGORITHM------------//

  int iterations = 200;
  int n_particles = 50;

  Float_t max_signal[n_particles];
  Float_t max_signal_er[n_particles];
  Float_t max_significance[n_particles];
  Float_t gbest[4];
  Int_t   is_new_gbest = 0;
  Int_t   is_new_pbest = 0;
  Float_t gbest_significance = -1;
  vector<Float_t> gbest_vector;

  //Search space Boundaries {min, max}
  Float_t bounds[3][2] = 
  {  
    {2.2, 3.0},          // dMetZPhi bounds
    {90.0, 110.0},       // met_tst bounds
    {0.5, 0.85}            // MetOHT bounds
  };

  for (int index = 0; index < 4; index++) {gbest[index] = -1;}
  for (int index = 0; index < n_particles; index++) {max_significance[index] = -1;}
  for (int index = 0; index < n_particles; index++) {max_signal[index] = -1; max_signal_er[index] = -1;}

  TCanvas *c1 = new TCanvas("c1", "Significance vs Iterations", 800, 600);
  TCanvas *c2 = new TCanvas("c2", "Search Space met 2D", 800, 600);
  TCanvas *c3 = new TCanvas("c3", "Search Space metoht 2D", 800, 600);
  TCanvas *c4 = new TCanvas("c4", "Search Space 3D", 800, 600);

  TH2F *hist2d_metoht = new TH2F("hist2d_metoht", "Search Space Histogram 2D", 30, bounds[0][0], bounds[0][1], 30, bounds[2][0], bounds[2][1]);
  TH2F *hist2d_met = new TH2F("hist2d_met", "Search Space Histogram 2D", 30, bounds[0][0], bounds[0][1], 30, bounds[1][0], bounds[1][1]);
  TH3F *hist3d = new TH3F("hist3d", "Search Space Histogram 3D", 20, bounds[0][0], bounds[0][1], 20, bounds[1][0], bounds[1][1], 20, bounds[2][0], bounds[2][1]);

  // Position uniform distributions
  uniform_real_distribution<> uni_dMetZPhi(bounds[0][0], bounds[0][1]);
  uniform_real_distribution<> uni_met_tst(bounds[1][0], bounds[1][1]);
  uniform_real_distribution<> uni_MetOHT(bounds[2][0], bounds[2][1]);
cout << " BOUND " << bounds[1][0] << endl << endl;
  // Velocity uniform distributions
  uniform_real_distribution<> uni_dMetZPhi2(-fabs(bounds[0][1] - bounds[0][0])/2, fabs(bounds[0][1] - bounds[0][0])/2); // bounds (-|min-max|, |min-max|)
  uniform_real_distribution<> uni_met_tst2(-fabs(bounds[1][1] - bounds[1][0])/2, fabs(bounds[1][1] - bounds[1][0])/2);
  uniform_real_distribution<> uni_MetOHT2(-fabs(bounds[2][1] - bounds[2][0])/2, fabs(bounds[2][1] - bounds[2][0])/2);

  random_device rd;
  mt19937 gen(rd());

  // Initialize particle swarm before first iterations
  vector<vector<particle>> swarm(iterations, vector<particle>(n_particles)); // Define the swarm

  for (int j = 0; j < n_particles; j++)
  {
    particle &particle = swarm[0][j]; // Define the particle
    // Initialize position
    particle.position[0] = uni_dMetZPhi(gen);
    particle.position[1] = uni_met_tst(gen);
    particle.position[2] = uni_MetOHT(gen);
    // Initialize pbest
    particle.pbest[0] = particle.position[0];
    particle.pbest[1] = particle.position[1];
    particle.pbest[2] = particle.position[2];
    // Initialize velocity
    particle.velocity[0] = uni_dMetZPhi2(gen);
    particle.velocity[1] = uni_met_tst2(gen);
    particle.velocity[2] = uni_MetOHT2(gen);
  }

  // Load root files
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

  
  for (int i = 0; i < iterations; i++)
  {
    for (int j = 0; j < n_particles; j++)
    {
      particle &particle = swarm[i][j]; // Define the particle

      cout << "  ____________________________________________________" << endl << endl;
      cout << "              ITERATION:  " << i + 1 << "   PARTICLE:  " << j + 1 << "   " << endl;
      cout << "  ____________________________________________________" << endl << endl << endl;
      
      cout << "    ================== DATA ==================    " << endl << endl;
      cout << "    DATA:";


      vector<Float_t> n_data = event_counter(particle, tree_data);


      //Break and try to prevent events from getting to 0
      if (n_data[0] == 0)
      {
        cout << "    No Events - Terminate " << endl;
        return;

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

      if (S != 0 && B != 0)
      {
        particle.significance = sqrt(2 * ((S + B) * log(1 + (S / B)) - S));
        particle.signal[0] = events_signal;
        particle.signal[1] = events_signal_er;
        
      }
      
      //Plot every particle update
      hist2d_met->Fill(particle.position[0], particle.position[1]);
      hist2d_metoht->Fill(particle.position[0], particle.position[2]);
      hist3d->Fill(particle.position[0], particle.position[1], particle.position[2]);

      if (particle.significance >= max_significance[j] && particle.significance != 0)
      { 
        max_significance[j] = particle.significance;
        max_signal[j] = particle.signal[0];
        max_signal_er[j] = particle.signal[1];
        

        is_new_pbest = is_new_pbest + 1;
  
        for (int index = 0; index < sizeof(particle.position) / sizeof(particle.position[0]); ++index)
        {
          particle.pbest[index] = particle.position[index];
        }
        // //Plot the particles with new pbest
        // hist2d_metst->Fill(particle.position[0], particle.position[1]);
        // hist2d_metoht->Fill(particle.pbest[0], particle.pbest[2]);
        // hist3d->Fill(particle.pbest[0], particle.pbest[1], particle.pbest[2]);


        if (particle.significance > gbest_significance)
        {
          gbest_significance = particle.significance;
          
          is_new_gbest = is_new_gbest + 1;

          for (int index = 0; index < sizeof(particle.pbest) / sizeof(particle.pbest[0]); ++index)
          {
            gbest[index] = particle.pbest[index];
          }
        }
      }
      else if (particle.significance < max_significance[j] && particle.significance != 0)
      {
        for (int index = 0; index < sizeof(particle.pbest) / sizeof(particle.pbest[0]); ++index)
        {
          particle.pbest[index] = swarm[i - 1][j].pbest[index];
        }
      }
      else 
      {
        cout << "   Unable to calculate particle significance  " << particle.significance <<  endl;
      }

      print_particle(particle, gbest, is_new_gbest, is_new_pbest, gbest_significance, max_significance, max_signal, max_signal_er, i, j);
    }

    // Update the swarm for the next iteration
    if (i < iterations - 1)
    {
      update_swarm(swarm, bounds, gbest, i, iterations, n_particles);
    }

    // Save the global best for plotting
    gbest_vector.push_back(gbest_significance);

    // Graph significance vs iterations


    c1->cd();
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
    graph->SetMarkerColor(kRed);
    graph->SetLineColor(kBlack);


    for (int i = 0; i < gbest_vector.size(); ++i)
    {
      graph->SetPoint(i, i + 1, gbest_vector[i]);
    }

    graph->Draw("ALP");

    Float_t current_best = gbest_vector[gbest_vector.size() - 1];
    Float_t x1 = 0.58, y1 = 0.12;
    Float_t x2 = 0.89, y2 = 0.28;

    char label[50];
    sprintf(label, "Global Best: (%.3f, %.3f, %.2f)", gbest[0], gbest[1], gbest[2]);

    // Adjust legend size to fit the parameter value
    TLegend *legend = new TLegend(x1, y1, x2, y2);
    legend->SetMargin(0.02);
    legend->SetTextSize(0.028);
    legend->AddEntry((TObject *)0, Form("Best Significance: %.3f", current_best), "");
    legend->AddEntry((TObject *)0, label, "");
    legend->Draw();


    c2->cd();
    c2->SetRightMargin(0.15);
    hist2d_met->GetXaxis()->SetTitle(" dMetZphi");
    hist2d_met->GetYaxis()->SetTitle(" met_tst");
    hist2d_met->GetZaxis()->SetTitle(" Entries");
    hist2d_met->SetTitle(" ");
    hist2d_met->SetStats(0);
    gStyle->SetPalette(kRainBow);
    // hist2d_met->SetContour(1000);
    hist2d_met->Draw("COLZ");
    // hist2d_met->Draw("LEGO2Z");

    c3->cd();
    c3->SetRightMargin(0.15);
    hist2d_metoht->GetXaxis()->SetTitle(" dMetZphi");
    hist2d_metoht->GetYaxis()->SetTitle(" MetOHT");
    hist2d_metoht->GetZaxis()->SetTitle(" Entries");
    hist2d_metoht->SetTitle(" ");
    hist2d_metoht->SetStats(0);
    gStyle->SetPalette(kRainBow);
    // hist2d_metoht->SetContour(1000);
    hist2d_metoht->Draw("COLZ");
    // hist2d_metoht->Draw("LEGO2Z");
    

    c4->cd();
    hist3d->GetXaxis()->SetTitle(" dMetZphi");
    hist3d->GetYaxis()->SetTitle(" met_tst");
    hist3d->GetZaxis()->SetTitle(" metOHT");
    hist3d->SetTitle(" ");
    hist3d->SetStats(0);
    // hist3d->SetContour(1000);
    // c4->SetRightMargin(0.15);
    gStyle->SetPalette(kRainBow);
    // hist3d->Draw("LEGO3Z");
    hist3d->Draw("BOX2 Z");


    c1->SetGrid();
    c1->Draw();
    c1->SaveAs("significance_vs_iterations.png");


    c2->Draw();
    if ((i+1) % 5 == 0)
    {
      string framename = "./3var/frames/met/met_iter_" + to_string(i+1) + ".png"; 
      c2->SaveAs(framename.c_str());
    }
    c2->SaveAs("search_space_met_2d.png");


    c3->Draw();
    if ((i+1) % 5 == 0)
    {
      string framename = "./3var/frames/metoht/metoht_iter_" + to_string(i+1) + ".png"; 
      c3->SaveAs(framename.c_str());
    }
    c3->SaveAs("search_space_metoht_2d.png");


    c4->Draw();
    c4->SaveAs("search_space_plot_3d.png");


  }



  cout << "   --------------------------------------------------------------------------" << endl << endl;
  cout << "   Max significance:      "  << gbest_significance << endl;
  cout << "   Best position:        (" << gbest [0] << ", " << gbest[1] << ", " << gbest[2] << ") " << endl;
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
