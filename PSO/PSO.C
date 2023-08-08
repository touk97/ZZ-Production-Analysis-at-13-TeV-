//ROOT headers
#include <TFile.h>
#include <TROOT.h>
#include <TTree.h>
#include <TH1.h>


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
  double dLepR;
  double dMetZPhi;
  double met_tst;
  double metOHT;
  double significance;
  double max_significance;
  double position[4];
  double pbest[4];
  double velocity[4];
};



void update_particle(vector<vector<particle>> &swarm, double *gbest, int i, int iterations, int n_particles)
{
  uniform_real_distribution<> uni_dist(0., 1.);
  random_device rd;
  mt19937 gen(rd());

  float_t w_min = 0.8;
  float_t w_max = 1.2;
  float_t w = w_max - (w_max - w_min) * i / iterations;
  float_t c1 = 0.5;
  float_t c2 = 0.5;
  float_t r1 = uni_dist(gen);
  float_t r2 = uni_dist(gen);

  for (int j = 0; j < n_particles; j++)
  {
    for (int k = 0; k < 4; k++)
    {
      swarm[i][j].velocity[k] = w * swarm[i - 1][j].velocity[k] + c1 * r1 * (swarm[i - 1][j].pbest[k] - swarm[i - 1][j].position[k]) + c2 * r2 * (gbest[k] - swarm[i - 1][j].position[k]);
      swarm[i][j].position[k] = swarm[i - 1][j].position[k] + swarm[i][j].velocity[k];
    }
  }
  return;
}



vector<Float_t> Counter(particle &particle, TTree *tree)
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
     
    if (dLepR < particle.position[0] && dMetZPhi > particle.position[1] && met_tst > particle.position[2] && MetOHT > particle.position[3])
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

void PSO()
{
  // Timer start
  auto start = std::chrono::high_resolution_clock::now();

  // Output log file
  ofstream logFile("./PSO.txt");

  DualStreamBuffer dualBuffer(std::cout.rdbuf(), logFile.rdbuf());

  std::streambuf *oldBuffer = std::cout.rdbuf(&dualBuffer);

  string filepath = "../../data/SAMPLES/SR/";
  cout << endl << endl << endl;
  cout << "   ------------------------------------------   " << endl;
  cout << "   FILEPATH:   " << filepath << endl;
  cout << "   ------------------------------------------   " << endl
       << endl;

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

  int iterations = 2;
  int n_particles = 2;

  double gbest_significance = - 1;
  double gbest[4];

  for (int i = 0; i < 4; i++) 
  {
    gbest[i] = -1;
  }



  // Initialization
  uniform_real_distribution<> uni_dist(0.1, 1.);
  uniform_real_distribution<> uni_dist2(-0.9, 1.1);

  random_device rd;
  mt19937 gen(rd());

  

  
  vector<vector<particle>> swarm(iterations, vector<particle>(n_particles)); //Define the swarm

  //Initialize the swarm
  for (int j = 0; j < n_particles; j++)
  {

    particle &particle = swarm[0][j]; // Define the particle

    // Initialize position
    particle.dLepR = uni_dist(gen) * (2.5 - 1.5);
    particle.dMetZPhi = uni_dist(gen) * (3. - 2.);
    particle.met_tst = uni_dist(gen) * (130. - 90.);
    particle.metOHT = uni_dist(gen) * (0.8 - 0.4);
    particle.position[0] = particle.dLepR;
    particle.position[1] = particle.dMetZPhi;
    particle.position[2] = particle.met_tst;
    particle.position[3] = particle.metOHT;

    // Initialize pbest
    particle.pbest[0] = particle.position[0];
    particle.pbest[1] = particle.position[1];
    particle.pbest[2] = particle.position[2];
    particle.pbest[3] = particle.position[3];

    // Initialize velocity
    particle.dLepR = uni_dist2(gen) * (2.5 - 1.5);
    particle.dMetZPhi = uni_dist2(gen) * (3. - 2.);
    particle.met_tst = uni_dist2(gen) * (130. - 90.);
    particle.metOHT = uni_dist2(gen) * (0.8 - 0.4);
    particle.velocity[0] = particle.dLepR;
    particle.velocity[1] = particle.dMetZPhi;
    particle.velocity[2] = particle.met_tst;
    particle.velocity[3] = particle.metOHT;

    particle.significance = -1.;
    particle.max_significance = -1.;

    // swarm.push_back(particle);
  }

  for (int i = 0; i < iterations; i++)
  {

    for (int j = 0; j < n_particles; j++)
    {
      update_particle(swarm, gbest, i, iterations, n_particles);

      particle &particle = swarm[i][j]; // Define the particle
      particle.position[0] = particle.dLepR;
      particle.position[1] = particle.dMetZPhi;
      particle.position[2] = particle.met_tst;
      particle.position[3] = particle.metOHT;

      particle.velocity[0] = particle.dLepR;
      particle.velocity[1] = particle.dMetZPhi;
      particle.velocity[2] = particle.met_tst;
      particle.velocity[3] = particle.metOHT;
    }

    cout << "   -------------------" << endl << endl;
    for (int j = 0; j < n_particles; j++)
    {
  
      particle &particle = swarm[i][j]; // Define the particle
  
      // cout << "   SWARM:  "  << swarm[i] << endl << endl;
      cout << "   ---------------------------------------------------" << endl;
      cout << "                  ITERATION:  " << i << "   PARTICLE:     "  << j << "   " << endl;
      cout << "   ---------------------------------------------------" << endl << endl << endl;
    
  
      cout << "   ================== DATA ==================    " << endl << endl;
      cout << "   DATA:";
      vector<Float_t> n_data = Counter(particle, tree_data);

      while (n_data[0] < 1)  //Re-initialize until the phase space produces non-zero events
      {
        // Initialize position
        particle.dLepR = uni_dist(gen) * (2.5 - 1.5);
        particle.dMetZPhi = uni_dist(gen) * (3. - 2.);
        particle.met_tst = uni_dist(gen) * (130. - 90.);
        particle.metOHT = uni_dist(gen) * (0.8 - 0.4);
        particle.position[0] = particle.dLepR;
        particle.position[1] = particle.dMetZPhi;
        particle.position[2] = particle.met_tst;
        particle.position[3] = particle.metOHT;

        //Initialize pbest
        particle.pbest[0] = particle.position[0] ;
        particle.pbest[1] = particle.position[1] ;
        particle.pbest[2] = particle.position[2] ;
        particle.pbest[3] = particle.position[3] ;

        // Initialize velocity
        particle.dLepR = uni_dist2(gen) * (2.5 - 1.5);
        particle.dMetZPhi = uni_dist2(gen) * (3. - 2.);
        particle.met_tst = uni_dist2(gen) * (130. - 90.);
        particle.metOHT = uni_dist2(gen) * (0.8 - 0.4);
        particle.velocity[0] = particle.dLepR;
        particle.velocity[1] = particle.dMetZPhi;
        particle.velocity[2] = particle.met_tst;
        particle.velocity[3] = particle.metOHT;
        n_data = Counter(particle, tree_data);
      }

      if (n_data[0] < 1)
      {
        cout << "--- no events ---" << endl << endl;
        cout << "-----------------" << endl << endl;
        continue; 
      }
    
      cout << "   ================== SIGNAL ==================    " << endl << endl;
      cout << "   llvv:";
      vector<Float_t> n_llvv = Counter(particle, tree_llvv);
      cout << "   llvvjj:";
      vector<Float_t> n_llvvjj = Counter(particle, tree_llvvjj);
  
      cout << "   ================== WZ ==================    " << endl << endl;
      cout << "   WZ:";
      vector<Float_t> n_WZ = Counter(particle, tree_WZ);
  
      cout << "   ================== Zjets ==================    " << endl << endl;
      cout << "   Z_jets_ee:";
      vector<Float_t> n_Zjets_ee = Counter(particle, tree_Z_jets_ee);
      cout << "   Z_jets_mumu:";
      vector<Float_t> n_Zjets_mumu = Counter(particle, tree_Z_jets_mumu);
  
      cout << "   ================== top ==================    " << endl << endl;
      cout << "   Top:";
      vector<Float_t> n_top = Counter(particle, tree_top);
      cout << "   ttbarV_ttbarVV:";
      vector<Float_t> n_ttbarV_ttbarVV = Counter(particle, tree_ttbarV_ttbarVV);
      cout << "   Wt:";
      vector<Float_t> n_Wt = Counter(particle, tree_Wt);
  
      cout << "   ================== WW ==================    " << endl << endl;
      cout << "   WW:";
      vector<Float_t> n_WW = Counter(particle, tree_WW);
  
      cout << "   ================== Othr ==================    " << endl << endl;
      cout << "   llll:";
      vector<Float_t> n_llll = Counter(particle, tree_llll);
      cout << "   llqq:";
      vector<Float_t> n_llqq = Counter(particle, tree_llqq);
      cout << "   VVV:";
      vector<Float_t> n_VVV = Counter(particle, tree_VVV);
      cout << "   W_jets:";
      vector<Float_t> n_Wjets = Counter(particle, tree_W_jets);
      cout << "   Ztt:";
      vector<Float_t> n_Ztt = Counter(particle, tree_Ztt);
      cout << "   WZ_jj:";
      vector<Float_t> n_WZjj = Counter(particle, tree_WZ_jj);
      cout << "   lllljj:";
      vector<Float_t> n_lllljj = Counter(particle, tree_lllljj);
      cout << "   llvvjj_WW:";
      vector<Float_t> n_llvvjj_WW = Counter(particle, tree_llvvjj_WW);
  
  
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
      
      if (particle.significance > particle.max_significance)
      {
        particle.max_significance = particle.significance;
        for (int index = 0; index < 4; ++index)
        {
          particle.pbest[index] = particle.position[index];
        }

        if (particle.max_significance > gbest_significance)
        {
          gbest_significance = particle.max_significance;
          cout << "   New search space max significance:   " << gbest_significance << endl; 
          for (int index = 0; index < 4; ++index)
          {
            gbest[index] = particle.pbest[index];
          }
          cout << "   New search space best position:      (" << gbest[0] << ", " << gbest[1] << ", " << gbest[2] << ", " << gbest[3] << ") " << endl << endl; 
        }
      }
  
      cout << endl << "   SIGNAL:  " << events_signal << "+-" << events_signal_er << endl;
      cout << "   Z:         " << swarm[i][j].significance << endl;
    }
  }

  cout << "----------------------------------------------------" << endl << endl;
  cout << "   Max significance is:  " << gbest_significance << endl;
  cout << "   Best position:        (" << gbest[0] << ", " << gbest[1] << ", " << gbest[2] << ", " << gbest[3] << ") " << endl;
  cout << "----------------------------------------------------" << endl << endl;


  // Timer stop
  auto end = chrono::high_resolution_clock::now();
  chrono::duration<float> duration = end - start;
  cout << endl << "   Script executed in: " << int(duration.count() / 60.0) << " minutes" << " and " << int((duration.count() / 60.0 - int(duration.count() / 60.0))*60) << " s" <<  endl << endl;


  // For the log file
  std::cout.rdbuf(oldBuffer); // Restore the original stream buffer for std::cout
  logFile.close(); // Close the log file

  return;
  
}
