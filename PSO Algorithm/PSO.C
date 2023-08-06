//ROOT headers
#include <TFile.h>
#include <TROOT.h>
#include <TTree.h>
#include <TH1.h>


//Cpp headers
#include <iostream>
#include <vector>
#include <random>
#include <cmath>
using namespace std;


struct particle
{
  double dLepR;
  double dMetZPhi;
  double met_tst;
  double metOHT;
  double pbest;
  double significance;
  double max_signif;
  double position[4];
  double velocity[4];
};


void print_particle(const particle& p) {
    cout << "    dLepR: " << p.dLepR << endl;
    cout << "    dMetZPhi: " << p.dMetZPhi << endl;
    cout << "    met_tst: " << p.met_tst << endl;
    cout << "    metOHT: " << p.metOHT << endl;
    cout << "    significance: " << p.significance << endl;
    cout << "    max_signif: " << p.max_signif << endl;
    cout << "    position[0]: " << p.position[0] << endl;
    cout << "    position[1]: " << p.position[1] << endl;
    cout << "    position[2]: " << p.position[2] << endl;
    cout << "    position[3]: " << p.position[3] << endl;
    cout << "    velocity[0]: " << p.velocity[0] << endl;
    cout << "    velocity[1]: " << p.velocity[1] << endl;
    cout << "    velocity[2]: " << p.velocity[2] << endl;
    cout << "    velocity[3]: " << p.velocity[3] << endl;
    cout << endl;
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
    if (dLepR < particle.dLepR && dMetZPhi > particle.dMetZPhi && met_tst > particle.met_tst && MetOHT > particle.metOHT)
    {
      signal = signal + weight;
      signaler = signaler + weight * weight;
    }
  }

  cout << "    ENTRIES = " << tree->GetEntries() << endl << endl;
  cout << "          N = " << signal << "+-" << sqrt(signaler) << endl << endl;

  events.push_back(signal);
  events.push_back(sqrt(signaler));

  return events;
}

void PSO()
{

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

  int n_particles = 2;
  int iterations = 10;

  float_t w = 0.8;
  float_t c1 = 0.5;
  float_t c2 = 0.5;


  // Initialization
  uniform_real_distribution<> dist_uni(0., 1.);
  uniform_real_distribution<> dist_uni2(-1., 1.);
  random_device rd;
  mt19937 gen(rd());


  std::vector<particle> swarm(n_particles); //Define the swarm
  for (int i = 0; i < n_particles; i++)
  {
    
    particle particle;      //Define the particle 

    //Initialize position
    particle.dLepR = dist_uni(gen);
    particle.dMetZPhi = dist_uni(gen);
    particle.met_tst = dist_uni(gen);
    particle.metOHT = dist_uni(gen);
    particle.position[0] = particle.dLepR;
    particle.position[1] = particle.dMetZPhi;
    particle.position[2] = particle.met_tst;
    particle.position[3] = particle.metOHT;
    
    //Initialize velocity
    particle.dLepR = dist_uni2(gen);
    particle.dMetZPhi = dist_uni2(gen);
    particle.met_tst = dist_uni2(gen);
    particle.metOHT = dist_uni2(gen);
    particle.velocity[0] = particle.dLepR;
    particle.velocity[1] = particle.dMetZPhi;
    particle.velocity[2] = particle.met_tst;
    particle.velocity[3] = particle.metOHT;

    swarm.push_back(particle);

    


    cout << "PARTICLE:  "  << i << "   " << endl << endl;

    print_particle(particle);
    

    cout << "   ================== DATA ==================    " << endl << endl;
    cout << "   DATA:";
    vector<Float_t> n_data = Counter(particle, tree_data);

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
    Float_t Z = sqrt(2 * ((S + B) * log(1 + (S / B)) - S));
    

    cout << endl << "SIGNAL:  " << events_signal << endl << "Z:   " << Z << endl;

    
  }

  return;

}



