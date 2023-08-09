
#include <TRandom3.h>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <TH1.h>
#include <TH2.h>
#include <vector>
#include <TGraphErrors.h>
#include <TMarker.h>
#include <TPaveText.h>
#include <TLatex.h>
#include "TProfile.h"
#include "TMath.h"
#include "TStyle.h"
#include <TLatex.h>
#include <TPad.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TLine.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <chrono> //Timer
using namespace std;



vector<Float_t> Counter(TTree *tree, TH1F *hist, string directory)
{

  TH1::SetDefaultSumw2(kTRUE);

  Int_t nentries = (Int_t)tree->GetEntries();

  Double_t M2Lep = 0.; 
  Double_t met_tst = 0.;
  // Double_t met_signif = 0.;
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

  // tree->SetBranchAddress("M2Lep", &M2Lep);
  tree->SetBranchAddress("met_tst", &met_tst);
  // tree->SetBranchAddress("met_signif", &met_signif);
  tree->SetBranchAddress("dMetZPhi", &dMetZPhi);
  tree->SetBranchAddress("MetOHT", &MetOHT);
  tree->SetBranchAddress("dLepR", &dLepR);
  tree->SetBranchAddress("M2Lep", &M2Lep);
  tree->SetBranchAddress("leading_pT_lepton", &leading_pT_lepton);
  tree->SetBranchAddress("subleading_pT_lepton", &subleading_pT_lepton);
  tree->SetBranchAddress("Z_pT", &Z_pT);
  tree->SetBranchAddress("n_jets", &n_jets);
  tree->SetBranchAddress("n_bjets", &n_bjets);
  tree->SetBranchAddress("detajj", &detajj);
  tree->SetBranchAddress("mjj", &mjj);
  tree->SetBranchAddress("leading_jet_pt", &leading_jet_pt);
  tree->SetBranchAddress("second_jet_pt", &second_jet_pt);
  tree->SetBranchAddress("event_3CR", &event_3CR);
  tree->SetBranchAddress("event_type", &event_type);
  tree->SetBranchAddress("global_weight", &weight);

  Double_t signal = 0.;
  Double_t signaler = 0.;
  // Float_t Norm = 1;


  

  // Loop over events
  for (Int_t i = 0; i < nentries; i++)
  {
    tree->GetEntry(i);
    if (directory == "SR")
    {

      if (event_3CR == 0 && (event_type == 0 || event_type == 1) &&
          leading_pT_lepton > 30 && subleading_pT_lepton > 20 && M2Lep > 80 && M2Lep < 100 && n_bjets < 1 &&
          dLepR < 1.8 && dMetZPhi > 2.7 && met_tst > 110 && MetOHT > 0.65)
      {
        // Inclusive
        signal = signal + weight;
        signaler = signaler + weight * weight;
        hist->Fill(Z_pT, weight);
      }
      // {
        // // Exclusive
        // if (n_jets > 1 && mjj > 0 && leading_jet_pt > 30 && second_jet_pt > 30)
        // // if (met_tst > 70 && n_jets > 1 && mjj > 100 && leading_jet_pt > 30 && second_jet_pt > 30 && detajj > 1 && MetOHT > 0.3 && dMetZPhi > 2.2 && n_bjets < 1 && dLepR < 2.2 )
        // {
        //   signal = signal + weight;
        //   signaler = signaler + weight * weight;
        //   hist->Fill(Z_pTt, weight);
        // }
      // }
    }
    else 
    {

      //Inclusive
      signal = signal + weight;
      signaler = signaler + weight * weight;
      hist->Fill(Z_pT, weight);

      // // Exclusive
      // if (n_jets > 1 && mjj > 100 && leading_jet_pt > 30 && second_jet_pt > 30)
      // // if (met_tst > 70 && n_jets > 1 && mjj > 100 && leading_jet_pt > 30 && second_jet_pt > 30 && detajj > 1 && MetOHT > 0.3 && dMetZPhi > 2.2 && n_bjets < 1 && dLepR < 2.2 )
      // {
      //   signal = signal + weight;
      //   signaler = signaler + weight * weight;
      //   hist->Fill(Z_pT, weight);
      // }
    }
  }

  
  cout << "    ENTRIES = " << tree->GetEntries() << endl << endl; 
  cout << "          N = " << signal << "+-" << sqrt(signaler) << endl << endl; 


  events.push_back(signal);
  events.push_back(sqrt(signaler));


  return events;


}
//
//
//
//
vector<Float_t> ZCounter(TTree *tree, TH1F *hist0, TH1F *hist1, TH1F *hist2, string directory)
{

  TH1::SetDefaultSumw2(kTRUE);

  Int_t nentries = (Int_t)tree->GetEntries();

  Double_t M2Lep; 
  Double_t met_tst;
  // Double_t met_signif;
  Double_t dMetZPhi;
  // Float_t frac_pT;
  Double_t MetOHT;
  Double_t dPhiJ100met;
  Double_t dLepR;
  Double_t n_bjets;
  Double_t leading_pT_lepton;
  Double_t subleading_pT_lepton;
  Double_t Z_pT;
  Double_t n_jets;
  Double_t detajj;
  Double_t mjj;
  Double_t leading_jet_pt;
  Double_t second_jet_pt;
  Double_t event_3CR;
  Double_t event_type;
  Double_t weight;
  vector<Float_t> events;
  events.clear();

  // tree->SetBranchAddress("M2Lep", &M2Lep);
  tree->SetBranchAddress("met_tst", &met_tst);
  // tree->SetBranchAddress("met_signif", &met_signif);
  tree->SetBranchAddress("dMetZPhi", &dMetZPhi);
  tree->SetBranchAddress("MetOHT", &MetOHT);
  tree->SetBranchAddress("dLepR", &dLepR);
  tree->SetBranchAddress("M2Lep", &M2Lep);
  tree->SetBranchAddress("leading_pT_lepton", &leading_pT_lepton);
  tree->SetBranchAddress("subleading_pT_lepton", &subleading_pT_lepton);
  tree->SetBranchAddress("Z_pT", &Z_pT);
  tree->SetBranchAddress("n_jets", &n_jets);
  tree->SetBranchAddress("n_bjets", &n_bjets);
  tree->SetBranchAddress("detajj", &detajj);
  tree->SetBranchAddress("mjj", &mjj);
  tree->SetBranchAddress("leading_jet_pt", &leading_jet_pt);
  tree->SetBranchAddress("second_jet_pt", &second_jet_pt);
  tree->SetBranchAddress("event_3CR", &event_3CR);
  tree->SetBranchAddress("event_type", &event_type);
  tree->SetBranchAddress("global_weight", &weight);

  Float_t signal0 = 0;
  Float_t signal1 = 0;
  Float_t signal2 = 0;
  Float_t signaler0 = 0;
  Float_t signaler1 = 0;
  Float_t signaler2 = 0;
  // Float_t Norm = 1;

  // Loop over events

  for (Int_t i = 0; i < nentries; i++)
  {
    tree->GetEntry(i);
    if (directory == "SR")
    {
      if (event_3CR == 0 && (event_type == 0 || event_type == 1) &&
          leading_pT_lepton > 30 && subleading_pT_lepton > 20 && M2Lep > 80 && M2Lep < 100 && n_bjets < 1 &&
          dLepR < 1.8 && dMetZPhi > 2.7 && met_tst > 110 && MetOHT > 0.65)
      {
        if (n_jets < 1)
        {
          signal0 = signal0 + weight;
          signaler0 = signaler0 + weight * weight;
          hist0->Fill(Z_pT, weight);
        }
        else if (n_jets > 0 && n_jets < 2)
        {
          signal1 = signal1 + weight;
          signaler1 = signaler1 + weight * weight;
          hist1->Fill(Z_pT, weight);
        }
        else if (n_jets > 1)
        {
          signal2 = signal2 + weight;
          signaler2 = signaler2 + weight * weight;
          hist2->Fill(Z_pT, weight);
        }
      }
    }
    else
    {
      // Inclusive
      if (n_jets < 1)
      {
        signal0 = signal0 + weight;
        signaler0 = signaler0 + weight * weight;
        hist0->Fill(Z_pT, weight);
      }
      else if (n_jets > 0 && n_jets < 2)
      {
        signal1 = signal1 + weight;
        signaler1 = signaler1 + weight * weight;
        hist1->Fill(Z_pT, weight);
      }
      else if (n_jets > 1)
      {
        signal2 = signal2 + weight;
        signaler2 = signaler2 + weight * weight;
        hist2->Fill(Z_pT, weight);
      }
      // // Exclusive
      // if (n_jets > 1 && mjj > 100 && leading_jet_pt > 30 && second_jet_pt > 30)
      // // if (met_tst > 70 && n_jets > 1 && mjj > 100 && leading_jet_pt > 30 && second_jet_pt > 30 && detajj > 1 && MetOHT > 0.3 && dMetZPhi > 2.2 && n_bjets < 1 && dLepR < 2.2 )
      // {
      //   signal = signal + weight;
      //   signaler = signaler + weight * weight;
      //   hist->Fill(Z_pT, weight);
      // }
    }
  }

  cout << "     ENTRIES = " << tree->GetEntries() << endl << endl; 
  cout << "          N0 = " << signal0 << "+-" << sqrt(signaler0) << endl << endl; 
  cout << "          N1 = " << signal1 << "+-" << sqrt(signaler1) << endl << endl; 
  cout << "          N2 = " << signal2 << "+-" << sqrt(signaler2) << endl << endl; 


  events.push_back(signal0);
  events.push_back(sqrt(signaler0));
  events.push_back(signal1);
  events.push_back(sqrt(signaler1));
  events.push_back(signal2);
  events.push_back(sqrt(signaler2));


  return events;


}
//
//
//
//
//
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
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
// MAIN
void zjets_splitted_unsc()
{

  // Timer start
  auto start = std::chrono::high_resolution_clock::now();

  //Output log file
  ofstream logFile("../cro/zjets_splitted_unsc/zjets_splitted_unsc.txt");

  DualStreamBuffer dualBuffer(std::cout.rdbuf(), logFile.rdbuf());

  std::streambuf *oldBuffer = std::cout.rdbuf(&dualBuffer);


  gROOT->SetBatch(kTRUE);


  //-------------------------ANALYSIS-------------------------//


  vector<string> directories = {"3lCR", "emCR_B", "emCR_A", "Zjets0", "Zjets1", "Zjets2", "SR"};
  vector<string> filenames = {"DATA", "WZ", "Z_jets_ee", "Z_jets_mumu",
                              "top", "ttbarV_ttbarVV", "Wt", "WW",
                              "llll", "llqq", "VVV", "W_jets", "Ztt",
                              "lllljj", "llvv", "llvvjj", "llvvjj_WW",
                              "WZ_jj"};


  Float_t xbins[22] = {70, 100, 130, 160, 190, 220, 250, 280, 310, 340, 370, 400, 430, 460, 490, 520, 550, 580, 610, 640, 670, 700}; 
  // Float_t xbins[36] = {70, 80, 90, 100, 110, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 280, 310, 340, 370, 400, 430, 460, 490, 520, 550, 580, 610, 640, 670, 700}; //pTZ
  // Float_t xbins[29] = {0, 30, 50, 60, 70, 80, 90, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 345, 375, 410, 450, 500, 580, 700, 800, 1000, 1200};  //sTZ
  // Float_t xbins[27] = {100, 200, 220, 240, 260, 280, 300, 320, 345, 375, 410, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 1000, 1100, 1200, 1300, 1400, 1500};  //STjj
  // Float_t xbins[21] = {100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1050, 1100};  //mjj
  // Float_t xbins[31] = {-8, -7.5, -7.0, -6.5, -6.0, -5.5, -5.0, -4.5, -4, -3.5, -3, -2.5,-2, -1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 7.0, 9.0};  //detajj
  


  for (string &directory : directories)
  {
    string filepath = "../data/SAMPLES/" + directory + "/";
    cout << endl << endl << endl;
    cout << "   ------------------------------------------   " << endl;
    cout << "   FILEPATH:   " << filepath << endl;
    cout << "   ------------------------------------------   " << endl << endl;

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

    TH1::SetDefaultSumw2(kTRUE);

    // Data
    TH1F *hist_data = new TH1F("hist_data", " ", sizeof(xbins) / sizeof(xbins[0]) - 1, xbins);

    // Signal
    TH1F *hist_signal = new TH1F("Hist_signal", "", sizeof(xbins) / sizeof(xbins[0]) - 1, xbins);
    TH1F *hist_llvv = new TH1F("Hist_llvv", " ", sizeof(xbins) / sizeof(xbins[0]) - 1, xbins);
    TH1F *hist_llvvjj = new TH1F("Hist_llvvjj", " ", sizeof(xbins) / sizeof(xbins[0]) - 1, xbins);

    // WZ
    TH1F *hist_WZ = new TH1F("Hist_WZ", " ", sizeof(xbins) / sizeof(xbins[0]) - 1, xbins);

    // Zjets = Z_jets_ee + Z_jets_mumu
    TH1F *hist_Zjets = new TH1F("hist_Zjets", " ", sizeof(xbins) / sizeof(xbins[0]) - 1, xbins);
    TH1F *hist_Zjets0 = new TH1F("hist_Zjets0", " ", sizeof(xbins) / sizeof(xbins[0]) - 1, xbins);
    TH1F *hist_Zjets1 = new TH1F("hist_Zjets1", " ", sizeof(xbins) / sizeof(xbins[0]) - 1, xbins);
    TH1F *hist_Zjets2 = new TH1F("hist_Zjets2", " ", sizeof(xbins) / sizeof(xbins[0]) - 1, xbins);
    TH1F *hist_Zjets_ee = new TH1F("hist_Zjets_ee", " ", sizeof(xbins) / sizeof(xbins[0]) - 1, xbins);
    TH1F *hist_Zjets_mumu = new TH1F("hist_Zjets_mumu", " ", sizeof(xbins) / sizeof(xbins[0]) - 1, xbins);
    TH1F *hist_Zjets_ee0 = new TH1F("hist_Zjets_ee0", " ", sizeof(xbins) / sizeof(xbins[0]) - 1, xbins);
    TH1F *hist_Zjets_ee1 = new TH1F("hist_Zjets_ee1", " ", sizeof(xbins) / sizeof(xbins[0]) - 1, xbins);
    TH1F *hist_Zjets_ee2 = new TH1F("hist_Zjets_ee2", " ", sizeof(xbins) / sizeof(xbins[0]) - 1, xbins);
    TH1F *hist_Zjets_mumu0 = new TH1F("hist_Zjets_mumu0", " ", sizeof(xbins) / sizeof(xbins[0]) - 1, xbins);
    TH1F *hist_Zjets_mumu1 = new TH1F("hist_Zjets_mumu1", " ", sizeof(xbins) / sizeof(xbins[0]) - 1, xbins);
    TH1F *hist_Zjets_mumu2 = new TH1F("hist_Zjets_mumu2", " ", sizeof(xbins) / sizeof(xbins[0]) - 1, xbins);

    // top = top + ttbarV_ttbar_VV + Wt
    TH1F *hist_top = new TH1F("Hist_top", " ", sizeof(xbins) / sizeof(xbins[0]) - 1, xbins);
    TH1F *hist_ttbarV_ttbarVV = new TH1F("Hist_ttbarV_ttbarVV", " ", sizeof(xbins) / sizeof(xbins[0]) - 1, xbins);
    TH1F *hist_Wt = new TH1F("Hist_Wt", " ", sizeof(xbins) / sizeof(xbins[0]) - 1, xbins);

    // WW
    TH1F *hist_WW = new TH1F("Hist_WW", " ", sizeof(xbins) / sizeof(xbins[0]) - 1, xbins);

    // othr = llll, llqq, VVV, Wjets, Ztt
    TH1F *hist_othr = new TH1F("Hist_othr", " ", sizeof(xbins) / sizeof(xbins[0]) - 1, xbins);
    TH1F *hist_llll = new TH1F("Hist_llll", " ", sizeof(xbins) / sizeof(xbins[0]) - 1, xbins);
    TH1F *hist_llqq = new TH1F("Hist_llqq", " ", sizeof(xbins) / sizeof(xbins[0]) - 1, xbins);
    TH1F *hist_VVV = new TH1F("Hist_VVV", " ", sizeof(xbins) / sizeof(xbins[0]) - 1, xbins);
    TH1F *hist_W_jets = new TH1F("Hist_W_jets", " ", sizeof(xbins) / sizeof(xbins[0]) - 1, xbins);
    TH1F *hist_Ztt = new TH1F("Hist_Ztt", " ", sizeof(xbins) / sizeof(xbins[0]) - 1, xbins);
    TH1F *hist_WZ_jj = new TH1F("Hist_WZ_jj", " ", sizeof(xbins) / sizeof(xbins[0]) - 1, xbins);
    TH1F *hist_lllljj = new TH1F("Hist_lllljj", " ", sizeof(xbins) / sizeof(xbins[0]) - 1, xbins);
    TH1F *hist_llvvjj_WW = new TH1F("Hist_llvvjj_WW", " ", sizeof(xbins) / sizeof(xbins[0]) - 1, xbins);

    TH1F *hist_signal_sig = new TH1F("hist_signal_sig ", " ", sizeof(xbins) / sizeof(xbins[0]) - 1, xbins);

 
 
    cout << "   ================== DATA ==================    " << endl << endl;
    cout << "   DATA:";
    vector<Float_t> n_data = Counter(tree_data, hist_data, directory);

    cout << "   ================== SIGNAL ==================    " << endl << endl;
    cout << "   llvv:";
    vector<Float_t> n_llvv = Counter(tree_llvv, hist_llvv, directory);
    cout << "   llvvjj:";
    vector<Float_t> n_llvvjj = Counter(tree_llvvjj, hist_llvvjj, directory);

    cout << "   ================== WZ ==================    " << endl << endl;
    cout << "   WZ:";
    vector<Float_t> n_WZ = Counter(tree_WZ, hist_WZ, directory);

    cout << "   ================== Zjets ==================    " << endl << endl;
    cout << "   Z_jets_ee:";
    vector<Float_t> n_Zjets_ee = ZCounter(tree_Z_jets_ee, hist_Zjets_ee0, hist_Zjets_ee1, hist_Zjets_ee2, directory);
    cout << "   Z_jets_mumu:";
    vector<Float_t> n_Zjets_mumu = ZCounter(tree_Z_jets_mumu, hist_Zjets_mumu0, hist_Zjets_mumu1, hist_Zjets_mumu2, directory);

    cout << "   ================== top ==================    " << endl << endl;
    cout << "   Top:";
    vector<Float_t> n_top = Counter(tree_top, hist_top, directory);
    cout << "   ttbarV_ttbarVV:";
    vector<Float_t> n_ttbarV_ttbarVV = Counter(tree_ttbarV_ttbarVV, hist_ttbarV_ttbarVV, directory);
    cout << "   Wt:";
    vector<Float_t> n_Wt = Counter(tree_Wt, hist_Wt, directory);

    cout << "   ================== WW ==================    " << endl << endl;
    cout << "   WW:";
    vector<Float_t> n_WW = Counter(tree_WW, hist_WW, directory);

    cout << "   ================== Othr ==================    " << endl << endl;
    cout << "   llll:";
    vector<Float_t> n_llll = Counter(tree_llll, hist_llll, directory);
    cout << "   llqq:";
    vector<Float_t> n_llqq = Counter(tree_llqq, hist_llqq, directory);
    cout << "   VVV:";
    vector<Float_t> n_VVV = Counter(tree_VVV, hist_VVV, directory);
    cout << "   W_jets:";
    vector<Float_t> n_Wjets = Counter(tree_W_jets, hist_W_jets, directory);
    cout << "   Ztt:";
    vector<Float_t> n_Ztt = Counter(tree_Ztt, hist_Ztt, directory);
    cout << "   WZ_jj:";
    vector<Float_t> n_WZjj = Counter(tree_WZ_jj, hist_WZ_jj, directory);
    cout << "   lllljj:";
    vector<Float_t> n_lllljj = Counter(tree_lllljj, hist_lllljj, directory);
    cout << "   llvvjj_WW:";
    vector<Float_t> n_llvvjj_WW = Counter(tree_llvvjj_WW, hist_llvvjj_WW, directory);




    //Plotting

    TCanvas *c1 = new TCanvas("c1", "pTZ_SR", 1400, 600, 700, 700);
    TCanvas *c2 = new TCanvas("c2", "pTZ_3lCR", 1400, 600, 700, 700);
    TCanvas *c3 = new TCanvas("c3", "pTZ_emuB", 1400, 600, 700, 700);
    TCanvas *c4 = new TCanvas("c4", "pTZ_emuA", 1400, 600, 700, 700);
    TCanvas *c5 = new TCanvas("c5", "pTZ_Zjets0", 1400, 600, 700, 700);
    TCanvas *c6 = new TCanvas("c6", "pTZ_Zjets1", 1400, 600, 700, 700);
    TCanvas *c7 = new TCanvas("c7", "pTZ_Zjets2", 1400, 600, 700, 700);

    if (directory == "SR")
    {
      c1->cd();
    }
    else if (directory == "3lCR")
    {
      c2->cd();
    }
    else if (directory == "emCR_B")
    {
      c3->cd();
    }
    else if (directory == "emCR_A")
    {
      c4->cd();
    }
    else if (directory == "Zjets0")
    {
      c5->cd();
    }
    else if (directory == "Zjets1")
    {
      c6->cd();
    }
    else if (directory == "Zjets2")
    {
      c7->cd();
    }

    TPad *pad1 = new TPad("pad1", "pad1", 0.01, 0.23, 1., 1.);
    TPad *pad2 = new TPad("pad2", "pad2", 0.01, 0.01, 1., 0.25);

    pad1->SetBorderSize(0);
    pad1->SetBottomMargin(0.035);
    pad2->SetBottomMargin(0.35);
    pad2->SetTopMargin(0.0);
    pad2->SetBorderSize(0);
    pad1->Draw();
    pad2->Draw();

    gStyle->SetOptTitle(kTRUE);
    gStyle->SetLegendFont(1);
    gStyle->SetPalette(20);

    // Data
    hist_data->SetMarkerStyle(20);
    hist_data->SetMarkerSize(1.);
    hist_data->SetMarkerColor(kBlack);
    hist_data->SetLineColor(kBlack);

    // SIGNAL = llvv + llvvjj
    hist_signal->SetFillColor(TColor::GetColor("#DFFF00")); // BRIGHT YELLOW
    hist_signal->SetLineColor(kBlack);
    hist_signal->SetLineWidth(1);
    hist_signal->Add(hist_llvv);
    hist_signal->Add(hist_llvvjj);

    // WZ
    hist_WZ->SetFillColor(TColor::GetColor("#FFBF00")); // MUSTARD
    hist_WZ->SetLineColor(kBlack);
    hist_WZ->SetLineWidth(1);

    // WW
    hist_WW->SetFillColor(TColor::GetColor("#40E0D0")); // TURQUAZ
    hist_WW->SetLineColor(kBlack);
    hist_WW->SetLineWidth(1);

    // top = top + ttbarV_ttbar_VV + Wt
    hist_top->SetFillColor(TColor::GetColor("#6495ED")); // LIGHT BLUE
    hist_top->SetLineColor(kBlack);
    hist_top->SetLineWidth(1);
    hist_top->Add(hist_ttbarV_ttbarVV);
    hist_top->Add(hist_Wt);

    // Zjets = Zjets_ee + Zjets_mumu
    hist_Zjets->SetFillColor(TColor::GetColor("#DE3163")); // DARK PINK
    hist_Zjets->SetLineColor(kBlack);
    hist_Zjets->SetLineWidth(1);

    hist_Zjets0->Add(hist_Zjets_ee0);
    hist_Zjets0->Add(hist_Zjets_mumu0);
    hist_Zjets1->Add(hist_Zjets_ee1);
    hist_Zjets1->Add(hist_Zjets_mumu1);
    hist_Zjets2->Add(hist_Zjets_ee2);
    hist_Zjets2->Add(hist_Zjets_mumu2);

    // Othr = llll + llqq + VVV + Wjets + Ztt
    hist_othr->SetFillColor(TColor::GetColor("#50C878")); // KHAKI
    hist_othr->SetLineColor(kBlack);
    hist_othr->SetLineWidth(1);
    hist_othr->Add(hist_llll);
    hist_othr->Add(hist_llqq);
    hist_othr->Add(hist_VVV);
    hist_othr->Add(hist_W_jets);
    hist_othr->Add(hist_Ztt);
    hist_othr->Add(hist_WZ_jj);
    hist_othr->Add(hist_lllljj);
    hist_othr->Add(hist_llvvjj_WW);


    Float_t events_total;
    Float_t events_nonWZ;
    Float_t events_nonWZ_er;
    Float_t events_nontop;
    Float_t events_nontop_er;
    Float_t events_nonWW;
    Float_t events_nonWW_er;
    Float_t events_nonZjets;
    Float_t events_nonZjets_er;
    Float_t events_nonZjets0;
    Float_t events_nonZjets0_er;
    Float_t events_nonZjets1;
    Float_t events_nonZjets1_er;
    Float_t events_nonZjets2;
    Float_t events_nonZjets2_er;


    Float_t events_data = n_data[0];
    Float_t events_data_er = n_data[1];
    Float_t events_signal = n_llvv[0] + n_llvvjj[0];
    Float_t events_signal_er = sqrt(pow(n_llvv[1], 2) + pow(n_llvvjj[1], 2));
    Float_t events_bkg;
    Float_t events_bkg_er;
    Float_t events_WZ = n_WZ[0];
    Float_t events_WZ_er = n_WZ[1];
    Float_t events_Zjets;
    Float_t events_Zjets_er;
    Float_t events_Zjets0 = n_Zjets_ee[0] + n_Zjets_mumu[0];
    Float_t events_Zjets0_er = sqrt(pow(n_Zjets_ee[1], 2) + pow(n_Zjets_mumu[1], 2));
    Float_t events_Zjets1 = n_Zjets_ee[2] + n_Zjets_mumu[2];
    Float_t events_Zjets1_er = sqrt(pow(n_Zjets_ee[3], 2) + pow(n_Zjets_mumu[3], 2));
    Float_t events_Zjets2 = n_Zjets_ee[4] + n_Zjets_mumu[4];
    Float_t events_Zjets2_er = sqrt(pow(n_Zjets_ee[5], 2) + pow(n_Zjets_mumu[5], 2));
    Float_t events_top = n_top[0] + n_ttbarV_ttbarVV[0] + n_Wt[0];
    Float_t events_top_er = sqrt(pow(n_top[1],2) + pow(n_ttbarV_ttbarVV[1], 2) + pow(n_Wt[1], 2));
    Float_t events_WW = n_WW[0];
    Float_t events_WW_er = n_WW[1];
    Float_t events_othr = n_llll[0] + n_llqq[0] + n_VVV[0] + n_Wjets[0] + n_Ztt[0] + n_WZjj[0] + n_lllljj[0] + n_llvvjj_WW[0];
    Float_t events_othr_er = sqrt(pow(n_llll[1], 2) + pow(n_llqq[1], 2) + pow(n_VVV[1], 2) + pow(n_Wjets[1], 2) + pow(n_Ztt[1], 2) + pow(n_WZjj[1], 2) + pow(n_lllljj[1], 2) + pow(n_llvvjj_WW[1], 2));


    //Merge Zjets before start plotting
    hist_Zjets->Add(hist_Zjets0);
    hist_Zjets->Add(hist_Zjets1);
    hist_Zjets->Add(hist_Zjets2);

    Double_t Z_max;
    Double_t Bin_max;
    if (directory == "SR")
    {
       // Stacking with a specific order
       hist_Zjets->Add(hist_othr);
       hist_top->Add(hist_Zjets);
       hist_WW->Add(hist_top);
       hist_WZ->Add(hist_WW);
       cout << "   SIGNAL/BKG = " << hist_signal->Integral(-5000, 5000) / hist_WZ->Integral(-5000, 5000) << endl << endl;

       

       for (int bin = 1; bin < hist_signal->GetSize(); ++bin)
       {
        Double_t B = hist_WZ->GetBinContent(bin);
        Double_t S = hist_signal->GetBinContent(bin);

        if (B > 0 && S > 0)
        {
          Double_t Z = sqrt(2 * ((S + B) * log(1 + (S / B)) - S));
          // cout << "B=" << B << ", S=" << S << ", Z_bin=" << Z_bin << endl;
          hist_signal_sig->SetBinContent(bin, Z);

          if (Z > Z_max)
          {
            Z_max = Z;
            Bin_max = hist_signal_sig->GetBinLowEdge(bin);
          }
          
         
        }
       }

      //Signal
      events_bkg = events_WZ + events_top + events_WW + events_Zjets0 + events_Zjets1 + events_Zjets2 + events_othr;
      events_bkg_er = sqrt(pow(events_WZ_er, 2) + pow(events_top_er, 2) + pow(events_WW_er, 2) + pow(events_Zjets0_er, 2) + pow(events_Zjets1_er, 2) + pow(events_Zjets2_er, 2) + pow(events_othr_er, 2));
      events_signal = events_data - events_bkg;
      events_signal_er = sqrt(pow(events_data_er, 2) + pow(events_bkg_er, 2));

      //Print calculated events for every region
      
       cout << "------------------------------------------------------------------" << endl << endl;
       cout << endl << endl  << "   SIGNAL   =  " << events_signal << " +- " << events_signal_er  << endl;
       cout << "   SIGNAL/BKG = " << events_signal/events_bkg << endl;
       cout << "   Maximum significance is:  " << Z_max << "   for bin:   " << Bin_max << endl << endl;
       cout << "_________________________________" << endl << endl;
       cout << "   Data: " << "     " << events_data << " +- " << events_data_er << endl << endl;
       cout << "   WZ: " << "       " << events_WZ << " +- " << events_WZ_er << endl << endl;
       cout << "   Top: " << "      " << events_top << " +- " << events_top_er << endl << endl;
       cout << "   WW: " << "       " << events_WW << " +- " << events_WW_er << endl << endl;
       cout << "   Zjets0: " << "   " << events_Zjets0 << " +- " << events_Zjets0_er << endl << endl;
       cout << "   Zjets1: " << "   " << events_Zjets1 << " +- " << events_Zjets1_er << endl << endl;
       cout << "   Zjets2: " << "   " << events_Zjets2 << " +- " << events_Zjets2_er << endl << endl;
       cout << "   Other: " << "    "  << events_othr << " +- " << events_othr_er << endl << endl;
       cout << "------------------------------------------------------------------" << endl << endl;
    }
    else
    {
       // Stacking with a specific order
       hist_Zjets->Add(hist_othr);
       hist_top->Add(hist_Zjets);
       hist_WW->Add(hist_top);
       hist_WZ->Add(hist_WW);
       hist_signal->Add(hist_WZ);
    }


    //  cout << "   DATA:     " << "MEAN =     " << hist_data->GetMean() << endl;
    //  cout << "             " << "RMS =      " << hist_data->GetRMS() << endl;
    //  cout << "             " << "INTEGRAL = " << hist_data->Integral(-3000, 3000) << endl << endl;
    //  cout << "   MC:       " << "MEAN =     " << hist_signal->GetMean() << endl;
    //  cout << "             " << "RMS =      " << hist_signal->GetRMS() << endl;
    //  cout << "             " << "INTEGRAL = " << hist_signal->Integral(-3000, 3000) << endl << endl;


    //----------------------------------------PLOTS----------------------------------------//

    // PAD 1

    // Plotting according to stacking order
    pad1->cd();

    // hist_signal->Draw("hist");
    // hist_WZ->Draw("histsame");
    // hist_WW->Draw("histsame");
    // hist_top->Draw("histsame");
    // hist_Zjets->Draw("histsame");
    // hist_othr->Draw("histsame");
    // hist_data->Draw("same");
    
    if (directory == "SR")
    {
      hist_WZ->Draw("hist");
      hist_signal->Draw("histsame");
      hist_WZ->GetYaxis()->SetRangeUser(0, hist_signal->GetMaximum() * 1.4);
      hist_WZ->SetStats(0);
      hist_WZ->SetLineWidth(2);
      hist_WZ->SetLineColor(kBlue);
      hist_WZ->SetFillColorAlpha(TColor::GetColor("#6495ED"), 0.9);
      hist_WZ->SetLineWidth(2);
      hist_WZ->GetYaxis()->SetTitle("Events");
      hist_signal->SetLineColor(TColor::GetColor("#DE3163"));
      hist_signal->SetFillColorAlpha(TColor::GetColor("#DE3163"), 0.8);
      hist_signal->SetLineWidth(2);
      hist_signal->SetFillStyle(3244);
    }
    else
    {
      hist_signal->Draw("hist");
      hist_WZ->Draw("histsame");
      hist_WW->Draw("histsame");
      hist_top->Draw("histsame");
      hist_Zjets->Draw("histsame");
      hist_othr->Draw("histsame");
      hist_data->Draw("same");
    }

    hist_signal->GetXaxis()->SetTitleSize(0.03);
    hist_signal->GetXaxis()->SetTitleFont(42);
    hist_signal->GetYaxis()->SetTitle("Events");
    hist_signal->GetXaxis()->SetTitleOffset(1.1);
    hist_signal->GetYaxis()->SetTitleFont(42);
    hist_signal->GetYaxis()->SetRangeUser(0, hist_signal->GetMaximum() * 1.4);
    hist_signal->SetStats(0);

    pad1->RedrawAxis();

    // TLatex *tex1 = new TLatex(0.3, 0.8, "#intL dt = 138.9 fb^{-1}");
    TLatex *tex1 = new TLatex(0.6, 0.55, "#intL dt = 138.9 fb^{-1}");
    tex1->SetNDC();
    tex1->SetTextFont(42);
    tex1->SetTextSize(0.04);
    tex1->SetLineWidth(2);
    tex1->Draw();
    // TLatex *tex2 = new TLatex(0.3, 0.7, "#sqrt{s} = 13 TeV");
    TLatex *tex2 = new TLatex(0.6, 0.45, "#sqrt{s} = 13 TeV");
    tex2->SetNDC();
    tex2->SetTextFont(42);
    tex2->SetTextSize(0.04);
    tex2->SetLineWidth(2);
    tex2->Draw();

    // TLegend *leg = new TLegend(0.6, 0.65, 0.9, 0.85);
    // leg->SetBorderSize(0);
    // leg->SetTextSize(0.03);
    // leg->SetMargin(0.25);
    // leg->SetEntrySeparation(0.1);
    // leg->SetNColumns(2);
    // leg->AddEntry(hist_data, "Data", "p");
    // leg->AddEntry(hist_signal, "Signal", "f");
    // leg->AddEntry(hist_WZ, "WZ", "f");
    // leg->AddEntry(hist_WW, "WW", "f");
    // leg->AddEntry(hist_Zjets, "Z+jets", "f");
    // leg->AddEntry(hist_top, "top", "f");
    // leg->AddEntry(hist_othr, "othr", "f");
    // leg->Draw("same");

    if (directory == "SR")
    {
      TLegend *leg = new TLegend(0.6, 0.65, 0.9, 0.85);
      leg->SetBorderSize(0);
      leg->SetTextSize(0.03);
      leg->SetMargin(0.25);
      leg->SetEntrySeparation(0.1);
      leg->AddEntry(hist_signal, "Signal", "f");
      leg->AddEntry(hist_WZ, "Backgournd", "f");
      leg->Draw("same");
    }
    else
    {
      TLegend *leg = new TLegend(0.6, 0.65, 0.9, 0.85);
      leg->SetBorderSize(0);
      leg->SetTextSize(0.03);
      leg->SetMargin(0.25);
      leg->SetEntrySeparation(0.1);
      leg->SetNColumns(2);
      leg->AddEntry(hist_data, "Data", "p");
      leg->AddEntry(hist_signal, "Signal", "f");
      leg->AddEntry(hist_WZ, "WZ", "f");
      leg->AddEntry(hist_WW, "WW", "f");
      leg->AddEntry(hist_Zjets, "Z+jets", "f");
      leg->AddEntry(hist_top, "top", "f");
      leg->AddEntry(hist_othr, "othr", "f");
      leg->Draw("same");
    }

    pad2->cd();
    

    TH1F *numerator = new TH1F("Numerator", " ", sizeof(xbins) / sizeof(xbins[0]) - 1, xbins);
    TH1F *denominator = new TH1F("Numerator", " ", sizeof(xbins) / sizeof(xbins[0]) - 1, xbins);

    if (directory == "SR")
    {
      hist_signal_sig->Draw("E0X0");
      hist_signal_sig->SetLineColor(kBlack);
      hist_signal_sig->SetMarkerStyle(20);
      hist_signal_sig->GetXaxis()->SetTitle("");
      hist_signal_sig->GetXaxis()->SetTitleSize(0.15);
      hist_signal_sig->GetXaxis()->SetTitleOffset(1.1);
      hist_signal_sig->GetXaxis()->SetLabelSize(0.14);

      hist_signal_sig->GetYaxis()->SetTitleOffset(0.4);
      hist_signal_sig->GetYaxis()->SetTitleSize(0.12);
      hist_signal_sig->GetYaxis()->SetLabelSize(0.10);
      hist_signal_sig->SetStats(0);
      pad2->SetGrid();
    }
    else if (directory == "Zjets0" || directory == "Zjets1" || directory == "Zjets2")
    {
      numerator->Add(hist_data);
      denominator->Add(hist_signal);
      numerator->GetYaxis()->SetRangeUser(-1, 4);
    }
    
    else
    {
      numerator->Add(hist_data);
      denominator->Add(hist_signal);
      // pad2->SetLogy();
      pad2->Update();
    }
    

    if (directory != "SR")
    {
      // numerator->Sumw2(1);
      // denominator->Sumw2(1);
      TLine *line = new TLine(xbins[0], 1, xbins[sizeof(xbins) / sizeof(xbins[0]) - 1], 1);
      line->SetLineStyle(2);
      line->SetLineWidth(2);
      line->Draw("same");
      pad2->RedrawAxis();

      numerator->SetMarkerStyle(20);
      numerator->GetXaxis()->SetTitle("");
      numerator->GetXaxis()->SetTitleSize(0.15);
      numerator->GetXaxis()->SetTitle("P^{Z}_{T} [GeV]");
      numerator->GetXaxis()->SetTitleOffset(1.1);
      numerator->GetXaxis()->SetLabelSize(0.14);

      numerator->GetYaxis()->SetTitleOffset(0.4);
      numerator->GetYaxis()->SetTitleSize(0.12);
      numerator->GetYaxis()->SetLabelSize(0.10);
      numerator->SetStats(0);

    // numerator->GetYaxis()->SetLabelSize(0.12);
    // numerator->GetYaxis()->SetTitle("Events");
    // numerator->GetYaxis()->SetTitleSize(0.05);

      numerator->SetTitle(" ");
      numerator->SetLineColor(kBlack);
      numerator->Divide(denominator);

      numerator->Draw();
    }

    TLine *line = new TLine(xbins[0], 1, xbins[sizeof(xbins) / sizeof(xbins[0]) - 1], 1);
    line->SetLineStyle(2);
    line->SetLineWidth(2);
    line->Draw("same");
    pad2->RedrawAxis();
    pad2->Update();

    if (directory == "SR")
    {
      pad1->cd();
      // TLatex *tex3 = new TLatex(0.28, 0.65, "Signal Region");
      TLatex *tex3 = new TLatex(0.6, 0.4, "Signal Region");
      tex3->SetNDC();
      tex3->SetTextFont(1);
      tex3->SetTextSize(0.04);
      tex3->SetLineWidth(1);
      tex3->Draw();
      hist_signal_sig->GetYaxis()->SetTitle("Significance");
      hist_signal_sig->GetXaxis()->SetTitle("P^{Z}_{T} [GeV]");
      pad1->Update();
      c1->Update();
      c1->SaveAs("../cro/zjets_splitted_unsc/ptz_SR_splitted_unsc.png");
    }
    else if (directory == "3lCR")
    {
      pad1->cd();
      // TLatex *tex3 = new TLatex(0.28, 0.65, "3l Control Region");
      TLatex *tex3 = new TLatex(0.6, 0.4, "3l Control Region");
      tex3->SetNDC();
      tex3->SetTextFont(1);
      tex3->SetTextSize(0.04);
      tex3->SetLineWidth(1);
      tex3->Draw();
      numerator->GetYaxis()->SetTitle("#frac{Data}{MC}");

      pad1->Update();
      c2->Update();
      c2->SaveAs("../cro/zjets_splitted_unsc/ptz_3lCR_splitted_unsc.png");
    }
    else if (directory == "emCR_B")
    {
      pad1->cd();
      // TLatex *tex3 = new TLatex(0.28, 0.65, "emuB Control Region");
      TLatex *tex3 = new TLatex(0.6, 0.4, "emuB Control Region");
      tex3->SetNDC();
      tex3->SetTextFont(1);
      tex3->SetTextSize(0.04);
      tex3->SetLineWidth(1);
      tex3->Draw();
      numerator->GetYaxis()->SetTitle("#frac{Data}{MC}");

      pad1->Update();
      c3->Update();
      c3->SaveAs("../cro/zjets_splitted_unsc/ptz_emCR_B_splitted_unsc.png");
    }
    else if (directory == "emCR_A")
    {
      pad1->cd();
      // TLatex *tex3 = new TLatex(0.58, 0.65, "emuA Control Region");
      TLatex *tex3 = new TLatex(0.6, 0.4, "emuA Control Region");
      tex3->SetNDC();
      tex3->SetTextFont(1);
      tex3->SetTextSize(0.04);
      tex3->SetLineWidth(1);
      tex3->Draw();
      numerator->GetYaxis()->SetTitle("#frac{Data}{MC}");

      pad1->Update();
      c4->Update();
      c4->SaveAs("../cro/zjets_splitted_unsc/ptz_emCR_A_splitted_unsc.png");
    }
    else if (directory == "Zjets0")
    {
      pad1->cd();
      // TLatex *tex3 = new TLatex(0.26, 0.65, "Zjets0 Control Region");
      TLatex *tex3 = new TLatex(0.6, 0.4, "Zjets0 Control Region");
      tex3->SetNDC();
      tex3->SetTextFont(1);
      tex3->SetTextSize(0.04);
      tex3->SetLineWidth(1);
      tex3->Draw();
      numerator->GetYaxis()->SetTitle("#frac{Data}{MC}");

      pad1->Update();
      c5->Update();
      c5->SaveAs("../cro/zjets_splitted_unsc/ptz_Zjets0_splitted_unsc.png");
    }
    else if (directory == "Zjets1")
    {
      pad1->cd();
      // TLatex *tex3 = new TLatex(0.26, 0.65, "Zjets1 Control Region");
      TLatex *tex3 = new TLatex(0.6, 0.4, "Zjets1 Control Region");
      tex3->SetNDC();
      tex3->SetTextFont(1);
      tex3->SetTextSize(0.04);
      tex3->SetLineWidth(1);
      tex3->Draw();
      numerator->GetYaxis()->SetTitle("#frac{Data}{MC}");

      pad1->Update();
      c6->Update();
      c6->SaveAs("../cro/zjets_splitted_unsc/ptz_Zjets1_splitted_unsc.png");
    }
    else if (directory == "Zjets2")
    {
      pad1->cd();
      // TLatex *tex3 = new TLatex(0.26, 0.65, "Zjets2 Control Region");
      TLatex *tex3 = new TLatex(0.6, 0.4, "Zjets2 Control Region");
      tex3->SetNDC();
      tex3->SetTextFont(1);
      tex3->SetTextSize(0.04);
      tex3->SetLineWidth(1);
      tex3->Draw();
      numerator->GetYaxis()->SetTitle("#frac{Data}{MC}");

      pad1->Update();
      c7->Update();
      c7->SaveAs("../cro/zjets_splitted_unsc/ptz_Zjets2_splitted_unsc.png");

    }


    //To avoid memory leak
    delete c1, c2, c3, c4, c5, c6, c7;
    delete file_data;
    delete file_llvv;
    delete file_llvvjj;
    delete file_WZ;
    delete file_Z_jets_ee;
    delete file_Z_jets_mumu;
    delete file_top;
    delete file_ttbarV_ttbarVV;
    delete file_Wt;
    delete file_WW;
    delete file_llll;
    delete file_llqq;
    delete file_VVV;
    delete file_W_jets;
    delete file_Ztt;
    delete file_WZ_jj;
    delete file_lllljj;
    delete file_llvvjj_WW;

  }
  // Timer stop
    auto end = chrono::high_resolution_clock::now();

    chrono::duration<float> duration = end - start;

    cout << endl << "   Script executed in: " << int(duration.count() / 60.0) << " minutes" << " and " << int((duration.count() / 60.0 - int(duration.count() / 60.0))*60) << " s" <<  endl << endl;

    // For the log file
    std::cout.rdbuf(oldBuffer); // Restore the original stream buffer for std::cout

    logFile.close(); // Close the log file

    return;
}
