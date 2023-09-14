#include <TRandom3.h>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <TH1.h>
#include <TH2.h>
#include "TProfile.h"
#include "TMath.h"
#include "TStyle.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLegendEntry.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
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

// gStyle->SetLegendFont(42);
// gStyle->SetPalette(1);

std::vector<float> Counter(TTree *tree, TH1F *hist, Float_t met_tst_value, Float_t dLepR_value, Float_t dMetZPhi_value, Float_t MetOHT_value)
{

  TH1::SetDefaultSumw2(kTRUE);

  Int_t nentries = (Int_t)tree->GetEntries();

  Double_t M2Lep = 0.;
  Double_t M2Lep_signif;
  Double_t met_tst = 0.;
  Double_t met_signif = 0.;
  Double_t dMetZPhi = 0.;
  Double_t frac_pT = 0.;
  Double_t MetOHT = 0.;
  Double_t dPhiJ100met = 0.;
  Double_t dLepR = 0.;
  Double_t n_bjets = 0.;
  Double_t n_bjets_signif = 0.;
  Double_t n_jets = 0.;
  Double_t leading_pT_lepton = 0;
  Double_t subleading_pT_lepton = 0;
  Double_t detajj = 0;
  Double_t mjj = 0;
  Double_t leading_jet_pt = 0;
  Double_t second_jet_pt = 0;
  Double_t event_3CR = 0.;
  Double_t event_type = 0.;
  Double_t weight = 1.;

  std::vector<float> events;
  events.clear();

  vector<TString> branches = {
      // "M2Lep",
      "met_tst",
      // "met_signif",
      "dMetZPhi",
      "MetOHT",
      "dLepR",
      // "M2Lep",
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

  for (const auto &branch : branches)
  {
    tree->SetBranchStatus(branch, 1);
  }

  tree->SetBranchAddress("M2Lep", &M2Lep);
  tree->SetBranchAddress("met_tst", &met_tst);
  tree->SetBranchAddress("met_signif", &met_signif);
  tree->SetBranchAddress("dMetZPhi", &dMetZPhi);
  tree->SetBranchAddress("MetOHT", &MetOHT);
  tree->SetBranchAddress("dLepR", &dLepR);
  tree->SetBranchAddress("leading_pT_lepton", &leading_pT_lepton);
  tree->SetBranchAddress("subleading_pT_lepton", &subleading_pT_lepton);
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

  // Loop over events
  for (int i = 0; i < nentries; i++)
  {
    tree->GetEntry(i);

    //SR SAMPLE -> FID.VOLUME: 80 < M2Lep < 100, met_tst > 70, dLepR < 1.8, dMetZPhi > 2.2, nbjets < 1, (event_3CR == 0 && (event_type == 0 || event_type == 1)
    if (dLepR < dLepR_value && dMetZPhi > dMetZPhi_value && met_tst > met_tst_value && MetOHT > MetOHT_value)

    // if (event_3CR == 0 && (event_type == 0 || event_type == 1) &&
    //       leading_pT_lepton > 30 && subleading_pT_lepton > 20  && M2Lep > 80 && M2Lep < 100 && n_bjets < 1 &&
    //       dLepR < 1.8 && dMetZPhi > 2.7 && met_tst > 110 && MetOHT > 0.65)

    {
      signal = signal + weight;              // signal yield is sum of weights
      signaler = signaler + weight * weight; // keep Sum(weights^2) for calculating the error on the signal yield
      hist->Fill(met_tst, weight);
      // Fills histogram but with wieght
    }
  }

  cout << "   N = " << signal << " +- " << sqrt(signaler) << endl << endl; // signal yield and error on this yield

  events.push_back(signal);
  events.push_back(sqrt(signaler));

  return events;
}
//
//
//
//
//
//
//
//
int bro_plot_SR()
{

  // Timer start
  auto start = std::chrono::high_resolution_clock::now();

  // Output log file
  ofstream logFile("./bro_plot.txt");

  DualStreamBuffer dualBuffer(std::cout.rdbuf(), logFile.rdbuf());

  std::streambuf *oldBuffer = std::cout.rdbuf(&dualBuffer);

  gStyle->SetLegendFont(42);
  gStyle->SetPalette(1);
  Float_t dLepR_max = -1.0;
  Float_t dMetZphi_max = -1.0;
  Float_t met_tst_max = -1.0;
  Float_t MetOHT_max = -1.0;
  Float_t Z_max = -1.0;
  Float_t iteration = 1;

  //SR SAMPLE -> FID.VOLUME: 80 < M2Lep < 100, met_tst > 70, dLepR < 1.8, dMetZPhi > 2.2, nbjets < 1, 
  //(event_3CR == 0 && (event_type == 0 || event_type == 1)

  Float_t dLepR_value[5] = {1.6, 1.65, 1.7, 1.75, 1.8};
  Float_t dMetZphi_value[11] = {2.2, 2.25, 2.3, 2.35, 2.4, 2.45, 2.5, 2.55, 2.6, 2.7, 2.8};
  Float_t met_tst_value[1] = {90};
  Float_t MetOHT_value[5] = {0.6, 0.65, 0.7, 0.75, 0.8};

  // Float_t dLepR_value[4] = {1.8, 2.0, 2.2, 2.5};
  // Float_t dMetZphi_value[4] = {2.2, 2.5, 2.7, 3.0};
  // Float_t met_tst_value[1] = {90};
  // Float_t MetOHT_value[4] = {0.5, 0.6, 0.65, 0.7};

  // Float_t dLepR_value[4] = {1.6, 1.7, 1.8, 1.9};
  // Float_t dMetZphi_value[5] = {2.0, 2.1, 2.2, 2.3, 2.4};
  // Float_t met_tst_value[1] = {90};
  // Float_t MetOHT_value[5] = {0.6, 0.65, 0.7, 0.75, 0.8};

  // Float_t dLepR_value[2] = {1.8, 2.0};
  // Float_t dMetZphi_value[2] = {2.2, 2.5};
  // Float_t met_tst_value[1] = {90};
  // Float_t MetOHT_value[2] = {0.5, 0.6};

  // Float_t dLepR_value[10] = {1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5};
  // Float_t dMetZphi_value[9] = {2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0};
  // Float_t met_tst_value[1] = {80};
  // Float_t MetOHT_value[7] = {0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8};

  // Float_t dLepR_value[9] = {1.62, 1.64, 1.66, 1.68, 1.7, 1.72, 1.74, 1.76, 1.78};
  // Float_t dMetZphi_value[9] = {2.62, 2.64, 2.66, 2.68, 2.7, 2.72, 2.74, 2.76, 2.78};
  // Float_t met_tst_value[1] = {80};
  // Float_t MetOHT_value[9] = {0.66, 0.67, 0.68, 0.69, 0.7, 0.71, 0.72, 0.73, 0.74};

  string filepath = "/home/touk/Desktop/touk/master/thesis/data/SAMPLES/SR/";

  // Float_t xbins[24] = {80, 90, 100, 110, 120, 130, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 345, 375, 410, 450, 500, 580, 700, 1000};

  Float_t xbins[29] = {80, 85, 90, 95, 100, 105, 110, 115, 120, 125, 130, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 345, 375, 410, 450, 500, 580, 700, 1000};

  Float_t binmin = 80;
  Float_t binmax = 500;

  // Data File
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

  for (Int_t i = 0; i < sizeof(dLepR_value) / sizeof(dLepR_value[0]); i++)
  {
    for (Int_t j = 0; j < sizeof(dMetZphi_value) / sizeof(dMetZphi_value[0]); j++)
    {
      for (Int_t k = 0; k < sizeof(MetOHT_value) / sizeof(MetOHT_value[0]); k++)
      {
        for (Int_t l = 0; l < sizeof(met_tst_value) / sizeof(met_tst_value[0]); l++)
        {

          cout << endl
               << endl
               << "   --------------------------------------   " << endl;
          cout << "              ITERATION #:  " << iteration << endl;
          cout << "   --------------------------------------   " << endl
               << endl;

          // Data
          TH1F *hist_data = new TH1F("hist_data", " ", sizeof(xbins) / sizeof(xbins[0]) - 1, xbins);

          // Signal
          TH1F *hist_signal = new TH1F("hist_signal", "", sizeof(xbins) / sizeof(xbins[0]) - 1, xbins);
          TH1F *hist_llvv = new TH1F("hist_llvv", " ", sizeof(xbins) / sizeof(xbins[0]) - 1, xbins);
          TH1F *hist_llvvjj = new TH1F("hist_llvvjj", " ", sizeof(xbins) / sizeof(xbins[0]) - 1, xbins);

          // WZ
          TH1F *hist_WZ = new TH1F("hist_WZ", " ", sizeof(xbins) / sizeof(xbins[0]) - 1, xbins);

          // Zjets = Z_jets_ee + Z_jets_mumu
          TH1F *hist_Zjets = new TH1F("hist_Zjets", " ", sizeof(xbins) / sizeof(xbins[0]) - 1, xbins);
          TH1F *hist_Zjets_ee = new TH1F("hist_Zjets_ee", " ", sizeof(xbins) / sizeof(xbins[0]) - 1, xbins);
          TH1F *hist_Zjets_mumu = new TH1F("hist_Zjets_mumu", " ", sizeof(xbins) / sizeof(xbins[0]) - 1, xbins);

          // top = top + ttbarV_ttbar_VV + Wt
          TH1F *hist_top = new TH1F("hist_top", " ", sizeof(xbins) / sizeof(xbins[0]) - 1, xbins);
          TH1F *hist_ttbarV_ttbarVV = new TH1F("hist_ttbarV_ttbarVV", " ", sizeof(xbins) / sizeof(xbins[0]) - 1, xbins);
          TH1F *hist_Wt = new TH1F("hist_Wt", " ", sizeof(xbins) / sizeof(xbins[0]) - 1, xbins);

          // WW
          TH1F *hist_WW = new TH1F("hist_WW", " ", sizeof(xbins) / sizeof(xbins[0]) - 1, xbins);

          // othr = llll, llqq, VVV, Wjets, Ztt
          TH1F *hist_othr = new TH1F("hist_othr", " ", sizeof(xbins) / sizeof(xbins[0]) - 1, xbins);
          TH1F *hist_llll = new TH1F("hist_llll", " ", sizeof(xbins) / sizeof(xbins[0]) - 1, xbins);
          TH1F *hist_llqq = new TH1F("hist_llqq", " ", sizeof(xbins) / sizeof(xbins[0]) - 1, xbins);
          TH1F *hist_VVV = new TH1F("hist_VVV", " ", sizeof(xbins) / sizeof(xbins[0]) - 1, xbins);
          TH1F *hist_W_jets = new TH1F("hist_W_jets", " ", sizeof(xbins) / sizeof(xbins[0]) - 1, xbins);
          TH1F *hist_Ztt = new TH1F("hist_Ztt", " ", sizeof(xbins) / sizeof(xbins[0]) - 1, xbins);
          TH1F *hist_WZ_jj = new TH1F("hist_WZ_jj", " ", sizeof(xbins) / sizeof(xbins[0]) - 1, xbins);
          TH1F *hist_lllljj = new TH1F("hist_lllljj", " ", sizeof(xbins) / sizeof(xbins[0]) - 1, xbins);
          TH1F *hist_llvvjj_WW = new TH1F("hist_llvvjj_WW", " ", sizeof(xbins) / sizeof(xbins[0]) - 1, xbins);

          TH1F *hist_signal_sig = new TH1F("hist_signal_sig ", " ", sizeof(xbins) / sizeof(xbins[0]) - 1, xbins);

          // Event Yields
          cout << "   ================== DATA ==================    " << endl
               << endl;
          cout << "   DATA:";
          vector<Float_t> n_data = Counter(tree_data, hist_data, met_tst_value[l], dLepR_value[i], dMetZphi_value[j], MetOHT_value[k]);

          cout << "   ================== SIGNAL ==================    " << endl
               << endl;
          cout << "   llvv:";
          vector<Float_t> n_llvv = Counter(tree_llvv, hist_llvv, met_tst_value[l], dLepR_value[i], dMetZphi_value[j], MetOHT_value[k]);
          cout << "   llvvjj:";
          vector<Float_t> n_llvvjj = Counter(tree_llvvjj, hist_llvvjj, met_tst_value[l], dLepR_value[i], dMetZphi_value[j], MetOHT_value[k]);

          cout << "   ================== WZ ==================    " << endl
               << endl;
          cout << "   WZ:";
          vector<Float_t> n_WZ = Counter(tree_WZ, hist_WZ, met_tst_value[l], dLepR_value[i], dMetZphi_value[j], MetOHT_value[k]);

          cout << "   ================== Zjets ==================    " << endl
               << endl;
          cout << "   Z_jets_ee:";
          vector<Float_t> n_Zjets_ee = Counter(tree_Z_jets_ee, hist_Zjets_ee, met_tst_value[l], dLepR_value[i], dMetZphi_value[j], MetOHT_value[k]);
          cout << "   Z_jets_mumu:";
          vector<Float_t> n_Zjets_mumu = Counter(tree_Z_jets_mumu, hist_Zjets_mumu, met_tst_value[l], dLepR_value[i], dMetZphi_value[j], MetOHT_value[k]);

          cout << "   ================== top ==================    " << endl
               << endl;
          cout << "   Top:";
          vector<Float_t> n_top = Counter(tree_top, hist_top, met_tst_value[l], dLepR_value[i], dMetZphi_value[j], MetOHT_value[k]);
          cout << "   ttbarV_ttbarVV:";
          vector<Float_t> n_ttbarV_ttbarVV = Counter(tree_ttbarV_ttbarVV, hist_ttbarV_ttbarVV, met_tst_value[l], dLepR_value[i], dMetZphi_value[j], MetOHT_value[k]);
          cout << "   Wt:";
          vector<Float_t> n_Wt = Counter(tree_Wt, hist_Wt, met_tst_value[l], dLepR_value[i], dMetZphi_value[j], MetOHT_value[k]);

          cout << "   ================== WW ==================    " << endl
               << endl;
          cout << "   WW:";
          vector<Float_t> n_WW = Counter(tree_WW, hist_WW, met_tst_value[l], dLepR_value[i], dMetZphi_value[j], MetOHT_value[k]);

          cout << "   ================== Othr ==================    " << endl
               << endl;
          cout << "   llll:";
          vector<Float_t> n_llll = Counter(tree_llll, hist_llll, met_tst_value[l], dLepR_value[i], dMetZphi_value[j], MetOHT_value[k]);
          cout << "   llqq:";
          vector<Float_t> n_llqq = Counter(tree_llqq, hist_llqq, met_tst_value[l], dLepR_value[i], dMetZphi_value[j], MetOHT_value[k]);
          cout << "   VVV:";
          vector<Float_t> n_VVV = Counter(tree_VVV, hist_VVV, met_tst_value[l], dLepR_value[i], dMetZphi_value[j], MetOHT_value[k]);
          cout << "   W_jets:";
          vector<Float_t> n_Wjets = Counter(tree_W_jets, hist_W_jets, met_tst_value[l], dLepR_value[i], dMetZphi_value[j], MetOHT_value[k]);
          cout << "   Ztt:";
          vector<Float_t> n_Ztt = Counter(tree_Ztt, hist_Ztt, met_tst_value[l], dLepR_value[i], dMetZphi_value[j], MetOHT_value[k]);
          cout << "   WZ_jj:";
          vector<Float_t> n_WZjj = Counter(tree_WZ_jj, hist_WZ_jj, met_tst_value[l], dLepR_value[i], dMetZphi_value[j], MetOHT_value[k]);
          cout << "   lllljj:";
          vector<Float_t> n_lllljj = Counter(tree_lllljj, hist_lllljj, met_tst_value[l], dLepR_value[i], dMetZphi_value[j], MetOHT_value[k]);
          cout << "   llvvjj_WW:";
          vector<Float_t> n_llvvjj_WW = Counter(tree_llvvjj_WW, hist_llvvjj_WW, met_tst_value[l], dLepR_value[i], dMetZphi_value[j], MetOHT_value[k]);

          Float_t events_data = n_data[0];
          Float_t events_data_er = n_data[1];
          Float_t signal_mc = n_llvv[0] + n_llvvjj[0];
          Float_t signal_mc_er = sqrt(pow(n_llvv[1], 2) + pow(n_llvvjj[1], 2));
          Float_t events_WZ = n_WZ[0];
          Float_t events_WZ_er = n_WZ[1];
          Float_t events_Zjets = n_Zjets_ee[0] + n_Zjets_mumu[0];
          Float_t events_Zjets_er = sqrt(pow(n_Zjets_ee[1], 2) + pow(n_Zjets_mumu[1], 2));
          Float_t events_top = n_top[0] + n_ttbarV_ttbarVV[0] + n_Wt[0];
          Float_t events_top_er = sqrt(pow(n_top[1], 2) + pow(n_ttbarV_ttbarVV[1], 2) + pow(n_Wt[1], 2));
          Float_t events_WW = n_WW[0];
          Float_t events_WW_er = n_WW[1];
          Float_t events_othr = n_llll[0] + n_llqq[0] + n_VVV[0] + n_Wjets[0] + n_Ztt[0] + n_WZjj[0] + n_lllljj[0] + n_llvvjj_WW[0];
          Float_t events_othr_er = sqrt(pow(n_llll[1], 2) + pow(n_llqq[1], 2) + pow(n_VVV[1], 2) + pow(n_Wjets[1], 2) + pow(n_Ztt[1], 2) + pow(n_WZjj[1], 2) + pow(n_lllljj[1], 2) + pow(n_llvvjj_WW[1], 2));
          Float_t events_bkg = events_WZ + events_top + events_WW + events_Zjets + events_othr;
          Float_t events_bkg_er = sqrt(pow(events_WZ_er, 2) + pow(events_top_er, 2) + pow(events_WW_er, 2) + pow(events_Zjets_er, 2) + pow(events_othr_er, 2));

          cout << "   Bkg = " << events_bkg << " +- " << events_bkg_er << endl;
          cout << "   S/B = " << signal_mc / events_bkg << endl;
          Double_t S = signal_mc;
          Double_t B = events_bkg;
          Double_t Z = sqrt(2 * ((S + B) * log(1 + (S / B)) - S));
          cout << "   Significance = " << Z << endl << endl;

          // Merge signal contributions
          hist_signal->Add(hist_llvv);
          hist_signal->Add(hist_llvvjj);

          // top = top + ttbarV_ttbar_VV + Wt
          hist_top->Add(hist_ttbarV_ttbarVV);
          hist_top->Add(hist_Wt);

          // Zjets = Zjets_ee + Zjets_mumu
          hist_Zjets->Add(hist_Zjets_ee);
          hist_Zjets->Add(hist_Zjets_mumu);

          // Othr = llll + llqq + VVV + Wjets + Ztt
          hist_othr->Add(hist_llll);
          hist_othr->Add(hist_llqq);
          hist_othr->Add(hist_VVV);
          hist_othr->Add(hist_W_jets);
          hist_othr->Add(hist_Ztt);
          hist_othr->Add(hist_WZ_jj);
          hist_othr->Add(hist_lllljj);
          hist_othr->Add(hist_llvvjj_WW);

          hist_Zjets->Add(hist_othr);
          hist_top->Add(hist_Zjets);
          hist_WW->Add(hist_top);
          hist_WZ->Add(hist_WW);

          Double_t Z_bin;

          for (int bin = 1; bin < hist_WZ->GetSize(); ++bin)
          {
            Double_t B = hist_WZ->Integral(bin, hist_WZ->GetNbinsX());
            Double_t S = hist_signal->Integral(bin, hist_signal->GetNbinsX());

            if (B > 0 && S > 0)
            {
              Z_bin = sqrt(2 * ((S + B) * log(1 + (S / B)) - S));

              // cout << hist_WZ->GetBinLowEdge(bin) << "   " << Z_bin << "   " << Z << endl;
              hist_signal_sig->SetBinContent(bin, Z_bin);

              if (Z_bin > Z_max)
              {
                Z_max = Z_bin;
                dLepR_max = dLepR_value[i];
                dMetZphi_max = dMetZphi_value[j];
                MetOHT_max = MetOHT_value[k];
                met_tst_max = hist_signal_sig->GetBinLowEdge(bin);
                cout << "   New maximum significance = " << Z_max << endl << endl;
                cout << "   New optimal set of values are:  dLepR = " << dLepR_max << ",  dMetZphi = " << dMetZphi_max << ",  met_tst = " << met_tst_max << ",  MetOHT = " << MetOHT_max << endl << endl;
              }
            }
          }

          iteration += 1;

          delete hist_signal;
          delete hist_data;
          delete hist_llvv;
          delete hist_llvvjj;
          delete hist_Zjets;
          delete hist_Zjets_ee;
          delete hist_Zjets_mumu;
          delete hist_WZ;
          delete hist_WW;
          delete hist_top;
          delete hist_ttbarV_ttbarVV;
          delete hist_Wt;
          delete hist_othr;
          delete hist_llll;
          delete hist_llqq;
          delete hist_VVV;
          delete hist_W_jets;
          delete hist_Ztt;
          delete hist_WZ_jj;
          delete hist_lllljj;
          delete hist_llvvjj_WW;
          delete hist_signal_sig;

        }
      }
    }
  }

  cout << endl << "   Optimal Set:  dLepR = " << dLepR_max << ",  dMetZphi = " << dMetZphi_max << ",  met_tst = " << met_tst_max << ",  MetOHT = " << MetOHT_max << endl << endl;
  cout << endl << "   Maximum Significance:   " << Z_max << endl << endl;

  // Timer stop
  auto end = chrono::high_resolution_clock::now();
  chrono::duration<float> duration = end - start;
  cout << endl << "   Script executed in: " << int(duration.count() / 60.0) << " minutes" << " and " << int((duration.count() / 60.0 - int(duration.count() / 60.0)) * 60) << " s" << endl << endl;

  // For the log file
  std::cout.rdbuf(oldBuffer); // Restore the original stream buffer for std::cout
  logFile.close();            // Close the log file

  return 0;
}