
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
  Double_t met_signif = 0.;
  Double_t dMetZPhi = 0.;
  Float_t frac_pT = 0.;
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
      "met_signif",
      "dMetZPhi",
      "MetOHT",
      "dLepR",
      "M2Lep",
      "leading_pT_lepton",
      "subleading_pT_lepton",
      "Z_pT",
      "n_jets",
      "n_bjets",
      "detajj",
      "mjj",
      "leading_jet_pt",
      "second_jet_pt",
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

  frac_pT = abs(met_tst - Z_pT) / Z_pT;

  Double_t signal = 0.;
  Double_t signaler = 0.;
  // Float_t Norm = 1;


  // Loop over events
  for (Int_t i = 0; i < nentries; i++)
  {
    tree->GetEntry(i);
    if (directory == "SR")
    {
      //SR SAMPLE -> FID.VOLUME: 80 < M2Lep < 100, met_tst > 70, dLepR < 1.8, dMetZPhi > 2.2, nbjets < 1, (event_3CR == 0 && (event_type == 0 || event_type == 1)
      if (dLepR < 1.8 && dMetZPhi > 2.7 && met_tst > 110 && MetOHT > 0.65 )
      {
        // Inclusive
        signal = signal + weight;
        signaler = signaler + weight * weight;
        hist->Fill(met_tst, weight);
      }
    }
    else 
    {
      //Inclusive
      signal = signal + weight;
      signaler = signaler + weight * weight;
      hist->Fill(met_tst, weight);
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
  Double_t met_signif;
  Double_t dMetZPhi;
  Float_t frac_pT;
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

  vector<TString> branches = {
      "M2Lep",
      "met_tst",
      "met_signif",
      "dMetZPhi",
      "MetOHT",
      "dLepR",
      "M2Lep",
      // "leading_pT_lepton",
      // "subleading_pT_lepton",
      "Z_pT",
      "n_jets",
      "n_bjets",
      "detajj",
      "mjj",
      // "leading_jet_pt",
      // "second_jet_pt",
      "event_3CR",
      "event_type",
      "global_weight"};


  tree->SetBranchStatus("*", 0);

  for (const auto& branch : branches) {
    tree->SetBranchStatus(branch, 1);
  }

  tree->SetBranchAddress("met_tst", &met_tst);
  tree->SetBranchAddress("met_signif", &met_signif);
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
  frac_pT = abs(met_tst - Z_pT) / Z_pT;

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
      //SR SAMPLE -> FID.VOLUME: 80 < M2Lep < 100, met_tst > 70, dLepR < 1.8, dMetZPhi > 2.2, nbjets < 1, (event_3CR == 0 && (event_type == 0 || event_type == 1)
      if (dLepR < 1.8 && dMetZPhi > 2.7 && met_tst > 110 && MetOHT > 0.65)
      {
        if (n_jets < 1)
        {
          signal0 = signal0 + weight;
          signaler0 = signaler0 + weight * weight;
          hist0->Fill(met_tst, weight);
        }
        else if (n_jets > 0 && n_jets < 2)
        {
          signal1 = signal1 + weight;
          signaler1 = signaler1 + weight * weight;
          hist1->Fill(met_tst, weight);
        }
        else if (n_jets > 1)
        {
          signal2 = signal2 + weight;
          signaler2 = signaler2 + weight * weight;
          hist2->Fill(met_tst, weight);
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
        hist0->Fill(met_tst, weight);
      }
      else if (n_jets > 0 && n_jets < 2)
      {
        signal1 = signal1 + weight;
        signaler1 = signaler1 + weight * weight;
        hist1->Fill(met_tst, weight);
      }
      else if (n_jets > 1)
      {
        signal2 = signal2 + weight;
        signaler2 = signaler2 + weight * weight;
        hist2->Fill(met_tst, weight);
      }
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
void hist_norm(TH1F *hist)
{
  Float_t integral = hist->Integral(1, hist->GetNbinsX());
  hist->Scale(1.0 / integral);
}
//
//
//
//
//
void plot_info(string directory, Int_t right_corner, TH1F *hist_data, TH1F *hist_signal, TH1F *hist_signal_reco, TH1F *hist_WZ, TH1F *hist_WW, TH1F *hist_Zjets, TH1F *hist_top, TH1F *hist_othr)
{

  TLatex *tex1;
  TLatex *tex2;
  TLatex *tex3;
  TLatex *tex4;
  TLegend *leg;

  if (right_corner == 1)
  {
    tex1 = new TLatex(0.65, 0.65, "#intL dt = 138.9 fb^{-1}");
    tex2 = new TLatex(0.65, 0.55, "#sqrt{s} = 13 TeV");
    leg = new TLegend(0.65, 0.75, 0.9, 0.85);

    if (directory == "SR")
    {
      tex3 = new TLatex(0.65, 0.45, "Signal Region");
    }
    else if (directory == "3lCR")
    {
      tex3 = new TLatex(0.6, 0.45, "3l Control Region"); 
    }
    else if (directory == "emCR_B")
    {
      tex3 = new TLatex(0.6, 0.45, "emuB Control Region"); 
    }
    else if (directory == "emCR_A")
    {
      tex3 = new TLatex(0.6, 0.45, "emuA Control Region"); 
    }
    else if (directory == "Zjets0")
    {
      tex3 = new TLatex(0.6, 0.45, "Zjets0 Control Region"); 
    }
    else if (directory == "Zjets1")
    {
      tex3 = new TLatex(0.6, 0.45, "Zjets1 Control Region");
    }
    else if (directory == "Zjets2")
    {
      tex3 = new TLatex(0.6, 0.45, "Zjets2 Control Region");
    }

    tex4 = new TLatex(0.65, 0.4, "Cuts Applied");
  }
  else
  {
    tex1 = new TLatex(0.15, 0.65, "#intL dt = 138.9 fb^{-1}");
    tex2 = new TLatex(0.15, 0.55, "#sqrt{s} = 13 TeV");
    leg = new TLegend(0.15, 0.75, 0.35, 0.85); 

    if (directory == "SR")
    {
      tex3 = new TLatex(0.15, 0.45, "Signal Region");
    }
    else if (directory == "3lCR")
    {
      tex3 = new TLatex(0.15, 0.45, "3l Control Region");
    }
    else if (directory == "emCR_B")
    {
      tex3 = new TLatex(0.15, 0.45, "emuB Control Region"); 
    }
    else if (directory == "emCR_A")
    {
      tex3 = new TLatex(0.15, 0.45, "emuA Control Region"); 
    }
    else if (directory == "Zjets0")
    {
      tex3 = new TLatex(0.15, 0.45, "Zjets0 Control Region"); 
    }
    else if (directory == "Zjets1")
    {
      tex3 = new TLatex(0.15, 0.45, "Zjets1 Control Region"); 
    }
    else if (directory == "Zjets2")
    {
      tex3 = new TLatex(0.15, 0.45, "Zjets2 Control Region");
    }

    tex4 = new TLatex(0.15, 0.4, "Cuts Applied");
  }

  if (directory == "SR")
  {
    leg->AddEntry(hist_signal_reco, "Signal", "f");
    leg->AddEntry(hist_WZ, "Backgournd", "f");
  }
  else
  {
      leg->SetNColumns(2);
      leg->AddEntry(hist_data, "Data", "p");
      leg->AddEntry(hist_signal, "Signal", "f");
      leg->AddEntry(hist_WZ, "WZ", "f");
      leg->AddEntry(hist_WW, "WW", "f");
      leg->AddEntry(hist_Zjets, "Z+jets", "f");
      leg->AddEntry(hist_top, "top", "f");
      leg->AddEntry(hist_othr, "othr", "f");
  }

  tex1->SetNDC();
  tex1->SetTextFont(42);
  tex1->SetTextSize(0.04);
  tex1->SetLineWidth(2);
  tex1->Draw();

  tex2->SetNDC();
  tex2->SetTextFont(42);
  tex2->SetTextSize(0.04);
  tex2->SetLineWidth(2);
  tex2->Draw();

  tex3->SetNDC();
  tex3->SetTextFont(1);
  tex3->SetTextSize(0.04);
  tex3->SetLineWidth(1);
  tex3->Draw();

  leg->SetBorderSize(0);
  leg->SetTextSize(0.03);
  leg->SetMargin(0.25);
  leg->SetEntrySeparation(0.1);
  leg->Draw("same");

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

  return;

}
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
void zjets_splitted_sc()
{

  // Timer start
  auto start = std::chrono::high_resolution_clock::now();

  //Output log file
  ofstream logFile("../../cro/zjets_splitted_sc/zjets_splitted_sc.txt");

  DualStreamBuffer dualBuffer(std::cout.rdbuf(), logFile.rdbuf());

  std::streambuf *oldBuffer = std::cout.rdbuf(&dualBuffer);


  gROOT->SetBatch(kTRUE);


  //-------------------------ANALYSIS-------------------------//


  vector<string> directories = {"emCR_B", "3lCR", "Zjets2", "Zjets1", "Zjets0", "emCR_A", "SR"};
  vector<string> filenames = {"DATA", "WZ", "Z_jets_ee", "Z_jets_mumu",
                              "top", "ttbarV_ttbarVV", "Wt", "WW",
                              "llll", "llqq", "VVV", "W_jets", "Ztt",
                              "lllljj", "llvv", "llvvjj", "llvvjj_WW",
                              "WZ_jj"};

  
  //Binning according to variable plotting

  Float_t xbins[21] = {70, 100, 130, 160, 190, 220, 250, 280, 310, 340, 370, 400, 430, 460, 490, 520, 550, 580, 610, 640, 700}; //met_tst
  // Float_t xbins[23] = {40, 70, 100, 130, 160, 190, 220, 250, 280, 310, 340, 370, 400, 430, 460, 490, 520, 550, 580, 610, 640, 670, 700}; //pTZ
  // Float_t xbins[13] = {70, 80, 82, 84, 86, 88, 90, 92, 94, 96, 98, 100, 110}; //M2Lep
  // Float_t xbins[18] = {0, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0}; //dLepR
  // Float_t xbins[14] = {2.0, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4}; //dMetZPhi
  // Float_t xbins[20] = {0.0, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.5}; //frac_pT
  // Float_t xbins[15] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.4, 1.6}; //MetOHT
  // Float_t xbins[16] = {0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30}; //met_signif
  // Float_t xbins[36] = {70, 80, 90, 100, 110, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 280, 310, 340, 370, 400, 430, 460, 490, 520, 550, 580, 610, 640, 670, 700}; //pTZ
  // Float_t xbins[29] = {0, 30, 50, 60, 70, 80, 90, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 345, 375, 410, 450, 500, 580, 700, 800, 1000, 1200};  //sTZ
  // Float_t xbins[27] = {100, 200, 220, 240, 260, 280, 300, 320, 345, 375, 410, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 1000, 1100, 1200, 1300, 1400, 1500};  //STjj
  // Float_t xbins[21] = {100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1050, 1100};  //mjj
  // Float_t xbins[31] = {-8, -7.5, -7.0, -6.5, -6.0, -5.5, -5.0, -4.5, -4, -3.5, -3, -2.5,-2, -1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 7.0, 9.0};  //detajj
  
  Float_t sf_3lCR, sf_3lCR_er;
  Float_t sf_emuB, sf_emuB_er;
  Float_t sf_emuA, sf_emuA_er;
  Float_t sf_Zjets0, sf_Zjets0_er;
  Float_t sf_Zjets1, sf_Zjets1_er;
  Float_t sf_Zjets2, sf_Zjets2_er;
  Float_t sf_signal, sf_signal_er;


  for (string &directory : directories)
  {
    string filepath = "/home/touk/Desktop/touk/master/thesis/data/SAMPLES/" + directory + "/";

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

    TH1F *hist_signal_reco = new TH1F("Hist_signal_reco", " ", sizeof(xbins) / sizeof(xbins[0]) - 1, xbins);
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
    Float_t signal_mc = n_llvv[0] + n_llvvjj[0];
    Float_t signal_mc_er = sqrt(pow(n_llvv[1], 2) + pow(n_llvvjj[1], 2));
    Float_t signal_reco;
    Float_t signal_reco_er;
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


    //Scaling factors calculations 

    if (directory == "emCR_B")
    {
      events_nontop = signal_mc + events_WZ + events_WW + events_Zjets0 + events_Zjets1 + events_Zjets2 + events_othr;
      events_nontop_er = sqrt(pow(signal_mc_er, 2) + pow(events_WZ_er, 2) + pow(events_WW_er, 2) + pow(events_Zjets0_er, 2) + pow(events_Zjets1_er, 2) + pow(events_Zjets2_er, 2) + pow(events_othr_er, 2));
      sf_emuB = (events_data - events_nontop) / events_top;
      sf_emuB_er = sqrt(pow((events_data - events_nontop), 2) * pow(events_top_er, 2)/pow(events_top, 4) + (pow(events_data_er, 2) + pow(events_nontop_er, 2))/pow(events_top, 2));

      cout << "   Events top / Events Non-top (%) =  " << events_top * 100 / (events_top + events_nontop) << endl << endl;
      cout << "   TOP SCALING FACTOR =  " << sf_emuB << " +- " << sf_emuB_er << endl << endl;

    }
    else if (directory == "3lCR")
    {
      //Correcting top events
      events_top = events_top * sf_emuB;
      events_top_er = sqrt( pow(sf_emuB, 2) * pow(events_top_er, 2) + pow(events_top, 2) * pow(sf_emuB_er, 2));

      events_nonWZ = signal_mc + events_WW + events_top + events_Zjets0 + events_Zjets1 + events_Zjets2 + events_othr;
      events_nonWZ_er = sqrt(pow(signal_mc_er, 2) + pow(events_WW_er, 2) + pow(events_top_er, 2) + pow(events_Zjets0_er, 2) + pow(events_Zjets1_er, 2) + pow(events_Zjets2_er, 2) + pow(events_othr_er, 2));
      sf_3lCR = (events_data - events_nonWZ) / events_WZ; 
      sf_3lCR_er = sqrt(pow((events_data - events_nonWZ), 2) * pow(events_WZ_er, 2)/pow(events_WZ, 4)  +  (pow(events_data_er, 2) + pow(events_nonWZ_er, 2))/pow(events_WZ, 2));

      cout << "   Events WZ / Events Non-WZ (%) =  " << events_WZ * 100 / (events_WZ + events_nonWZ) << endl << endl;
      cout << "   WZ SCALING FACTOR =  " << sf_3lCR << " +- " << sf_3lCR_er << endl << endl;
    
    }
    else if (directory == "Zjets2")
    {
      //Correcting for top and WZ events
      events_top = events_top * sf_emuB;
      events_top_er = sqrt( pow(sf_emuB, 2) * pow(events_top_er, 2) + pow(events_top, 2) * pow(sf_emuB_er, 2));
      events_WZ = events_WZ*sf_3lCR;
      events_WZ_er = sqrt( pow(sf_3lCR, 2) * pow(events_WZ_er, 2) + pow(events_WZ, 2) * pow(sf_3lCR_er, 2));

      events_nonZjets2 = signal_mc + events_WZ + events_top + events_WW + events_Zjets0 + events_Zjets1 + events_othr;
      events_nonZjets2_er = sqrt(pow(signal_mc_er, 2) + pow(events_WZ_er, 2) + pow(events_top_er, 2) + pow(events_WW_er, 2) + pow(events_Zjets0_er, 2) + pow(events_Zjets1_er, 2) + pow(events_othr_er, 2));
      sf_Zjets2 = (events_data - events_nonZjets2) / events_Zjets2;
      sf_Zjets2_er = sqrt(pow((events_data - events_nonZjets2), 2) * pow(events_Zjets2_er, 2)/pow(events_Zjets2, 4) + (pow(events_data_er, 2) + pow(events_nonZjets2_er, 2))/pow(events_Zjets2, 2));

      cout << "   Events Zjets2 / Events Non-Zjets2 (%) =  " << events_Zjets2 * 100/ (events_Zjets2 + events_nonZjets2) << endl << endl;
      cout << "   Zjets2 SCALING FACTOR =  " << sf_Zjets2 << " +- " << sf_Zjets2_er << endl << endl;

    }
    else if (directory == "Zjets1")
    {

      //Correcting for top, WZ and Zjets2 events
      events_top = events_top * sf_emuB;
      events_top_er = sqrt( pow(sf_emuB, 2) * pow(events_top_er, 2) + pow(events_top, 2) * pow(sf_emuB_er, 2));
      events_WZ = events_WZ*sf_3lCR;
      events_WZ_er = sqrt( pow(sf_3lCR, 2) * pow(events_WZ_er, 2) + pow(events_WZ, 2) * pow(sf_3lCR_er, 2));
      events_Zjets2 = events_Zjets2 * sf_Zjets2;
      events_Zjets2_er = sqrt( pow(sf_Zjets2, 2) * pow(events_Zjets2_er, 2) + pow(events_Zjets2, 2) * pow(sf_Zjets2_er, 2));

      events_nonZjets1 = signal_mc + events_WZ + events_top + events_WW + events_Zjets0 + events_Zjets2 + events_othr;
      events_nonZjets1_er = sqrt(pow(signal_mc_er, 2) + pow(events_WZ_er, 2) + pow(events_top_er, 2) + pow(events_WW_er, 2) + pow(events_Zjets0_er, 2) + pow(events_Zjets2_er, 2) + pow(events_othr_er, 2));
      sf_Zjets1 = (events_data - events_nonZjets1) / events_Zjets1;
      sf_Zjets1_er = sqrt(pow((events_data - events_nonZjets1), 2) * pow(events_Zjets1_er, 2)/pow(events_Zjets1, 4) + (pow(events_data_er, 2) + pow(events_nonZjets1_er, 2))/pow(events_Zjets1, 2));

      cout << "   Events Zjets1 / Events Non-Zjets1 (%) =  " << events_Zjets1 * 100 / (events_Zjets1 + events_nonZjets1) << endl << endl;
      cout << "   Zjets1 SCALING FACTOR =  " << sf_Zjets1 << " +- " << sf_Zjets1_er << endl << endl;

    }
    else if (directory == "Zjets0")
    {

      //Correcting for top, WZ, Zjets2 and Zjets1 events
      events_top = events_top * sf_emuB;
      events_top_er = sqrt( pow(sf_emuB, 2) * pow(events_top_er, 2) + pow(events_top, 2) * pow(sf_emuB_er, 2));
      events_WZ = events_WZ*sf_3lCR;
      events_WZ_er = sqrt( pow(sf_3lCR, 2) * pow(events_WZ_er, 2) + pow(events_WZ, 2) * pow(sf_3lCR_er, 2));
      events_Zjets2 = events_Zjets2 * sf_Zjets2;
      events_Zjets2_er = sqrt( pow(sf_Zjets2, 2) * pow(events_Zjets2_er, 2) + pow(events_Zjets2, 2) * pow(sf_Zjets2_er, 2));
      events_Zjets1 = events_Zjets1 * sf_Zjets1;
      events_Zjets1_er = sqrt( pow(sf_Zjets1, 2) * pow(events_Zjets1_er, 2) + pow(events_Zjets1, 2) * pow(sf_Zjets1_er, 2));

      events_nonZjets0 = signal_mc + events_WZ + events_top + events_WW + events_Zjets1 + events_Zjets2 + events_othr;
      events_nonZjets0_er = sqrt(pow(signal_mc_er, 2) + pow(events_WZ_er, 2) + pow(events_top_er, 2) + pow(events_WW_er, 2) + pow(events_Zjets1_er, 2) + pow(events_Zjets2_er, 2)   + pow(events_othr_er, 2));
      sf_Zjets0 = (events_data - events_nonZjets0) / events_Zjets0;
      sf_Zjets0_er = sqrt(pow((events_data - events_nonZjets0), 2) * pow(events_Zjets0_er, 2)/pow(events_Zjets0, 4) + (pow(events_data_er, 2) + pow(events_nonZjets0_er, 2))/pow(events_Zjets0, 2));

      cout << "   Events Zjets0 / Events Non-Zjets0 (%) =  " << events_Zjets0 * 100/ (events_Zjets0 + events_nonZjets0) << endl << endl;
      cout << "   Zjets0 SCALING FACTOR =  " << sf_Zjets0 << " +- " << sf_Zjets0_er << endl << endl;

    }
    else if (directory == "emCR_A")
    {

      //Correcting for top, WZ, Zjets2, Zjets1 and Zjets0 events
      events_top = events_top * sf_emuB;
      events_top_er = sqrt( pow(sf_emuB, 2) * pow(events_top_er, 2) + pow(events_top, 2) * pow(sf_emuB_er, 2));
      events_WZ = events_WZ*sf_3lCR;
      events_WZ_er = sqrt( pow(sf_3lCR, 2) * pow(events_WZ_er, 2) + pow(events_WZ, 2) * pow(sf_3lCR_er, 2));
      events_Zjets2 = events_Zjets2 * sf_Zjets2;
      events_Zjets2_er = sqrt( pow(sf_Zjets2, 2) * pow(events_Zjets2_er, 2) + pow(events_Zjets2, 2) * pow(sf_Zjets2_er, 2));
      events_Zjets1 = events_Zjets1 * sf_Zjets1;
      events_Zjets1_er = sqrt( pow(sf_Zjets1, 2) * pow(events_Zjets1_er, 2) + pow(events_Zjets1, 2) * pow(sf_Zjets1_er, 2));
      events_Zjets0 = events_Zjets0 * sf_Zjets0;
      events_Zjets0_er = sqrt( pow(sf_Zjets0, 2) * pow(events_Zjets0_er, 2) + pow(events_Zjets0, 2) * pow(sf_Zjets0_er, 2));

      events_nonWW = signal_mc + events_WZ + events_top + events_Zjets0 + events_Zjets1 + events_Zjets2 + events_othr;
      events_nonWW_er = sqrt(pow(signal_mc_er, 2) + pow(events_WZ_er, 2) + pow(events_top_er, 2) + pow(events_Zjets0_er, 2) + pow(events_Zjets1_er, 2) + pow(events_Zjets2_er, 2) + pow(events_othr_er, 2));
      sf_emuA = (events_data - events_nonWW) / events_WW;
      sf_emuA_er = sqrt(pow((events_data - events_nonWW), 2) * pow(events_WW_er, 2)/pow(events_WW, 4) + (pow(events_data_er, 2) + pow(events_nonWW_er, 2))/pow(events_WW, 2));

      cout << "   Events WW / Events Non-WW (%) =  " << events_WW * 100 / (events_WW + events_nonWW) << endl << endl;
      cout << "   WW SCALING FACTOR =  " << sf_emuA << " +- " << sf_emuA_er << endl << endl;

    }
    else if (directory == "SR")
    {
      //Correcting for top, WZ, Zjets2, Zjets1, Zjets0 and WW events
      events_top = events_top * sf_emuB;
      events_top_er = sqrt( pow(sf_emuB, 2) * pow(events_top_er, 2) + pow(events_top, 2) * pow(sf_emuB_er, 2));
      events_WZ = events_WZ*sf_3lCR;
      events_WZ_er = sqrt( pow(sf_3lCR, 2) * pow(events_WZ_er, 2) + pow(events_WZ, 2) * pow(sf_3lCR_er, 2));
      events_Zjets2 = events_Zjets2 * sf_Zjets2;
      events_Zjets2_er = sqrt( pow(sf_Zjets2, 2) * pow(events_Zjets2_er, 2) + pow(events_Zjets2, 2) * pow(sf_Zjets2_er, 2));
      events_Zjets1 = events_Zjets1 * sf_Zjets1;
      events_Zjets1_er = sqrt( pow(sf_Zjets1, 2) * pow(events_Zjets1_er, 2) + pow(events_Zjets1, 2) * pow(sf_Zjets1_er, 2));
      events_Zjets0 = events_Zjets0 * sf_Zjets0;
      events_Zjets0_er = sqrt( pow(sf_Zjets0, 2) * pow(events_Zjets0_er, 2) + pow(events_Zjets0, 2) * pow(sf_Zjets0_er, 2));
      events_WW = events_WW * sf_emuA;
      events_WW_er = sqrt( pow(sf_emuA, 2) * pow(events_WW_er, 2) + pow(events_WW, 2) * pow(sf_emuA_er, 2));
      
      //Signal
      events_bkg = events_WZ + events_top + events_WW + events_Zjets0 + events_Zjets1 + events_Zjets2 + events_othr;
      events_bkg_er = sqrt(pow(events_WZ_er, 2) + pow(events_top_er, 2) + pow(events_WW_er, 2) + pow(events_Zjets0_er, 2) + pow(events_Zjets1_er, 2) + pow(events_Zjets2_er, 2) + pow(events_othr_er, 2));
      signal_reco = events_data - events_bkg;
      signal_reco_er = sqrt(pow(events_data_er, 2) + pow(events_bkg_er, 2));
      sf_signal = signal_reco / signal_mc;
      sf_signal_er = sqrt(pow(signal_reco, 2) * pow(signal_mc_er, 2)/pow(signal_mc, 4) + pow(signal_reco_er, 2)/pow(signal_mc, 2));
      cout << "   Zjets1 SCALING FACTOR =  " << sf_Zjets1 << " +- " << sf_Zjets1_er << endl << endl;
    }
    
    //Scale histograms before adding and plotting
    if (directory == "emCR_B")
    {
       for (int bin = 1; bin < sizeof(xbins) / sizeof(xbins[0]); bin++)
       {
        hist_WZ->SetBinContent(bin, hist_WZ->GetBinContent(bin) * sf_3lCR);
       }
    }
    if (directory == "emCR_A")
    {
       for (int bin = 1; bin < sizeof(xbins) / sizeof(xbins[0]); bin++)
       {
        hist_WZ->SetBinContent(bin, hist_WZ->GetBinContent(bin) * sf_3lCR);
        hist_top->SetBinContent(bin, hist_top->GetBinContent(bin) * sf_emuB);
       }
    }
    if (directory == "Zjets0")
    {
       for (int bin = 1; bin < sizeof(xbins) / sizeof(xbins[0]); bin++)
       {
        hist_WZ->SetBinContent(bin, hist_WZ->GetBinContent(bin) * sf_3lCR);
        hist_top->SetBinContent(bin, hist_top->GetBinContent(bin) * sf_emuB);
        hist_WW->SetBinContent(bin, hist_WW->GetBinContent(bin) * sf_emuA);
       }
    }
    if (directory == "Zjets1")
    {
       for (int bin = 1; bin < sizeof(xbins) / sizeof(xbins[0]); bin++)
       {
        hist_WZ->SetBinContent(bin, hist_WZ->GetBinContent(bin) * sf_3lCR);
        hist_top->SetBinContent(bin, hist_top->GetBinContent(bin) * sf_emuB);
        hist_WW->SetBinContent(bin, hist_WW->GetBinContent(bin) * sf_emuA);
        hist_Zjets0->SetBinContent(bin, hist_Zjets0->GetBinContent(bin) * sf_Zjets0);
       }
    }
    else if (directory == "Zjets2")
    {
       for (int bin = 1; bin < sizeof(xbins) / sizeof(xbins[0]); bin++)
       {
        hist_WZ->SetBinContent(bin, hist_WZ->GetBinContent(bin) * sf_3lCR);
        hist_top->SetBinContent(bin, hist_top->GetBinContent(bin) * sf_emuB);
        hist_WW->SetBinContent(bin, hist_WW->GetBinContent(bin) * sf_emuA);
        hist_Zjets0->SetBinContent(bin, hist_Zjets0->GetBinContent(bin) * sf_Zjets0);
        hist_Zjets1->SetBinContent(bin, hist_Zjets1->GetBinContent(bin) * sf_Zjets1);
       }
    }
    else if (directory == "SR")
    {
       for (int bin = 1; bin < sizeof(xbins) / sizeof(xbins[0]); bin++)
       {
        hist_WZ->SetBinContent(bin, hist_WZ->GetBinContent(bin) * sf_3lCR);
        hist_top->SetBinContent(bin, hist_top->GetBinContent(bin) * sf_emuB);
        hist_WW->SetBinContent(bin, hist_WW->GetBinContent(bin) * sf_emuA);
        hist_Zjets0->SetBinContent(bin, hist_Zjets0->GetBinContent(bin) * sf_Zjets0);
        hist_Zjets1->SetBinContent(bin, hist_Zjets1->GetBinContent(bin) * sf_Zjets1);
        hist_Zjets2->SetBinContent(bin, hist_Zjets2->GetBinContent(bin) * sf_Zjets2);
       }
    }

    //Merge Zjets before start plotting
    hist_Zjets->Add(hist_Zjets0);
    hist_Zjets->Add(hist_Zjets1);
    hist_Zjets->Add(hist_Zjets2);

    //Calculate signal histogram
    
    Double_t Z_max;
    Double_t Bin_max;
    if (directory == "SR")
    {
       // Stacking with a specific order
       hist_Zjets->Add(hist_othr);
       hist_top->Add(hist_Zjets);
       hist_WW->Add(hist_top);
       hist_WZ->Add(hist_WW);
       cout << "   SIGNAL/BKG = " << (hist_data->Integral(1, hist_data->GetNbinsX()) - hist_WZ->Integral(1, hist_WZ->GetNbinsX())) / hist_WZ->Integral(1, hist_WZ->GetNbinsX()) << endl << endl;
      //  cout << "   SIGNAL/BKG = " << hist_signal->Integral(1, hist_signal->GetNbinsX()) / hist_WZ->Integral(1, hist_WZ->GetNbinsX()) << endl << endl;
      
      
      //Reconstructed sugnal events
      hist_signal_reco->Add(hist_data);
      hist_signal_reco->Add(hist_WZ, -1);

      for (int bin = 1; bin < hist_signal_reco->GetSize(); ++bin)
      {
       Double_t B = hist_WZ->GetBinContent(bin);
       Double_t S = hist_signal_reco->GetBinContent(bin);

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

      cout << "   SIGNAL INTEGRAL IS:   " << hist_signal_reco->Integral(1, hist_signal_reco->GetNbinsX()) << endl << endl;

      //Print calculated events for every region
    
      cout << "------------------------------------------------------------------" << endl << endl;
      cout << "   SIGNAL MC       =  " << signal_mc << " +- " << signal_mc_er  << endl;
      cout << "   SIGNAL RECO     =  " << signal_reco << " +- " << signal_reco_er  << endl;
      cout << "   SIGNAL RECO/BKG =  " << signal_reco/events_bkg << endl << endl;
      cout << "   Maximum significance is:  " << Z_max << " for bin:   " << Bin_max << endl << endl;
      cout << "_________________________________" << endl << endl;
      cout << "   Data: " << "     " << events_data << " +- " << events_data_er << endl << endl;
      cout << "   WZ: " << "       " << events_WZ << " +- " << events_WZ_er << "   |   mu_WZ = " << sf_3lCR << " +- " << sf_3lCR_er << endl << endl;
      cout << "   Top: " << "      " << events_top << " +- " << events_top_er << "   |   mu_top = " << sf_emuB << " +- " << sf_emuB_er << endl << endl;
      cout << "   WW: " << "       " << events_WW << " +- " << events_WW_er << "   |   mu_WW = " << sf_emuA << " +- " << sf_emuA_er << endl << endl;
      cout << "   Zjets0: " << "   " << events_Zjets0 << " +- " << events_Zjets0_er << "   |   mu_Zjets0 = " << sf_Zjets0 << " +- " << sf_Zjets0_er << endl << endl;
      cout << "   Zjets1: " << "   " << events_Zjets1 << " +- " << events_Zjets1_er << "   |   mu_Zjets1 = " << sf_Zjets1 << " +- " << sf_Zjets1_er << endl << endl;
      cout << "   Zjets2: " << "   " << events_Zjets2 << " +- " << events_Zjets2_er << "   |   mu_Zjets2 = " << sf_Zjets2 << " +- " << sf_Zjets2_er << endl << endl;
      cout << "   Other: " << "    " << events_othr << " +- " << events_othr_er << endl << endl;
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
      cout << "   DATA/MC = " << hist_data->Integral(1, hist_data->GetNbinsX())/ hist_signal->Integral(1, hist_signal->GetNbinsX()) << endl << endl;
    }

    //----------------------------------------PLOTS----------------------------------------//


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
    pad1->SetBottomMargin(0.045);
    pad2->SetBottomMargin(0.35);
    pad2->SetTopMargin(0.0);
    pad2->SetBorderSize(0);
    pad1->Draw();
    pad2->Draw();

    gStyle->SetOptTitle(kTRUE);
    gStyle->SetLegendFont(1);
    gStyle->SetPalette(20);

    //---------------PAD 1---------------//

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
      hist_signal_reco->Draw("histsame");

      //Normalise to unity
      hist_norm(hist_WZ);
      hist_norm(hist_signal_reco);
      
      if (hist_signal_reco->GetMaximum() > hist_WZ->GetMaximum())
      {
       hist_WZ->GetYaxis()->SetRangeUser(0, hist_signal_reco->GetMaximum() * 1.4);
      }
      else
      {
       hist_WZ->GetYaxis()->SetRangeUser(0, hist_WZ->GetMaximum() * 1.4);
      }
      hist_WZ->SetStats(0);
      hist_WZ->SetLineWidth(2);
      hist_WZ->SetLineColor(kBlue);
      hist_WZ->SetFillColorAlpha(TColor::GetColor("#6495ED"), 0.9);
      hist_WZ->SetLineWidth(2);
      hist_WZ->GetYaxis()->SetTitle("Events");
      hist_WZ->GetXaxis()->SetLabelSize(0);

      hist_signal_reco->SetLineColor(TColor::GetColor("#DE3163"));
      hist_signal_reco->SetFillColorAlpha(TColor::GetColor("#DE3163"), 0.8);
      hist_signal_reco->SetLineWidth(2);
      hist_signal_reco->SetFillStyle(3244);
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

      // //Normalise to unity
      // hist_norm(hist_signal);
      // hist_norm(hist_WZ);
      // hist_norm(hist_WW);
      // hist_norm(hist_top);
      // hist_norm(hist_Zjets);
      // hist_norm(hist_othr);
      // hist_norm(hist_data);

      hist_signal->GetXaxis()->SetTitleSize(0.03);
      hist_signal->GetXaxis()->SetTitleFont(42);
      hist_signal->GetXaxis()->SetTitleOffset(1.1);
      hist_signal->GetYaxis()->SetTitleFont(42);
      hist_signal->GetYaxis()->SetTitle("Events");
      hist_signal->GetYaxis()->SetRangeUser(0, hist_signal->GetMaximum() * 1.4);
      hist_signal->SetStats(0);
      hist_signal->GetXaxis()->SetLabelSize(0);
    }


    Int_t right_corner = 0; //For use to choose corner to plot information
    plot_info(directory, right_corner, hist_data, hist_signal, hist_signal_reco, hist_WZ, hist_WW, hist_Zjets, hist_top, hist_othr);

    pad1->RedrawAxis();



    //---------------PAD 2---------------//

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
      hist_signal_sig->GetYaxis()->SetTitle("Significance");
      hist_signal_sig->GetXaxis()->SetTitle("met_tst");
      pad2->SetGrid();
    }
    else
    {
      numerator->Add(hist_data);
      denominator->Add(hist_signal);
      // pad2->SetLogy();
    }
    

    if (directory != "SR")
    {
      // numerator->Sumw2(1);
      // denominator->Sumw2(1);
      TLine *line = new TLine(xbins[0], 1, xbins[sizeof(xbins) / sizeof(xbins[0]) - 1], 1);
      line->SetLineStyle(2);
      line->SetLineWidth(2);

      numerator->SetMarkerStyle(20);
      numerator->GetXaxis()->SetTitle("");
      numerator->GetXaxis()->SetTitleSize(0.15);
      numerator->GetXaxis()->SetTitle("met_tst");
      numerator->GetXaxis()->SetTitleOffset(1.1);
      numerator->GetXaxis()->SetLabelSize(0.14);

      numerator->GetYaxis()->SetTitleOffset(0.4);
      numerator->GetYaxis()->SetTitleSize(0.12);
      numerator->GetYaxis()->SetTitle("#frac{Data}{MC}");
      numerator->GetYaxis()->SetLabelSize(0.10);
      numerator->SetStats(0);

      numerator->SetTitle(" ");
      numerator->SetLineColor(kBlack);
      numerator->Divide(denominator);

      numerator->Draw();
      line->Draw("same");
      pad2->RedrawAxis();
    }

    
    
    if (directory == "SR")
    {
      c1->SaveAs("../../cro/zjets_splitted_sc/ptz_SR_splitted_sc.png");
    }
    else if (directory == "3lCR")
    {
      c2->SaveAs("../../cro/zjets_splitted_sc/ptz_3lCR_splitted_sc.png");
    }
    else if (directory == "emCR_B")
    {
      c3->SaveAs("../../cro/zjets_splitted_sc/ptz_emCR_B_splitted_sc.png");
    }
    else if (directory == "emCR_A")
    {
      c4->SaveAs("../../cro/zjets_splitted_sc/ptz_emCR_A_splitted_sc.png");
    }
    else if (directory == "Zjets0")
    {
      c5->SaveAs("../../cro/zjets_splitted_sc/ptz_Zjets0_splitted_sc.png");
    }
    else if (directory == "Zjets1")
    {
      c6->SaveAs("../../cro/zjets_splitted_sc/ptz_Zjets1_splitted_sc.png");
    }
    else if (directory == "Zjets2")
    {
      c7->SaveAs("../../cro/zjets_splitted_sc/ptz_Zjets2_splitted_sc.png");
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
  
  cout << endl << endl;
  cout << "   SCALING FACTORS:       μ_WZ = " <<  sf_3lCR << " +- " << sf_3lCR_er << endl << endl 
       << "                          μ_top = " << sf_emuB  << " +- " << sf_emuB_er << endl << endl 
       << "                          μ_WW = " << sf_emuA << " +- " << sf_emuA_er << endl << endl
       << "                          μ_Zjets0 = " << sf_Zjets0 << " +- " << sf_Zjets0_er << endl << endl  
       << "                          μ_Zjets1 = " << sf_Zjets1 << " +- " << sf_Zjets1_er << endl << endl 
       << "                          μ_Zjets2 = " << sf_Zjets2 << " +- " << sf_Zjets2_er << endl << endl
       << "                          μ_signal = " << sf_signal << " +- " << sf_signal_er << endl << endl;
          
          
  // Timer stop
  auto end = chrono::high_resolution_clock::now();
  chrono::duration<float> duration = end - start;
  cout << endl << "   Script executed in: " << int(duration.count() / 60.0) << " minutes" << " and " << int((duration.count() / 60.0 - int(duration.count() / 60.0))*60) << " s" <<  endl << endl;

  // For the log file
  std::cout.rdbuf(oldBuffer); // Restore the original stream buffer for std::cout
  logFile.close(); // Close the log file

  return;
  
}
