
#include <TRandom3.h>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <TF1.h>
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
  Double_t mTW = 0.;
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
      // "M2Lep",
      "leading_pT_lepton",
      "subleading_pT_lepton",
      "Z_pT",
      "n_jets",
      "n_bjets",
      "Z_rapidity",
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
  // tree->SetBranchAddress("mTW", &mTW);

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
     if (dMetZPhi > 2.6898 && met_tst > 101.54 && MetOHT > 0.7375)
      {
        // Inclusive
        signal = signal + weight;
        signaler = signaler + weight * weight;
        hist->Fill(dMetZPhi, weight);
      }
    }
    else if (directory == "Zjets0" || directory == "Zjets1" || directory == "Zjets2")
    {
      if ((met_tst < 80 && MetOHT > 0.3) || (met_tst > 80 && MetOHT > 0.3 && MetOHT < 0.45))
      {
        // Inclusive
        signal = signal + weight;
        signaler = signaler + weight * weight;
        hist->Fill(dMetZPhi, weight);
      }
    }
    // else
    else if (directory == "emCR_B" || directory == "emCR_A" || directory == "3lCR")
    {
        // Inclusive
        signal = signal + weight;
        signaler = signaler + weight * weight;
        hist->Fill(dMetZPhi, weight);
    }
    else if (directory == "VR")
    {


      if ((met_tst > 80 && met_tst < 90 && MetOHT > 0.45) || (met_tst > 90 && MetOHT > 0.45 && MetOHT < 0.5))        //VR3: 3517, 9.25%      

      {
        signal = signal + weight;
        signaler = signaler + weight * weight;
        hist->Fill(dMetZPhi, weight);
      }
      
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
  Double_t mTW;
  vector<Float_t> events;
  events.clear();

  vector<TString> branches = {
      "M2Lep",
      "met_tst",
      "met_signif",
      "dMetZPhi",
      "MetOHT",
      "dLepR",
      // "M2Lep",
      "leading_pT_lepton",
      "subleading_pT_lepton",
      "Z_pT",
      "n_jets",
      "n_bjets",
      "Z_rapidity",
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
  // tree->SetBranchAddress("mTW", &mTW);
  frac_pT = abs(met_tst - Z_pT) / Z_pT;

  Float_t signal0 = 0;
  Float_t signal1 = 0;
  Float_t signal2 = 0;
  Float_t signal_tot = 0;
  Float_t signal_tot_er = 0;
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
      if (dMetZPhi > 2.6898 && met_tst > 101.54 && MetOHT > 0.7375)
      {
        if (n_jets < 1)
        {
          signal0 = signal0 + weight;
          signaler0 = signaler0 + weight * weight;
          hist0->Fill(dMetZPhi, weight);
        }
        else if (n_jets > 0 && n_jets < 2)
        {
          signal1 = signal1 + weight;
          signaler1 = signaler1 + weight * weight;
          hist1->Fill(dMetZPhi, weight);
        }
        else if (n_jets > 1)
        {
          signal2 = signal2 + weight;
          signaler2 = signaler2 + weight * weight;
          hist2->Fill(dMetZPhi, weight);
        }
      }
    }
    else if (directory == "Zjets0" || directory == "Zjets1" || directory == "Zjets2")
    {
      if ((met_tst < 80 && MetOHT > 0.3) || (met_tst > 80 && MetOHT > 0.3 && MetOHT < 0.45))
      {
        // Inclusive
        if (n_jets < 1)
        {
          signal0 = signal0 + weight;
          signaler0 = signaler0 + weight * weight;
          hist0->Fill(dMetZPhi, weight);
        }
        else if (n_jets > 0 && n_jets < 2)
        {
          signal1 = signal1 + weight;
          signaler1 = signaler1 + weight * weight;
          hist1->Fill(dMetZPhi, weight);
        }
        else if (n_jets > 1)
        {
          signal2 = signal2 + weight;
          signaler2 = signaler2 + weight * weight;
          hist2->Fill(dMetZPhi, weight);
        }
      }
    }
    // else
    else if (directory == "emCR_B" || directory == "emCR_A" || directory == "3lCR")
    {
      // Inclusive
      if (n_jets < 1)
      {
        signal0 = signal0 + weight;
        signaler0 = signaler0 + weight * weight;
        hist0->Fill(dMetZPhi, weight);
      }
      else if (n_jets > 0 && n_jets < 2)
      {
        signal1 = signal1 + weight;
        signaler1 = signaler1 + weight * weight;
        hist1->Fill(dMetZPhi, weight);
      }
      else if (n_jets > 1)
      {
        signal2 = signal2 + weight;
        signaler2 = signaler2 + weight * weight;
        hist2->Fill(dMetZPhi, weight);
      }
    }
    else if (directory == "VR")
    {
      if ((met_tst > 80 && met_tst < 90 && MetOHT > 0.45) || (met_tst > 90 && MetOHT > 0.45 && MetOHT < 0.5))        //VR3: 3517, 9.25%
      {
        if (n_jets < 1)
        {
          signal0 = signal0 + weight;
          signaler0 = signaler0 + weight * weight;
          hist0->Fill(dMetZPhi, weight);
        }
        else if (n_jets > 0 && n_jets < 2)
        {
          signal1 = signal1 + weight;
          signaler1 = signaler1 + weight * weight;
          hist1->Fill(dMetZPhi, weight);
        }
        else if (n_jets > 1)
        {
          signal2 = signal2 + weight;
          signaler2 = signaler2 + weight * weight;
          hist2->Fill(dMetZPhi, weight);
        }
      }
    }
  }

  signal_tot = signal0 + signal1 + signal2;
  signal_tot_er = signaler0 + signaler1 + signaler2;



  cout << "     ENTRIES = " << tree->GetEntries() << endl << endl; 
  cout << "          N0 = " << signal0    << "+-" << sqrt(signaler0) << endl << endl; 
  cout << "          N1 = " << signal1    << "+-" << sqrt(signaler1) << endl << endl; 
  cout << "          N2 = " << signal2    << "+-" << sqrt(signaler2) << endl << endl; 


  events.push_back(signal0);
  events.push_back(sqrt(signaler0));
  events.push_back(signal1);
  events.push_back(sqrt(signaler1));
  events.push_back(signal2);
  events.push_back(sqrt(signaler2));
  events.push_back(signal_tot);
  events.push_back(sqrt(signal_tot_er));

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
void plot_info(string directory, Int_t right_corner, TF1 *fitFunc, TH1F *hist_data, TH1F *hist_signal, TH1F *hist_signal_obs, TH1F *hist_WZ, TH1F *hist_WW, TH1F *hist_Zjets, TH1F *hist_top, TH1F *hist_other, TH1F *hist_uncert_dummy)
{

  TLatex *tex1;
  TLatex *tex2;
  TLatex *tex3;
  TLatex *tex4;
  TLegend *leg;


  if (right_corner == 2)
  {
    tex1 = new TLatex(0.12, 0.82, "ATLAS");
    tex2 = new TLatex(0.12, 0.82, "#sqrt{s}  = 13 TeV, 140.1 fb^{-1}");    
    leg = new TLegend(0.44, 0.7, 0.87, 0.87);
    // string text = "Pre-fit";
    // string text = "Post-fit";
    string text = "Direct Est.";

    if (directory == "SR")
    {
      tex3 = new TLatex(0.12, 0.76, "Signal Region");
      tex4 = new TLatex(0.12, 0.71, "Data Blinded");
    }
    else if (directory == "3lCR")
    {
      tex3 = new TLatex(0.12, 0.76, ("3l CR, " + text).c_str());
    }
    else if (directory == "emCR_B")
    {
      tex3 = new TLatex(0.12, 0.76, ("emuB CR, " + text).c_str());
    }
    else if (directory == "emCR_A")
    {
      tex3 = new TLatex(0.12, 0.76, ("emuA CR, " + text).c_str());
    }
    else if (directory == "Zjets0")
    {
      tex3 = new TLatex(0.12, 0.76, ("Zjets0 CR, " + text).c_str());
    }
    else if (directory == "Zjets1")
    {
      tex3 = new TLatex(0.12, 0.76, ("Zjets1 CR, " + text).c_str());
    }
    else if (directory == "Zjets2")
    {
      tex3 = new TLatex(0.12, 0.76, ("Zjets2 CR, " + text).c_str());
    }
    else if (directory == "VR")
    {
      tex3 = new TLatex(0.12, 0.76, "Validation Region");
      tex4 = new TLatex(0.12, 0.70, (text).c_str());
    }
    

    // if (directory == "SR")
    // {
    //  tex4 = new TLatex(0.7, 0.55, "Data Blinded");
    // }
    // else
    // {
    //   tex4 = new TLatex(0.7, 0.55, "Pre-fit");
    // }
  }

  leg->SetNColumns(2);
  leg->AddEntry(hist_data, "Data", "p");
  leg->AddEntry(hist_signal, "ZZ", "f");
  leg->AddEntry(hist_WZ, "WZ", "f");
  leg->AddEntry(hist_WW, "WW", "f");
  leg->AddEntry(hist_top, "top", "f");
  leg->AddEntry(hist_Zjets, "Z+jets", "f");
  leg->AddEntry(hist_other, "other", "f");
  leg->AddEntry(hist_uncert_dummy, "Uncertainty", "f");
  leg->AddEntry(fitFunc, " Fit");
  
  tex1->SetNDC();
  tex1->SetTextSize(0.04);
  tex1->SetTextFont(72);
  tex1->SetLineWidth(2);
  // tex1->Draw();

  tex2->SetNDC();
  tex2->SetTextSize(0.04);
  tex2->SetLineWidth(2);
  tex2->Draw();

  tex3->SetNDC();
  tex3->SetTextSize(0.04);
  tex3->SetLineWidth(1);
  tex3->Draw();

  if (directory == "SR" || directory == "VR")
  {
    tex4->SetNDC();
    tex4->SetTextSize(0.04);
    tex4->SetLineWidth(1);
    tex4->Draw();
  }

  leg->SetBorderSize(0);
  leg->SetTextSize(0.04);
  leg->SetMargin(0.4);
  leg->SetEntrySeparation(0.25);
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
  // leg->AddEntry(hist_other, "other", "f");
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
  std::streambuf *primaryBuffer;   
  std::streambuf *secondaryBuffer; 
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
void regions_plotter()
{

  // Timer start
  auto start = std::chrono::high_resolution_clock::now();

  //Output log file
  ofstream logFile("./zjets_splitted_sc/met_signif_direct.txt");

  DualStreamBuffer dualBuffer(std::cout.rdbuf(), logFile.rdbuf());

  std::streambuf *oldBuffer = std::cout.rdbuf(&dualBuffer);


  gROOT->SetBatch(kTRUE);


  //-------------------------ANALYSIS-------------------------//


  vector<string> directories = {"emCR_B", "emCR_A", "3lCR", "Zjets2",  "Zjets1", "Zjets0", "SR", "VR"};
  vector<string> filenames = {"DATA", "WZ", "Z_jets_ee", "Z_jets_mumu",
                              "top", "ttbarV_ttbarVV", "Wt", "WW",
                              "llll", "llqq", "VVV", "W_jets", "Ztt",
                              "lllljj", "llvv", "llvvjj", "llvvjj_WW",
                              "WZ_jj"};

  
  //Binning according to variable plotting

  // Float_t xbins[21] = {70, 100, 130, 160, 190, 220, 250, 280, 310, 340, 370, 400, 430, 460, 490, 520, 550, 580, 610, 640, 700}; //met_tst
  // Float_t xbins[19] = {50, 70, 90, 110, 130, 150, 170, 190, 210, 230, 250, 270, 290, 310, 330, 350, 370, 390, 410}; //plot met_tst
  // Float_t xbins[23] = {50, 70, 80, 90, 100, 110, 120, 130, 140, 150, 170, 190, 210, 230, 250, 270, 290, 310, 330, 350, 370, 390, 410}; //plot met_tst
  // Float_t xbins[20] = {50, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250}; //plot met_tst
  // Float_t xbins[10] = {50, 70, 80, 90, 100, 110, 120, 130, 140, 150}; //plot met_tst
  // Float_t xbins[4] = {0, 70, 600, 1000}; //met_tst for nocut unsc
  // Float_t xbins[23] = {40, 70, 100, 130, 160, 190, 220, 250, 280, 310, 340, 370, 400, 430, 460, 490, 520, 550, 580, 610, 640, 670, 700}; //pTZ
  // Float_t xbins[20] = {30, 50, 70, 90, 110, 130, 150, 170, 190, 210, 230, 250, 270, 290, 310, 330, 350, 370, 390, 410}; //plot pTZ
  // Float_t xbins[13] = {70, 80, 82, 84, 86, 88, 90, 92, 94, 96, 98, 100, 110}; //M2Lep
  // Float_t xbins[12] = {0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2}; //dLepR
  Float_t xbins[13] = {2.0, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.4}; //dMetZPhi
  // Float_t xbins[8] = {2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4}; //dMetZPhi
  // Float_t xbins[4] = {2, 2.2, 3.2, 3.4 }; //dMetZPhi 1-bin
  // Float_t xbins[12] = {-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10}; //Njets
  // Float_t xbins[20] = {0.0, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.5}; //frac_pT
  // Float_t xbins[15] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.4, 1.6}; //MetOHT
  // Float_t xbins[8] = {0., 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5}; //plot MetOHT step 0.2
  // Float_t xbins[13] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.4}; //plot MetOHT step 0.1
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

    // other
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

    // other = llll, llqq, VVV, Wjets, Ztt
    TH1F *hist_other = new TH1F("Hist_other", " ", sizeof(xbins) / sizeof(xbins[0]) - 1, xbins);
    TH1F *hist_llll = new TH1F("Hist_llll", " ", sizeof(xbins) / sizeof(xbins[0]) - 1, xbins);
    TH1F *hist_llqq = new TH1F("Hist_llqq", " ", sizeof(xbins) / sizeof(xbins[0]) - 1, xbins);
    TH1F *hist_VVV = new TH1F("Hist_VVV", " ", sizeof(xbins) / sizeof(xbins[0]) - 1, xbins);
    TH1F *hist_W_jets = new TH1F("Hist_W_jets", " ", sizeof(xbins) / sizeof(xbins[0]) - 1, xbins);
    TH1F *hist_Ztt = new TH1F("Hist_Ztt", " ", sizeof(xbins) / sizeof(xbins[0]) - 1, xbins);
    TH1F *hist_WZ_jj = new TH1F("Hist_WZ_jj", " ", sizeof(xbins) / sizeof(xbins[0]) - 1, xbins);
    TH1F *hist_lllljj = new TH1F("Hist_lllljj", " ", sizeof(xbins) / sizeof(xbins[0]) - 1, xbins);
    TH1F *hist_llvvjj_WW = new TH1F("Hist_llvvjj_WW", " ", sizeof(xbins) / sizeof(xbins[0]) - 1, xbins);

    TH1F *hist_signal_obs = new TH1F("Hist_signal_obs", " ", sizeof(xbins) / sizeof(xbins[0]) - 1, xbins);
    TH1F *hist_signal_sig = new TH1F("hist_signal_sig ", " ", sizeof(xbins) / sizeof(xbins[0]) - 1, xbins);

    TH1F *hist_uncert_dummy = new TH1F("hist_uncert_dummy ", " ", sizeof(xbins) / sizeof(xbins[0]) - 1, xbins);

    TF1 *fitFunc = new TF1("fitFunc", "pol1");

 
    cout << "   ================== DATA ==================    " << endl << endl;
    cout << "   DATA:";
    vector<Float_t> n_data = Counter(tree_data, hist_data, directory);

    cout << "   ================== SIGNAL ==================    " << endl << endl;
    cout << "   llvv   :";
    vector<Float_t> n_llvv = Counter(tree_llvv, hist_llvv, directory);
    cout << "   llvvjj :";
    vector<Float_t> n_llvvjj = Counter(tree_llvvjj, hist_llvvjj, directory);
    cout << "   QCD+EWK:" << n_llvv[0] + n_llvvjj[0] << "+-" << sqrt(pow(n_llvv[1],2) + pow(n_llvvjj[1],2)) << endl << endl;

    cout << "   ================== WZ ==================    " << endl << endl;
    cout << "   WZ:";
    vector<Float_t> n_WZ = Counter(tree_WZ, hist_WZ, directory);

    cout << "   ================== Zjets ==================    " << endl << endl;
    // cout << "   Z_jets_ee:";
    // vector<Float_t> n_Zjets_ee = Counter(tree_Z_jets_ee, hist_Zjets_ee, directory);
    // cout << "   Z_jets_mumu:";
    // vector<Float_t> n_Zjets_mumu = Counter(tree_Z_jets_mumu, hist_Zjets_mumu, directory);
    // cout << "   Zjets_total = " << n_Zjets_ee[0] + n_Zjets_mumu[0] << "+-" << sqrt(pow(n_Zjets_ee[1],2) + pow(n_Zjets_mumu[1],2)) <<  endl << endl;

    cout << "   Z_jets_ee:";
    vector<Float_t> n_Zjets_ee = ZCounter(tree_Z_jets_ee, hist_Zjets_ee0, hist_Zjets_ee1, hist_Zjets_ee2, directory);
    cout << "    ee_total = " << n_Zjets_ee[6] << "+-" << n_Zjets_ee[7] << endl << endl;
    cout << "   Z_jets_mumu:";
    vector<Float_t> n_Zjets_mumu = ZCounter(tree_Z_jets_mumu, hist_Zjets_mumu0, hist_Zjets_mumu1, hist_Zjets_mumu2, directory);
    cout << "   mumu_total = " << n_Zjets_mumu[6] << "+-" << n_Zjets_mumu[7] << endl << endl;
    cout << "   Zjets0_total = " << n_Zjets_ee[0] + n_Zjets_mumu[0] << "+-" << sqrt(pow(n_Zjets_ee[1],2) + pow(n_Zjets_mumu[1],2)) <<  endl << endl;
    cout << "   Zjets1_total = " << n_Zjets_ee[2] + n_Zjets_mumu[2] << "+-" << sqrt(pow(n_Zjets_ee[3],2) + pow(n_Zjets_mumu[3],2)) <<  endl << endl;
    cout << "   Zjets2_total = " << n_Zjets_ee[4] + n_Zjets_mumu[4] << "+-" << sqrt(pow(n_Zjets_ee[5],2) + pow(n_Zjets_mumu[5],2)) <<  endl << endl;
    cout << "   Zjets_total =  " << n_Zjets_ee[6] + n_Zjets_mumu[6] << "+-" << sqrt(pow(n_Zjets_ee[7],2) + pow(n_Zjets_mumu[7],2)) <<  endl << endl;

    cout << "   ================== top ==================    " << endl << endl;
    cout << "   Top:";
    vector<Float_t> n_top = Counter(tree_top, hist_top, directory);
    cout << "   ttbarV_ttbarVV:";
    vector<Float_t> n_ttbarV_ttbarVV = Counter(tree_ttbarV_ttbarVV, hist_ttbarV_ttbarVV, directory);
    cout << "   Wt:";
    vector<Float_t> n_Wt = Counter(tree_Wt, hist_Wt, directory);

    cout << "   Top:  " << n_top[0] + n_ttbarV_ttbarVV[0] + n_Wt[0] << " +- " << sqrt(pow(n_top[1], 2) + pow(n_ttbarV_ttbarVV[1], 2) + pow(n_Wt[1], 2)) << endl << endl;

    cout << "   ================== WW ==================    " << endl << endl;
    cout << "   WW:";
    vector<Float_t> n_WW = Counter(tree_WW, hist_WW, directory);

    cout << "   ================== other ==================    " << endl << endl;
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

    cout << "   Other:  " << n_llll[0] + n_llqq[0] + n_VVV[0] + n_Wjets[0] + n_Ztt[0] + n_WZjj[0] + n_lllljj[0] + n_llvvjj_WW[0] << " +- "
    << sqrt(pow(n_llll[1], 2) + pow(n_llqq[1], 2) + pow(n_VVV[1], 2) + pow(n_Wjets[1], 2) + pow(n_Ztt[1], 2) + pow(n_WZjj[1], 2) + pow(n_lllljj[1], 2) + pow(n_llvvjj_WW[1], 2)) 
    << endl << endl;
    
    cout << "   ----------------------------------------------------------"  << endl << endl;

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

    // other = llll + llqq + VVV + Wjets + Ztt
    hist_other->SetFillColor(TColor::GetColor("#50C878")); // KHAKI
    hist_other->SetLineColor(kBlack);
    hist_other->SetLineWidth(1);
    hist_other->Add(hist_llll);
    hist_other->Add(hist_llqq);
    hist_other->Add(hist_VVV);
    hist_other->Add(hist_W_jets);
    hist_other->Add(hist_Ztt);
    hist_other->Add(hist_WZ_jj);
    hist_other->Add(hist_lllljj);
    hist_other->Add(hist_llvvjj_WW);

    //Dummy uncertainty hist to show on legend
    hist_uncert_dummy->SetFillColor(kBlack);
    hist_uncert_dummy->SetFillStyle(3004);
    hist_uncert_dummy->SetLineColor(kBlack);   // Set the outline color
    hist_uncert_dummy->SetLineWidth(1);      // Set the outline width


    Float_t events_total;
    Float_t events_total_er;
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
    Float_t signal_obs;
    Float_t signal_obs_er;
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
    Float_t events_other = n_llll[0] + n_llqq[0] + n_VVV[0] + n_Wjets[0] + n_Ztt[0] + n_WZjj[0] + n_lllljj[0] + n_llvvjj_WW[0];
    Float_t events_other_er = sqrt(pow(n_llll[1], 2) + pow(n_llqq[1], 2) + pow(n_VVV[1], 2) + pow(n_Wjets[1], 2) + pow(n_Ztt[1], 2) + pow(n_WZjj[1], 2) + pow(n_lllljj[1], 2) + pow(n_llvvjj_WW[1], 2));

    Double_t S;
    Double_t B;
    Double_t Z;

    if (directory == "SR")
    {

      B = events_top + events_WW + events_WZ + events_Zjets0 + events_Zjets1 + events_Zjets2 + events_other;
      S = signal_mc;

      if (B > 0 && S > 0)
      {
        Z = sqrt(2 * ((S + B) * log(1 + (S / B)) - S));
      }

      cout << "   Significance for signal mc:   " << Z << endl<< endl;
    }

    

    //Scaling factors calculations 

    if (directory == "emCR_B")
    {
      events_nontop = signal_mc + events_WZ + events_WW + events_Zjets0 + events_Zjets1 + events_Zjets2 + events_other;
      events_nontop_er = sqrt(pow(signal_mc_er, 2) + pow(events_WZ_er, 2) + pow(events_WW_er, 2) + pow(events_Zjets0_er, 2) + pow(events_Zjets1_er, 2) + pow(events_Zjets2_er, 2) + pow(events_other_er, 2));
      sf_emuB = (events_data - events_nontop) / events_top;
      sf_emuB_er = sqrt(pow((events_data - events_nontop), 2) * pow(events_top_er, 2)/pow(events_top, 4) + (pow(events_data_er, 2) + pow(events_nontop_er, 2))/pow(events_top, 2));
      
      cout << "   Events Total MC =  " << events_top + events_nontop << endl << endl;
      cout << "   Events Total NR =  " << events_top + events_WW << endl << endl;
      cout << "   TOP SCALING FACTOR =  " << sf_emuB << " +- " << sf_emuB_er << endl << endl;
      cout << "   Events NR / Total (%) =  " << (events_top + events_WW) * 100 / (events_top + events_nontop) << endl << endl;
      

    }
    else if (directory == "emCR_A")
    {

      //Correcting for top events
      events_top = events_top * sf_emuB;
      events_top_er = sqrt( pow(sf_emuB, 2) * pow(events_top_er, 2) + pow(events_top, 2) * pow(sf_emuB_er, 2));


      events_nonWW = signal_mc + events_WZ + events_top + events_Zjets0 + events_Zjets1 + events_Zjets2 + events_other;
      events_nonWW_er = sqrt(pow(signal_mc_er, 2) + pow(events_WZ_er, 2) + pow(events_top_er, 2) + pow(events_Zjets0_er, 2) + pow(events_Zjets1_er, 2) + pow(events_Zjets2_er, 2) + pow(events_other_er, 2));
      sf_emuA = (events_data - events_nonWW) / events_WW;
      sf_emuA_er = sqrt(pow((events_data - events_nonWW), 2) * pow(events_WW_er, 2)/pow(events_WW, 4) + (pow(events_data_er, 2) + pow(events_nonWW_er, 2))/pow(events_WW, 2));

      cout << "   Events Total MC =  " << events_WW + events_nonWW << endl << endl;
      cout << "   Events Total NR =  " << events_top + events_WW << endl << endl;
      cout << "   Events NR / Total (%) =  " << (events_top + events_WW) * 100 / (events_WW + events_nonWW) << endl << endl;
      cout << "   WW SCALING FACTOR =  " << sf_emuA << " +- " << sf_emuA_er << endl << endl;
      

    }
    else if (directory == "3lCR")
    {
      //Correcting top and WW events
      events_top = events_top * sf_emuB;
      events_top_er = sqrt(pow(sf_emuB, 2) * pow(events_top_er, 2) + pow(events_top, 2) * pow(sf_emuB_er, 2));
      events_WW = events_WW * sf_emuA;
      events_WW_er = sqrt( pow(sf_emuA, 2) * pow(events_WW_er, 2) + pow(events_WW, 2) * pow(sf_emuA_er, 2));


      events_nonWZ = signal_mc + events_WW + events_top + events_Zjets0 + events_Zjets1 + events_Zjets2 + events_other;
      events_nonWZ_er = sqrt(pow(signal_mc_er, 2) + pow(events_WW_er, 2) + pow(events_top_er, 2) + pow(events_Zjets0_er, 2) + pow(events_Zjets1_er, 2) + pow(events_Zjets2_er, 2) + pow(events_other_er, 2));
      sf_3lCR = (events_data - events_nonWZ) / events_WZ; 
      sf_3lCR_er = sqrt(pow((events_data - events_nonWZ), 2) * pow(events_WZ_er, 2)/pow(events_WZ, 4)  +  (pow(events_data_er, 2) + pow(events_nonWZ_er, 2))/pow(events_WZ, 2));
      
      cout << "   Events Total MC =  " << events_WZ + events_nonWZ << endl << endl;
      cout << "   WZ SCALING FACTOR =  " << sf_3lCR << " +- " << sf_3lCR_er << endl << endl;
      cout << "   Events WZ / Total (%) =  " << events_WZ * 100 / (events_WZ + events_nonWZ) << endl << endl;
      
    
    }
    else if (directory == "Zjets2")
    {
      // Correcting for top, WW and WZ events
      events_top = events_top * sf_emuB;
      events_top_er = sqrt(pow(sf_emuB, 2) * pow(events_top_er, 2) + pow(events_top, 2) * pow(sf_emuB_er, 2));
      events_WW = events_WW * sf_emuA;
      events_WW_er = sqrt( pow(sf_emuA, 2) * pow(events_WW_er, 2) + pow(events_WW, 2) * pow(sf_emuA_er, 2));
      events_WZ = events_WZ*sf_3lCR;
      events_WZ_er = sqrt( pow(sf_3lCR, 2) * pow(events_WZ_er, 2) + pow(events_WZ, 2) * pow(sf_3lCR_er, 2));

      events_nonZjets2 = signal_mc + events_WZ + events_top + events_WW + events_Zjets0 + events_Zjets1 + events_other;
      events_nonZjets2_er = sqrt(pow(signal_mc_er, 2) + pow(events_WZ_er, 2) + pow(events_top_er, 2) + pow(events_WW_er, 2) + pow(events_Zjets0_er, 2) + pow(events_Zjets1_er, 2) + pow(events_other_er, 2));
      sf_Zjets2 = (events_data - events_nonZjets2) / events_Zjets2;
      sf_Zjets2_er = sqrt(pow((events_data - events_nonZjets2), 2) * pow(events_Zjets2_er, 2)/pow(events_Zjets2, 4) + (pow(events_data_er, 2) + pow(events_nonZjets2_er, 2))/pow(events_Zjets2, 2));
      
      cout << "   Events Total MC =  " << events_Zjets2 + events_nonZjets2 << endl << endl;
      cout << "   Zjets2 SCALING FACTOR =  " << sf_Zjets2 << " +- " << sf_Zjets2_er << endl << endl;
      cout << "   Events Zjets2 / Total (%) =  " << events_Zjets2 * 100/ (events_Zjets2 + events_nonZjets2) << endl << endl;
      
    }
    else if (directory == "Zjets1")
    {

      //Correcting for top, WW, WZ and Zjets2 events
      events_top = events_top * sf_emuB;
      events_top_er = sqrt(pow(sf_emuB, 2) * pow(events_top_er, 2) + pow(events_top, 2) * pow(sf_emuB_er, 2));
      events_WW = events_WW * sf_emuA;
      events_WW_er = sqrt( pow(sf_emuA, 2) * pow(events_WW_er, 2) + pow(events_WW, 2) * pow(sf_emuA_er, 2));
      events_WZ = events_WZ*sf_3lCR;
      events_WZ_er = sqrt( pow(sf_3lCR, 2) * pow(events_WZ_er, 2) + pow(events_WZ, 2) * pow(sf_3lCR_er, 2));
      events_Zjets2 = events_Zjets2 * sf_Zjets2;
      events_Zjets2_er = sqrt( pow(sf_Zjets2, 2) * pow(events_Zjets2_er, 2) + pow(events_Zjets2, 2) * pow(sf_Zjets2_er, 2));
    

      events_nonZjets1 = signal_mc + events_WZ + events_top + events_WW + events_Zjets0 + events_Zjets2 + events_other;
      events_nonZjets1_er = sqrt(pow(signal_mc_er, 2) + pow(events_WZ_er, 2) + pow(events_top_er, 2) + pow(events_WW_er, 2) + pow(events_Zjets0_er, 2) + pow(events_Zjets2_er, 2) + pow(events_other_er, 2));
      sf_Zjets1 = (events_data - events_nonZjets1) / events_Zjets1;
      sf_Zjets1_er = sqrt(pow((events_data - events_nonZjets1), 2) * pow(events_Zjets1_er, 2)/pow(events_Zjets1, 4) + (pow(events_data_er, 2) + pow(events_nonZjets1_er, 2))/pow(events_Zjets1, 2));
      
      cout << "   Events Total MC =  " << events_Zjets1 + events_nonZjets1 << endl << endl;
      cout << "   Zjets1 SCALING FACTOR =  " << sf_Zjets1 << " +- " << sf_Zjets1_er << endl << endl;
      cout << "   Events Zjets1 / Total (%) =  " << events_Zjets1 * 100 / (events_Zjets1 + events_nonZjets1) << endl << endl;
      
    }
    else if (directory == "Zjets0")
    {

      //Correcting for top, WW, WZ, Zjets2 and Zjets1 events
      events_top = events_top * sf_emuB;
      events_top_er = sqrt(pow(sf_emuB, 2) * pow(events_top_er, 2) + pow(events_top, 2) * pow(sf_emuB_er, 2));
      events_WW = events_WW * sf_emuA;
      events_WW_er = sqrt( pow(sf_emuA, 2) * pow(events_WW_er, 2) + pow(events_WW, 2) * pow(sf_emuA_er, 2));
      events_WZ = events_WZ*sf_3lCR;
      events_WZ_er = sqrt( pow(sf_3lCR, 2) * pow(events_WZ_er, 2) + pow(events_WZ, 2) * pow(sf_3lCR_er, 2));
      events_Zjets2 = events_Zjets2 * sf_Zjets2;
      events_Zjets2_er = sqrt( pow(sf_Zjets2, 2) * pow(events_Zjets2_er, 2) + pow(events_Zjets2, 2) * pow(sf_Zjets2_er, 2));
      events_Zjets1 = events_Zjets1 * sf_Zjets1;
      events_Zjets1_er = sqrt( pow(sf_Zjets1, 2) * pow(events_Zjets1_er, 2) + pow(events_Zjets1, 2) * pow(sf_Zjets1_er, 2));

      events_nonZjets0 = signal_mc + events_WZ + events_top + events_WW + events_Zjets1 + events_Zjets2 + events_other;
      events_nonZjets0_er = sqrt(pow(signal_mc_er, 2) + pow(events_WZ_er, 2) + pow(events_top_er, 2) + pow(events_WW_er, 2) + pow(events_Zjets1_er, 2) + pow(events_Zjets2_er, 2)   + pow(events_other_er, 2));
      sf_Zjets0 = (events_data - events_nonZjets0) / events_Zjets0;
      sf_Zjets0_er = sqrt(pow((events_data - events_nonZjets0), 2) * pow(events_Zjets0_er, 2)/pow(events_Zjets0, 4) + (pow(events_data_er, 2) + pow(events_nonZjets0_er, 2))/pow(events_Zjets0, 2));
 
      events_total = events_Zjets0 + events_nonZjets0;
      events_total_er =sqrt(pow(events_Zjets0_er, 2) + pow(events_nonZjets0_er, 2));
 
      cout << "   Events Total MC =  " << events_total << endl << endl;
      cout << "   Zjets0 SCALING FACTOR =  " << sf_Zjets0 << " +- " << sf_Zjets0_er << endl << endl;
      cout << "   Events Zjets0 / Total (%) =  " << events_Zjets0 * 100/ events_total << " +- " << 100 * sqrt(pow(events_Zjets0_er / events_total, 2) + pow(events_Zjets0 * events_total_er / pow(events_total, 2), 2)) << endl << endl;
      
    }
    else if (directory == "SR")
    {
      //Correcting for top, WW, WZ, Zjets2, Zjets1 and Zjets0 events
      events_top = events_top * sf_emuB;
      events_top_er = sqrt( pow(sf_emuB, 2) * pow(events_top_er, 2) + pow(events_top, 2) * pow(sf_emuB_er, 2));
      events_WW = events_WW * sf_emuA;
      events_WW_er = sqrt( pow(sf_emuA, 2) * pow(events_WW_er, 2) + pow(events_WW, 2) * pow(sf_emuA_er, 2));
      events_WZ = events_WZ*sf_3lCR;
      events_WZ_er = sqrt( pow(sf_3lCR, 2) * pow(events_WZ_er, 2) + pow(events_WZ, 2) * pow(sf_3lCR_er, 2));
      events_Zjets2 = events_Zjets2 * sf_Zjets2;
      events_Zjets2_er = sqrt( pow(sf_Zjets2, 2) * pow(events_Zjets2_er, 2) + pow(events_Zjets2, 2) * pow(sf_Zjets2_er, 2));
      events_Zjets1 = events_Zjets1 * sf_Zjets1;
      events_Zjets1_er = sqrt( pow(sf_Zjets1, 2) * pow(events_Zjets1_er, 2) + pow(events_Zjets1, 2) * pow(sf_Zjets1_er, 2));
      events_Zjets0 = events_Zjets0 * sf_Zjets0;
      events_Zjets0_er = sqrt( pow(sf_Zjets0, 2) * pow(events_Zjets0_er, 2) + pow(events_Zjets0, 2) * pow(sf_Zjets0_er, 2));
      
      //Signal
      events_bkg = events_WZ + events_top + events_WW + events_Zjets0 + events_Zjets1 + events_Zjets2 + events_other;
      events_bkg_er = sqrt(pow(events_WZ_er, 2) + pow(events_top_er, 2) + pow(events_WW_er, 2) + pow(events_Zjets0_er, 2) + pow(events_Zjets1_er, 2) + pow(events_Zjets2_er, 2) + pow(events_other_er, 2));
      
      // signal_obs = events_data - events_bkg;
      // signal_obs_er = sqrt(pow(events_data_er, 2) + pow(events_bkg_er, 2));

      events_data = signal_mc + events_bkg;
      // events_data_er = sqrt(events_data);
      events_data_er = sqrt(signal_mc + events_bkg);

      signal_obs = events_data - events_bkg;
      signal_obs_er = sqrt(pow(events_data_er, 2) + pow(events_bkg_er, 2));

      // signal_obs = (signal_mc + events_bkg) - events_bkg;
      // signal_obs_er = signal_mc_er;

      sf_signal = signal_obs / signal_mc;
      sf_signal_er = sqrt(pow(signal_obs, 2) * pow(signal_mc_er, 2)/pow(signal_mc, 4) + pow(signal_obs_er, 2)/pow(signal_mc, 2));


      cout << "   Events Total MC =  " << events_bkg + signal_mc << " +- " << sqrt(events_bkg_er * events_bkg_er + signal_mc_er * signal_mc_er) << endl << endl;
    }
    else if (directory == "VR")
    {
      // // MLE Simfit factors
      // sf_3lCR = 1.012;
      // sf_emuB = 1.013;
      // sf_emuA = 1.435;
      // sf_Zjets0 = 1.352;
      // sf_Zjets1 = 1.322;
      // sf_Zjets2 = 1.065;
      // sf_signal = 1.006;

      //Correcting for top, WW, WZ, Zjets2, Zjets1, Zjets0 and signal events
      events_top = events_top * sf_emuB;
      events_top_er = sqrt( pow(sf_emuB, 2) * pow(events_top_er, 2) + pow(events_top, 2) * pow(sf_emuB_er, 2));
      events_WW = events_WW * sf_emuA;
      events_WW_er = sqrt( pow(sf_emuA, 2) * pow(events_WW_er, 2) + pow(events_WW, 2) * pow(sf_emuA_er, 2));
      events_WZ = events_WZ*sf_3lCR;
      events_WZ_er = sqrt( pow(sf_3lCR, 2) * pow(events_WZ_er, 2) + pow(events_WZ, 2) * pow(sf_3lCR_er, 2));
      events_Zjets2 = events_Zjets2 * sf_Zjets2;
      events_Zjets2_er = sqrt( pow(sf_Zjets2, 2) * pow(events_Zjets2_er, 2) + pow(events_Zjets2, 2) * pow(sf_Zjets2_er, 2));
      events_Zjets1 = events_Zjets1 * sf_Zjets1;
      events_Zjets1_er = sqrt( pow(sf_Zjets1, 2) * pow(events_Zjets1_er, 2) + pow(events_Zjets1, 2) * pow(sf_Zjets1_er, 2));
      events_Zjets0 = events_Zjets0 * sf_Zjets0;
      events_Zjets0_er = sqrt( pow(sf_Zjets0, 2) * pow(events_Zjets0_er, 2) + pow(events_Zjets0, 2) * pow(sf_Zjets0_er, 2));
      signal_mc = signal_mc * sf_signal;
      signal_mc_er = sqrt( pow(sf_signal, 2) * pow(signal_mc_er, 2) + pow(signal_mc, 2) * pow(sf_signal_er, 2));

      
      //Signal
      events_bkg = events_WZ + events_top + events_WW + events_Zjets0 + events_Zjets1 + events_Zjets2 + events_other;
      events_bkg_er = sqrt(pow(events_WZ_er, 2) + pow(events_top_er, 2) + pow(events_WW_er, 2) + pow(events_Zjets0_er, 2) + pow(events_Zjets1_er, 2) + pow(events_Zjets2_er, 2) + pow(events_other_er, 2));
      signal_obs = events_data - events_bkg;
      signal_obs_er = sqrt(pow(events_data_er, 2) + pow(events_bkg_er, 2));
      

    }

    events_total = signal_mc + events_WZ + events_top + events_WW + events_Zjets0 + events_Zjets1 + events_Zjets2 + events_other;
    events_total_er = sqrt(pow(signal_mc_er, 2) + pow(events_WZ_er, 2) + pow(events_top_er, 2) + pow(events_WW_er, 2) + pow(events_Zjets0_er, 2) + pow(events_Zjets1_er, 2) + pow(events_Zjets2_er, 2) + pow(events_other_er, 2));

    cout << "   Signal / Total (%) =  " << signal_mc * 100 / events_total << " +- " << 100 * sqrt(pow(signal_mc_er / events_total, 2) + pow(signal_mc * events_total_er / pow(events_total, 2), 2)) << endl << endl;
    cout << "   Signal / Data (%) =  " << signal_mc * 100 / events_data << " +- " << 100 * sqrt(pow(signal_mc_er / events_data, 2) + pow(signal_mc * events_data_er / pow(events_data, 2), 2)) << endl << endl;
    cout << "   Data / MC          =  " << events_data  / events_total << " +- " << sqrt(pow(events_data_er / events_total, 2) + pow(events_data * events_total_er / pow(events_total, 2), 2)) << endl << endl;



    //prefit - postfit

    // MLE Simfit factors
    sf_3lCR = 1.012;
    sf_emuB = 1.013;
    sf_emuA = 1.435;
    sf_Zjets0 = 1.352;
    sf_Zjets1 = 1.322;
    sf_Zjets2 = 1.065;
    sf_signal = 1.006;

    // hist_top->SetBinContent(bin, hist_top->GetBinContent(bin) * sf_emuB);
    // hist_WW->SetBinContent(bin, hist_WW->GetBinContent(bin) * sf_emuA);
    // hist_WZ->SetBinContent(bin, hist_WZ->GetBinContent(bin) * sf_3lCR);
    // hist_Zjets2->SetBinContent(bin, hist_Zjets2->GetBinContent(bin) * sf_Zjets2);
    // hist_Zjets1->SetBinContent(bin, hist_Zjets1->GetBinContent(bin) * sf_Zjets1);
    // hist_Zjets0->SetBinContent(bin, hist_Zjets0->GetBinContent(bin) * sf_Zjets0);
    // hist_signal->SetBinContent(bin, hist_signal->GetBinContent(bin) * sf_signal);

    for (int bin = 1; bin < sizeof(xbins) / sizeof(xbins[0]); bin++)
    {
      if (directory == "emCR_A")
      {
        hist_top->SetBinContent(bin, hist_top->GetBinContent(bin) * sf_emuB);
      }
      if (directory == "3lCR")
      {
        hist_top->SetBinContent(bin, hist_top->GetBinContent(bin) * sf_emuB);
        hist_WW->SetBinContent(bin, hist_WW->GetBinContent(bin) * sf_emuA);
      }
      else if (directory == "Zjets2")
      {
        hist_top->SetBinContent(bin, hist_top->GetBinContent(bin) * sf_emuB);
        hist_WW->SetBinContent(bin, hist_WW->GetBinContent(bin) * sf_emuA);
        hist_WZ->SetBinContent(bin, hist_WZ->GetBinContent(bin) * sf_3lCR);
      }
      else if (directory == "Zjets1")
      {
        hist_top->SetBinContent(bin, hist_top->GetBinContent(bin) * sf_emuB);
        hist_WW->SetBinContent(bin, hist_WW->GetBinContent(bin) * sf_emuA);
        hist_WZ->SetBinContent(bin, hist_WZ->GetBinContent(bin) * sf_3lCR);
        hist_Zjets2->SetBinContent(bin, hist_Zjets2->GetBinContent(bin) * sf_Zjets2);        
      }
      else if (directory == "Zjets0")
      {
        hist_top->SetBinContent(bin, hist_top->GetBinContent(bin) * sf_emuB);
        hist_WW->SetBinContent(bin, hist_WW->GetBinContent(bin) * sf_emuA);
        hist_WZ->SetBinContent(bin, hist_WZ->GetBinContent(bin) * sf_3lCR);
        hist_Zjets2->SetBinContent(bin, hist_Zjets2->GetBinContent(bin) * sf_Zjets2);        
        hist_Zjets1->SetBinContent(bin, hist_Zjets1->GetBinContent(bin) * sf_Zjets1);
      }
      else if (directory == "SR")
      {
        hist_top->SetBinContent(bin, hist_top->GetBinContent(bin) * sf_emuB);
        hist_WW->SetBinContent(bin, hist_WW->GetBinContent(bin) * sf_emuA);
        hist_WZ->SetBinContent(bin, hist_WZ->GetBinContent(bin) * sf_3lCR);
        hist_Zjets2->SetBinContent(bin, hist_Zjets2->GetBinContent(bin) * sf_Zjets2);        
        hist_Zjets1->SetBinContent(bin, hist_Zjets1->GetBinContent(bin) * sf_Zjets1);
        hist_Zjets0->SetBinContent(bin, hist_Zjets0->GetBinContent(bin) * sf_Zjets0);
      }
    }

    //Merge Zjets before start plotting
    hist_Zjets->Add(hist_Zjets0);
    hist_Zjets->Add(hist_Zjets1);
    hist_Zjets->Add(hist_Zjets2);



    cout << "   Significance mc after scaling is:   " << Z << endl << endl;
    cout << "   RECO INTEGRAL IS:   " << hist_signal_obs->Integral(1, hist_signal_obs->GetNbinsX()) << endl << endl;
    //Print calculated events for every region
  
    cout << "------------------------------------------------------------------" << endl << endl;
    cout << "   SIGNAL MC       =  " << signal_mc << " +- " << signal_mc_er  << endl;
    cout << "   SIGNAL RECO     =  " << signal_obs << " +- " << signal_obs_er  << endl;
    cout << "   SIGNAL RECO/BKG =  " << signal_obs/events_bkg << endl << endl;
    cout << "_________________________________" << endl << endl;
    // cout << "   Signal: " << "     " << signal_mc << " +- " << signal_mc_er << endl << endl;
    cout << "   Data: " << "       " << events_data   << " +- " << events_data_er   << endl << endl;
    cout << "   WZ: " << "         " << events_WZ     << " +- " << events_WZ_er     << "   |   μ_WZ = " << sf_3lCR << " +- " << sf_3lCR_er << endl << endl;
    cout << "   Top: " << "        " << events_top    << " +- " << events_top_er    << "   |   μ_top = " << sf_emuB << " +- " << sf_emuB_er << endl << endl;
    cout << "   WW: " << "         " << events_WW     << " +- " << events_WW_er     << "   |   μ_WW = " << sf_emuA << " +- " << sf_emuA_er << endl << endl;
    cout << "   Zjets0: " << "     " << events_Zjets0 << " +- " << events_Zjets0_er << "   |   μ_Zjets0 = " << sf_Zjets0 << " +- " << sf_Zjets0_er << endl << endl;
    cout << "   Zjets1: " << "     " << events_Zjets1 << " +- " << events_Zjets1_er << "   |   μ_Zjets1 = " << sf_Zjets1 << " +- " << sf_Zjets1_er << endl << endl;
    cout << "   Zjets2: " << "     " << events_Zjets2 << " +- " << events_Zjets2_er << "   |   μ_Zjets2 = " << sf_Zjets2 << " +- " << sf_Zjets2_er << endl << endl;
    cout << "   Signal: " << "     " << signal_mc     << " +- " << signal_mc_er     << "   |    μ_signal = " << sf_signal << " +- " << sf_signal_er << endl << endl;
    cout << "   Other: "  << "     " << events_other  << " +- " << events_other_er  << endl << endl;
    cout << "------------------------------------------------------------------"    << endl << endl;

    // Stacking with a specific order
    hist_Zjets->Add(hist_other);
    hist_top->Add(hist_Zjets);
    hist_WW->Add(hist_top);
    hist_WZ->Add(hist_WW);
    hist_signal->Add(hist_WZ);    
    cout << "   DATA/MC Scaled = " << events_data / (events_bkg + signal_mc) << endl << endl;
    

    //----------------------------------------PLOTS----------------------------------------//


    TCanvas *c1 = new TCanvas("c1", "SR", 1400, 600, 1000, 1050);
    TCanvas *c2 = new TCanvas("c2", "3lCR", 1400, 600, 1000, 1050);
    TCanvas *c3 = new TCanvas("c3", "emuB", 1400, 600, 1000, 1050);
    TCanvas *c4 = new TCanvas("c4", "emuA", 1400, 600, 1000, 1050);
    TCanvas *c5 = new TCanvas("c5", "Zjets0", 1400, 600, 1000, 1050);
    TCanvas *c6 = new TCanvas("c6", "Zjets1", 1400, 600, 1000, 1050);
    TCanvas *c7 = new TCanvas("c7", "Zjets2", 1400, 600, 1000, 1050);
    TCanvas *c8 = new TCanvas("c8", "VR", 1400, 600, 1000, 1050);

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
    else if (directory == "VR")
    {
      c8->cd();
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
    gStyle->SetLegendFont(102);
    gStyle->SetTextFont(62);
    gStyle->SetPalette(20);


    //---------------PAD 1---------------//

    // Plotting according to stacking order
    pad1->cd();

    // hist_signal->Draw("hist");
    // hist_WZ->Draw("histsame");
    // hist_WW->Draw("histsame");
    // hist_top->Draw("histsame");
    // hist_Zjets->Draw("histsame");
    // hist_other->Draw("histsame");
    // hist_data->Draw("same");

    // hist_data->SetBinContent(3, -10);
    // hist_signal->SetBinContent(3, -10);
    // hist_WZ->SetBinContent(3, 0);
    // hist_WW->SetBinContent(3, 0);
    // hist_top->SetBinContent(3, 0);
    // hist_Zjets->SetBinContent(3, 0);
    // hist_other->SetBinContent(3, 0);



    hist_signal->Draw("hist");

    if (directory == "SR")
    {
      // For blinded data in SR
      hist_data = static_cast<TH1F *>(hist_signal->Clone("hist_data"));
      hist_data->SetMarkerStyle(20);
      hist_data->SetMarkerSize(1.);
      hist_data->SetMarkerColor(kBlack);
      hist_data->SetLineColor(kBlack);
      
    }

    hist_WZ->Draw("histsame");
    hist_WW->Draw("histsame");
    hist_top->Draw("histsame");
    hist_Zjets->Draw("histsame");
    hist_other->Draw("histsame");
    hist_data->Draw("same");
    

    for (int bin = 1; bin <= hist_signal->GetNbinsX(); ++bin)
    {

      float bin_error = hist_signal->GetBinError(bin);

      float xmin = hist_signal->GetBinLowEdge(bin);
      float xmax = hist_signal->GetBinLowEdge(bin + 1);

      float ymin = hist_signal->GetBinContent(bin) - bin_error;
      float ymax = hist_signal->GetBinContent(bin) + bin_error;
      


      TBox *box = new TBox(xmin, ymin, xmax, ymax);
      // box->SetFillStyle(3013);
      box->SetFillStyle(3004);
      box->SetLineColor(kBlack);
      // box->SetFillColor(kBlack);
      box->SetFillColorAlpha(kBlack, 1.);
      box->Draw("f");
      }

      // //Normalise to unity
      // hist_norm(hist_signal);
      // hist_norm(hist_WZ);
      // hist_norm(hist_WW);
      // hist_norm(hist_top);
      // hist_norm(hist_Zjets);
      // hist_norm(hist_other);
      // hist_norm(hist_data);

      hist_signal->GetXaxis()->SetTitleSize(0.05);
      hist_signal->GetXaxis()->SetTitleFont(42);
      hist_signal->GetXaxis()->SetTitleOffset(1.1);
      hist_signal->GetXaxis()->SetLabelSize(0);

      hist_signal->GetYaxis()->SetTitleFont(62);
      hist_signal->GetYaxis()->SetTitleOffset(1.4);
      hist_signal->GetYaxis()->SetTitleSize(0.04);
      hist_signal->GetYaxis()->SetTitle("Events");
      hist_signal->GetYaxis()->SetRangeUser(0, hist_data->GetMaximum() * 1.5);
      hist_signal->SetStats(0);
      


    Int_t right_corner = 2; //For use to choose corner to plot information
    plot_info(directory, right_corner, fitFunc, hist_data, hist_signal, hist_signal_obs, hist_WZ, hist_WW, hist_Zjets, hist_top, hist_other, hist_uncert_dummy);
    pad1->SetTicks();
    pad1->RedrawAxis();



    //---------------PAD 2---------------//

    pad2->cd();

    TH1F *numerator = new TH1F("Numerator", " ", sizeof(xbins) / sizeof(xbins[0]) - 1, xbins);
    TH1F *denominator = new TH1F("Numerator", " ", sizeof(xbins) / sizeof(xbins[0]) - 1, xbins);
    
    
    numerator->Add(hist_data);
    denominator->Add(hist_signal);
    // pad2->SetLogy();

    // numerator->Sumw2(1);
    // denominator->Sumw2(1);
    TLine *line = new TLine(xbins[0], 1, xbins[sizeof(xbins) / sizeof(xbins[0]) - 1], 1);
    line->SetLineStyle(2);
    line->SetLineWidth(2);
    
    numerator->SetMarkerStyle(20);
    numerator->GetXaxis()->SetTitle("");
    numerator->GetXaxis()->SetTitleSize(0.15);
    numerator->GetXaxis()->SetTitleFont(62);
    // numerator->GetXaxis()->SetTitle("#sigma(E^{miss}_{T})")<;
    // numerator->GetXaxis()->SetTitle("n_{jets}");
    // numerator->GetXaxis()->SetTitle("\\Delta R_{\\ell \\ell}");
    numerator->GetXaxis()->SetTitle("\\Delta \\phi (p_{T}^{Z}, E^{miss}_{T}) [rad]");
    // numerator->GetXaxis()->SetTitle("p_{T}^{l2} [GeV]");
    // numerator->GetXaxis()->SetTitle("p_{T}^{l1} [GeV]");
    // numerator->GetXaxis()->SetTitle("E^{miss}_{T} [GeV]");
    // numerator->GetXaxis()->SetTitle("p^{Z}_{T} [GeV]");
    // numerator->GetXaxis()->SetTitle("#eta^{Z}");
    // numerator->GetXaxis()->SetTitle("E^{miss}_{T}/H_{T}");
    numerator->GetXaxis()->SetTitleOffset(1.1);
    numerator->GetXaxis()->SetLabelSize(0.14);

    numerator->GetYaxis()->SetTitleOffset(0.4);
    numerator->GetYaxis()->SetTitleSize(0.12);
    numerator->GetYaxis()->SetNdivisions(505);
    numerator->GetYaxis()->SetRangeUser(0.7, 1.4);
    numerator->GetYaxis()->SetTitleFont(62);
    numerator->GetYaxis()->SetTitle("#frac{Data}{Prediction}");
    numerator->GetYaxis()->SetLabelSize(0.10);
    numerator->SetStats(0);

    numerator->SetTitle(" ");
    numerator->SetLineColor(kBlack);
    numerator->Divide(denominator);
    // numerator->SetBinContent(3, -10);

    //FIT FUNCTION

    TGraphErrors *graph = new TGraphErrors(numerator->GetNbinsX());
    for (int bin = 0; bin < numerator->GetNbinsX(); bin++)
    {
      double x = numerator->GetXaxis()->GetBinCenter(bin + 1);
      double y = numerator->GetBinContent(bin + 1);
      double ey = numerator->GetBinError(bin + 1);
      graph->SetPoint(bin, x, y);
      graph->SetPointError(bin, 0, ey); // Assuming the x-errors are negligible
    }

    
    fitFunc->SetRange(numerator->GetBinCenter(0), numerator->GetBinCenter(numerator->GetNbinsX()+1));
    fitFunc->SetLineWidth(2);

    graph->Fit(fitFunc, "R");

    double p0 = fitFunc->GetParameter(0);
    double p1 = fitFunc->GetParameter(1);
    double chi2 = fitFunc->GetChisquare();
    double ndf = fitFunc->GetNDF();

    numerator->Draw();
    line->Draw("same");
    fitFunc->Draw("same");
    

    pad2->SetTicks();
    pad2->RedrawAxis();

    //x Right
    double x_line = 0.64; 
    double x_chi2 = 0.64; 
    

    // //x Left
    // double x_line = 0.125; 
    // double x_chi2 = 0.125; 

    //y
    double y_line = 0.92;
    double y_chi2 = 0.82;

    
    

    TLatex tex_line;
    tex_line.SetNDC();
    tex_line.SetTextSize(0.1); 
    tex_line.SetTextFont(102);  

    TLatex tex_chi2;
    tex_chi2.SetNDC();
    tex_chi2.SetTextSize(0.1); 
    tex_chi2.SetTextFont(102);   

    char text_line[100];
    char text_chi2[100];

    sprintf(text_chi2, "#chi^{2}/NDF = %.1f/%.0f", chi2, ndf);

    if (p1 >= 0)
    {
      sprintf(text_line, "y = %.3f+%.3fx", p0, p1);
    }
    else
    {
      sprintf(text_line, "y = %.3f%.3fx", p0, p1);
    }
  
   
    tex_line.DrawLatex(x_line, y_line, text_line);
    tex_chi2.DrawLatex(x_chi2, y_chi2, text_chi2);


    // for (int bin = 1; bin <= numerator->GetNbinsX(); ++bin)
    // {

    //   float bin_error = numerator->GetBinError(bin);

    //   cout << numerator->GetBinContent(bin) - bin_error << endl
    //        << endl;

    //   float xmin = numerator->GetBinLowEdge(bin);
    //   float xmax = numerator->GetBinLowEdge(bin + 1);

    //   float ymin = numerator->GetBinContent(bin) - bin_error;
    //   float ymax = numerator->GetBinContent(bin) + bin_error;

    //   TBox *box = new TBox(xmin, ymin, xmax, ymax);
    //   // box->SetFillStyle(3013);
    //   box->SetFillStyle(3004);
    //   box->SetLineColor(kBlack);
    //   // box->SetFillColor(kBlack);
    //   box->SetFillColorAlpha(kBlack, 1.);
    //   box->Draw("f");
    // }

    

    if (directory == "SR")
    {
      c1->SaveAs("./zjets_splitted_sc/dmetzphi_SR-direct.png");
    }
    else if (directory == "3lCR")
    {
      c2->SaveAs("./zjets_splitted_sc/dmetzphi_3lCR-direct.png");
    }
    else if (directory == "emCR_B")
    {
      c3->SaveAs("./zjets_splitted_sc/dmetzphi_emuB-direct.png");
    }
    else if (directory == "emCR_A")
    {
      c4->SaveAs("./zjets_splitted_sc/dmetzphi_emuA-direct.png");
    }
    else if (directory == "Zjets0")
    {
      c5->SaveAs("./zjets_splitted_sc/dmetzphi_Zjets0-direct.png");
    }
    else if (directory == "Zjets1")
    {
      c6->SaveAs("./zjets_splitted_sc/dmetzphi_Zjets1-direct.png");
    }
    else if (directory == "Zjets2")
    {
      c7->SaveAs("./zjets_splitted_sc/dmetzphi_Zjets2-direct.png");
    }
    else if (directory == "VR")
    {
      c8->SaveAs("./zjets_splitted_sc/dmetzphi_VR-direct.png");
    }


    //To avoid memory leak
    delete c1;
    delete c2;
    delete c3;
    delete c4;
    delete c5;
    delete c6;
    delete c7;
    delete c8;
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
  cout << "   SCALING FACTORS:       μ_top = " << sf_emuB  << " +- " << sf_emuB_er << endl << endl
       << "                          μ_WW = " << sf_emuA << " +- " << sf_emuA_er << endl << endl
       << "                          μ_WZ = " <<  sf_3lCR << " +- " << sf_3lCR_er << endl << endl  
       << "                          μ_Zjets2 = " << sf_Zjets2 << " +- " << sf_Zjets2_er << endl << endl
       << "                          μ_Zjets1 = " << sf_Zjets1 << " +- " << sf_Zjets1_er << endl << endl
       << "                          μ_Zjets0 = " << sf_Zjets0 << " +- " << sf_Zjets0_er << endl << endl     
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
