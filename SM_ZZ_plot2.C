// For more details see https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/ZZjj2l2vRun2#Data_and_MC_samples
// Script for plotting 15/06/20 Dimitrii Krasnopevtsev
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
#include <TLatex.h>
#include <TPad.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TLine.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
using namespace std;

// gStyle->SetLegendFont(42);
// gStyle->SetPalette(1);

// Method to fill hitograms (MET and MET_significance) with events properly weighted
//
std::vector<float> DoReco(TTree *tree, TH1F *Hmet_pt, TH1F *Hmet_pt_signif)
{

  TH1::SetDefaultSumw2(kTRUE);

  Int_t nentries = (Int_t)tree->GetEntries();

  Float_t M2Lep = 0.;      // Mass of 2 leptons
  Float_t met_tst = 0.;    // Missing ET
  Float_t met_signif = 0.; // Missing ET over the sum of all ETs existing
  Float_t dMetZPhi = 0.;
  Float_t frac_pT = 0.;
  Float_t MetOHT = 0.;
  Float_t dPhiJ100met = 0.;
  Float_t dLepR = 0.;
  Int_t n_bjets = 0.;
  Int_t n_jets = 0.;
  Float_t leading_pT_lepton = 0;
  Float_t subleading_pT_lepton = 0;
  Float_t detajj = 0;
  Float_t mjj = 0;
  Float_t leading_jet_pt = 0;
  Float_t second_jet_pt = 0;
  Int_t event_3CR = 0.;  // Event with 3 leptons
  Int_t event_type = 0.; // Electron or muon
  Float_t weight;

  std::vector<float> events;
  events.clear();

  tree->SetBranchAddress("M2Lep", &M2Lep);
  tree->SetBranchAddress("met_tst", &met_tst);
  tree->SetBranchAddress("met_signif", &met_signif);
  tree->SetBranchAddress("dMetZPhi", &dMetZPhi);
  tree->SetBranchAddress("frac_pT", &frac_pT);
  tree->SetBranchAddress("MetOHT", &MetOHT);
  tree->SetBranchAddress("dPhiJ100met", &dPhiJ100met);
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

  double signal = 0.;
  double signaler = 0.;
  double Norm = 1;

  // Loop over events
  for (Int_t i = 0; i < nentries; i++)
  {
    tree->GetEntry(i);

    // ----- Define signal region -----  Change as needed:
    //
    // inclusive selection:
    //---------------------
    // RECO cut-based:
    // if (event_3CR==0 && (event_type==0 || event_type==1) &&  leading_pT_lepton>30 && subleading_pT_lepton>20 && M2Lep>80 && M2Lep<100 && n_bjets<1 && dLepR<1.8 && dMetZPhi>2.7 && met_tst > 110 && MetOHT>0.65)
    // if (event_3CR==0 && (event_type==0 || event_type==1) &&  leading_pT_lepton>30 && subleading_pT_lepton>20 && M2Lep>80 && M2Lep<100 && n_bjets<1 && dLepR<1.8 && dMetZPhi>2.7 && met_tst > 110 && MetOHT>0.65 && met_signif>10)
    // if (event_3CR==0 && (event_type==0 || event_type==1) &&  leading_pT_lepton>30 && subleading_pT_lepton>20 && M2Lep>80 && M2Lep<100 && n_bjets<1 && dLepR<1.8 && dMetZPhi>2.7 && met_tst > 110 && MetOHT>0.65 && met_signif>0) // met_tst/sqrt(met_tst/MetOHT)>10)
    // FIDUCIAL
    // if (event_3CR==0 && (event_type==0 || event_type==1) &&  leading_pT_lepton>30 && subleading_pT_lepton>20 && M2Lep>80 && M2Lep<100 && n_bjets<1 && dLepR<1.8 && dMetZPhi>2.7 && met_tst > 110 && MetOHT>0.65)
    if (event_3CR == 0 && (event_type == 0 || event_type == 1) && leading_pT_lepton > 30 && subleading_pT_lepton > 20 && M2Lep > 76 && M2Lep < 116 && n_bjets < 1 && dLepR < 1.9 && dMetZPhi > 2.6 && met_tst > 90 && MetOHT > 0.6)
    // FID-BDT:
    // if (event_3CR==0 && (event_type==0 || event_type==1) &&  leading_pT_lepton>30 && subleading_pT_lepton>20 && M2Lep>76 && M2Lep<106 && n_bjets<1 && dLepR<2.2 && dMetZPhi>1.3 && met_tst > 70 && MetOHT>0.3 && met_signif>0) // met_tst/sqrt(met_tst/MetOHT)>10)
    // FID-DNN:
    // if (event_3CR==0 && (event_type==0 || event_type==1) &&  leading_pT_lepton>30 && subleading_pT_lepton>20 && M2Lep>76 && M2Lep<106 && n_bjets<1 && dLepR<2.2 && dMetZPhi>1.3 && met_tst > 90 && MetOHT>0.1 && met_signif>0) // met_tst/sqrt(met_tst/MetOHT)>10)
    //
    // VBS selection:
    // --------------
    // if (event_3CR==0 && (event_type==0 || event_type==1) && met_tst>120 &&  met_signif>0 && dLepR<1.8 && dMetZPhi>2.5 && n_bjets<1 && MetOHT>0.6 && n_jets > 1)
    // if (event_3CR==0 && (event_type==0 || event_type==1) && met_tst>90 && met_signif>12 && dLepR<1.5 && dMetZPhi>2.6 && n_bjets<1 && MetOHT>0.6 )
    //(e.g., here: not an event in the 3-lepton Control-Region, event_type==0 or 1, MissingET > 90 GeV)
    // if (event_3CR==0 && (event_type==0 || event_type==1) && met_tst>90)
    // Nominal:
    // if (event_3CR==0 && (event_type==0 || event_type==1) &&  leading_pT_lepton>30  && subleading_pT_lepton>20 && M2Lep>80 && M2Lep<100 && n_jets>1 && n_bjets<1 && leading_jet_pt>30 && second_jet_pt>30 && detajj>1  && mjj>100 && dMetZPhi>2.2 && dLepR<1.8 && met_tst>150 && MetOHT>0.4)
    // Nominal + METSignif :
    // if (event_3CR==0 && (event_type==0 || event_type==1) &&  leading_pT_lepton>30  && subleading_pT_lepton>20 && M2Lep>80 && M2Lep<100 && n_jets>1 && n_bjets<1 && leading_jet_pt>30 && second_jet_pt>30 && detajj>1  && mjj>100 && dMetZPhi>2.2 && dLepR<1.8 && met_tst>150 && MetOHT>0.4 && met_signif>10 )
    // Inclusive + di-jet variables:
    // if (event_3CR==0 && (event_type==0 || event_type==1) &&  leading_pT_lepton>30  && subleading_pT_lepton>20 && M2Lep>80 && M2Lep<100 && n_jets>1 && n_bjets<1 && leading_jet_pt>30 && second_jet_pt>30 && detajj>1  && mjj>100 && dMetZPhi>2.2 && dLepR<1.8 && met_tst>110 && MetOHT>0.65 && met_signif>0 )

    {
      signal = signal + weight;              // signal yield is sum of weights
      signaler = signaler + weight * weight; // keep Sum(weights^2) for calculating the error on the signal yield
      Hmet_pt->Fill(met_tst, weight);
      Hmet_pt_signif->Fill(met_signif, weight);
    }
  }

  cout << "N=" << signal << "+-" << sqrt(signaler) << endl; // signal yield and error on this yield

  events.push_back(signal);
  events.push_back(sqrt(signaler));
  return events;
}

// the main method
//
int SM_ZZ_plot2()
{

  gStyle->SetLegendFont(42);
  gStyle->SetPalette(1);

  // TString sample = "/Results/trial2.root";
  // ---convert TString to string
  // std::string sample_string(sample.Data());
  // ---convert string to TString
  // TString sample_TString(sample_string)

  // Preliminary minitrees for ZZllvvjj analysis
  TFile *file_data = new TFile("data/data_fullRun2_MET90.root");
  TTree *tree_data = file_data->Get<TTree>("tree_PFLOW");

  TFile *file_signalQCD = new TFile("data/mc16ade_QCD_EWKcor_ZZllvv_MET90.root"); // QCD signal with EWK corrections
  TTree *tree_signalQCD = file_signalQCD->Get<TTree>("tree_PFLOW");

  TFile *file_signalEWK = new TFile("data/mc16ade_EWK_ZZllvvjj_MET90.root"); // EW signal
  TTree *tree_signalEWK = file_signalEWK->Get<TTree>("tree_PFLOW");

  TFile *file_Zjets = new TFile("data/mc16ade_Zjets_MET90.root");
  TTree *tree_Zjets = file_Zjets->Get<TTree>("tree_PFLOW");

  TFile *file_WZ = new TFile("data/mc16ade_WZ_MET90.root");
  TTree *tree_WZ = file_WZ->Get<TTree>("tree_PFLOW");

  TFile *file_WW = new TFile("data/mc16ade_WW_Pow_MET90.root");
  TTree *tree_WW = file_WW->Get<TTree>("tree_PFLOW");

  TFile *file_Wt = new TFile("data/mc16ade_Wt_MET90.root"); // Single top EW process
  TTree *tree_Wt = file_Wt->Get<TTree>("tree_PFLOW");

  TFile *file_tt = new TFile("data/mc16ade_TOP_MET90.root");
  TTree *tree_tt = file_tt->Get<TTree>("tree_PFLOW");

  TFile *file_trib = new TFile("data/mc16ade_VVV_MET90.root");
  TTree *tree_trib = file_trib->Get<TTree>("tree_PFLOW");

  TFile *file_othr = new TFile("data/mc16ade_ttV_ttVV_MET90.root");
  TTree *tree_othr = file_othr->Get<TTree>("tree_PFLOW");

  TH1::SetDefaultSumw2(kTRUE);

  Float_t xbins[21] = {90, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 345, 375, 410, 450, 500, 580, 700, 1000}; // Unequal bins

  TH1F *Hmet_pt_sigQCD = new TH1F("Hmet_pt_sigQCD", "Hmet_pt_sigQCD", 20, xbins);
  TH1F *Hmet_pt_sigEWK = new TH1F("Hmet_pt_sigEWK", "Hmet_pt_sigEWK", 20, xbins);
  TH1F *Hmet_pt_data = new TH1F("Hmet_pt_data", "Hmet_pt_data ", 20, xbins);
  TH1F *Hmet_pt_Zjets = new TH1F("Hmet_pt_Zjets", "Hmet_pt_Zjets ", 20, xbins);
  TH1F *Hmet_pt_WZ = new TH1F("Hmet_pt_WZ", "Hmet_pt_WZ ", 20, xbins);
  TH1F *Hmet_pt_WW = new TH1F("Hmet_pt_WW", "Hmet_pt_WW ", 20, xbins);
  TH1F *Hmet_pt_Wt = new TH1F("Hmet_pt_Wt", "Hmet_pt_Wt ", 20, xbins);
  TH1F *Hmet_pt_tt = new TH1F("Hmet_pt_tt", "Hmet_pt_tt ", 20, xbins);
  TH1F *Hmet_pt_trib = new TH1F("Hmet_pt_trib", "Hmet_pt_trib ", 20, xbins);
  TH1F *Hmet_pt_othr = new TH1F("Hmet_pt_othr", "Hmet_pt_othr ", 20, xbins);
  TH1F *Hmet_pt_signal_sing = new TH1F("Hmet_pt_signal_sing", "Hmet_pt_signal_sing", 20, xbins);

  // Significance histograms

  TH1F *Hmet_pt_signif_sigQCD = new TH1F("Hmet_pt_signif_sigQCD", "Hmet_pt_signif_sigQCD", 35, 0, 35);
  TH1F *Hmet_pt_signif_sigEWK = new TH1F("Hmet_pt_signif_sigEWK", "Hmet_pt_signif_sigEWK", 35, 0, 35);
  TH1F *Hmet_pt_signif_data = new TH1F("Hmet_pt_signif_data ", "Hmet_pt_signif_data ", 35, 0, 35);
  TH1F *Hmet_pt_signif_Zjets = new TH1F("Hmet_pt_signif_Zjets ", "Hmet_pt_signif_Zjets ", 35, 0, 35);
  TH1F *Hmet_pt_signif_WZ = new TH1F("Hmet_pt_signif_WZ ", "Hmet_pt_signif_WZ ", 35, 0, 35);
  TH1F *Hmet_pt_signif_WW = new TH1F("Hmet_pt_signif_WW ", "Hmet_pt_signif_WW ", 35, 0, 35);
  TH1F *Hmet_pt_signif_Wt = new TH1F("Hmet_pt_signif_Wt ", "Hmet_pt_signif_Wt ", 35, 0, 35);
  TH1F *Hmet_pt_signif_tt = new TH1F("Hmet_pt_signif_tt ", "Hmet_pt_signif_tt ", 35, 0, 35);
  TH1F *Hmet_pt_signif_trib = new TH1F("Hmet_pt_signif_trib ", "Hmet_pt_signif_trib", 35, 0, 35);
  TH1F *Hmet_pt_signif_othr = new TH1F("Hmet_pt_signif_othr ", "Hmet_pt_signif_othr", 35, 0, 35);

  TH1F *Hmet_pt_signif_signal_sing = new TH1F("Hmet_pt_signif_signal_sing ", "Hmet_pt_signif_signal_sing", 35, 0, 35);
  TH1F *Hmet_pt_signif_signal_pur = new TH1F("Hmet_pt_signif_signal_pur ", "Hmet_pt_signif_signal_pur", 35, 0, 35);

  std::vector<float> events_data;
  std::vector<float> events_NsigQCD;
  std::vector<float> events_NsigEWK;
  std::vector<float> events_NZjets;
  std::vector<float> events_NWZ;
  std::vector<float> events_NWW;
  std::vector<float> events_NWt;
  std::vector<float> events_Ntt;
  std::vector<float> events_Ntrib;
  std::vector<float> events_Nothr;

  // Event Yields
  cout << "===Data===" << endl;
  events_data = DoReco(tree_data, Hmet_pt_data, Hmet_pt_signif_data, 0); // DoReco Fills the vectors
  cout << "===Signal QCD ZZ===" << endl;
  events_NsigQCD = DoReco(tree_signalQCD, Hmet_pt_sigQCD, Hmet_pt_signif_sigQCD, 1);
  cout << "===Signal EWK ZZ===" << endl;
  events_NsigEWK = DoReco(tree_signalEWK, Hmet_pt_sigEWK, Hmet_pt_signif_sigEWK, 1);
  cout << "===Background===" << endl;
  cout << "Zjets" << endl;
  events_NZjets = DoReco(tree_Zjets, Hmet_pt_Zjets, Hmet_pt_signif_Zjets, 1);
  cout << "WZ" << endl;
  events_NWZ = DoReco(tree_WZ, Hmet_pt_WZ, Hmet_pt_signif_WZ, 1);
  cout << "tt" << endl;
  events_Ntt = DoReco(tree_tt, Hmet_pt_tt, Hmet_pt_signif_tt, 1);
  cout << "WW" << endl;
  events_NWW = DoReco(tree_WW, Hmet_pt_WW, Hmet_pt_signif_WW, 1);
  cout << "Wt" << endl;
  events_NWt = DoReco(tree_Wt, Hmet_pt_Wt, Hmet_pt_signif_Wt, 1);
  cout << "trib" << endl;
  events_Ntrib = DoReco(tree_trib, Hmet_pt_trib, Hmet_pt_signif_trib, 1);
  cout << "othr" << endl;
  events_Nothr = DoReco(tree_othr, Hmet_pt_othr, Hmet_pt_signif_othr, 1);

  // Some useful outputs
  double totalBKG = events_NZjets.at(0) + events_Ntt.at(0) + events_NWt.at(0) + events_NWW.at(0) + events_NWZ.at(0) + events_Ntrib.at(0) + events_Nothr.at(0);
  double totalBKG_er = sqrt(events_NZjets.at(1) * events_NZjets.at(1) + events_Ntt.at(1) * events_Ntt.at(1) + events_NWt.at(1) * events_NWt.at(1) + events_NWW.at(1) * events_NWW.at(1) +
                            events_NWZ.at(1) * events_NWZ.at(1) + events_Ntrib.at(1) * events_Ntrib.at(1) + events_Nothr.at(1) * events_Nothr.at(1));
  cout << "total Bkg=" << totalBKG << "+-" << totalBKG_er << endl;
  cout << "S/B=" << (events_NsigQCD.at(0) + events_NsigEWK.at(0)) / totalBKG << endl;
  Double_t S = (events_NsigQCD.at(0) + events_NsigEWK.at(0));
  Double_t B = totalBKG;
  Double_t Z = sqrt(2 * ((S + B) * log(1 + (S / B)) - S));
  cout << " Significance = " << Z << endl;

  // Plotting

  TH1F *Hmet_pt_sigQCD_new = (TH1F *)Hmet_pt_sigQCD->Clone("Hmet_pt_sigQCD_new");
  TH1F *Hmet_pt_sigEWK_new = (TH1F *)Hmet_pt_sigEWK->Clone("Hmet_pt_sigEWK_new");

  TCanvas *c1 = new TCanvas("c1", "MET", 1400., 600., 600, 600);

  TPad *pad_met = new TPad("pad_met", "This is pad_met", 0.01, 0.30, 1., 1.);
  TPad *pad_met2 = new TPad("pad_met2", "This is pad_met2", 0.01, 0.01, 1., 0.30);

  pad_met->SetBorderSize(0);
  pad_met->SetBottomMargin(0.02);
  pad_met->Draw();
  pad_met2->SetBottomMargin(0.35);
  pad_met2->SetTopMargin(0.0);
  pad_met2->SetBorderSize(0);
  pad_met2->Draw();
  pad_met->cd();

  // Draw contributions in stacked histogtams
  //
  Hmet_pt_sigQCD->SetFillColor(38);
  Hmet_pt_sigEWK->SetFillColor(6);
  Hmet_pt_Zjets->SetFillColor(8);
  Hmet_pt_WZ->SetFillColor(4);
  Hmet_pt_WW->SetFillColor(5);
  Hmet_pt_Wt->SetFillColor(7);
  Hmet_pt_tt->SetFillColor(42);
  Hmet_pt_othr->SetFillColor(30);
  Hmet_pt_trib->SetFillColor(46);

  Hmet_pt_trib->Add(Hmet_pt_othr);
  Hmet_pt_Wt->Add(Hmet_pt_trib);
  Hmet_pt_WW->Add(Hmet_pt_Wt);
  Hmet_pt_WZ->Add(Hmet_pt_WW);
  Hmet_pt_tt->Add(Hmet_pt_WZ);
  Hmet_pt_Zjets->Add(Hmet_pt_tt);
  Hmet_pt_sigEWK->Add(Hmet_pt_Zjets);
  Hmet_pt_sigQCD->Add(Hmet_pt_sigEWK);

  Hmet_pt_sigQCD->Draw("hist");
  Hmet_pt_sigEWK->Draw("histsame");
  Hmet_pt_Zjets->Draw("histsame");
  Hmet_pt_tt->Draw("histsame");
  Hmet_pt_WZ->Draw("histsame");
  Hmet_pt_WW->Draw("histsame");
  Hmet_pt_Wt->Draw("histsame");
  Hmet_pt_trib->Draw("histsame");
  Hmet_pt_othr->Draw("histsame");
  Hmet_pt_data->Draw("sameE0X0");

  Hmet_pt_sigQCD->GetYaxis()->SetTitleSize(0.06);
  Hmet_pt_sigQCD->GetYaxis()->SetTitleOffset(1.11);
  Hmet_pt_sigQCD->GetYaxis()->SetLabelSize(0.05);
  Hmet_pt_sigQCD->GetXaxis()->SetLabelSize(0.00);
  Hmet_pt_sigQCD->GetYaxis()->SetTitle("Events");

  TLatex *tex = new TLatex(0.2, 0.8, "#intL dt = 138.9 fb^{-1}, #sqrt{s} = 13 TeV, ee+#mu#mu, mc16ade"); // TeX format
  tex->SetNDC();
  tex->SetTextFont(42);
  tex->SetLineWidth(2);
  tex->Draw();

  TLegend *leg = new TLegend(1, 0.3, 0.9, 0.75, NULL, "brNDC");
  TLegendEntry *leg_entry;
  leg_entry = leg->AddEntry(Hmet_pt_data, "Data", "lp");
  leg_entry = leg->AddEntry(Hmet_pt_sigQCD, "ZZQCD", "f");
  leg_entry = leg->AddEntry(Hmet_pt_sigEWK, "ZZEWK", "f");
  leg_entry = leg->AddEntry(Hmet_pt_Zjets, "Zjets", "f");
  leg_entry = leg->AddEntry(Hmet_pt_tt, "TOP", "f");
  leg_entry = leg->AddEntry(Hmet_pt_WZ, "WZ", "f");
  leg_entry = leg->AddEntry(Hmet_pt_WW, "WW", "f");
  leg_entry = leg->AddEntry(Hmet_pt_Wt, "Wt", "f");
  leg_entry = leg->AddEntry(Hmet_pt_trib, "VVV", "f");
  leg_entry = leg->AddEntry(Hmet_pt_othr, "other", "f");
  leg->SetLineColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

  pad_met->RedrawAxis();
  pad_met2->cd();

  // to plot significance

  Double_t Z_bin = 0;
  for (int bin = 1; bin < Hmet_pt_Zjets->GetSize(); ++bin)
  {
    Double_t B = Hmet_pt_Zjets->GetBinContent(bin);
    Double_t S = Hmet_pt_sigQCD_new->GetBinContent(bin) + Hmet_pt_sigEWK_new->GetBinContent(bin);

    if (B > 0 && S > 0)
    {
      Double_t Z_bin = sqrt(2 * ((S + B) * log(1 + (S / B)) - S));
      // cout << "B=" << B << ", S=" << S << ", Z_bin=" << Z_bin << endl;
      Hmet_pt_signal_sing->SetBinContent(bin, Z_bin);
    }
  }

  TH1F *h01 = (TH1F *)Hmet_pt_signal_sing->Clone("h01");
  h01->SetLineColor(4);
  h01->SetMarkerColor(4);
  h01->Draw("E0X0");

  // do plot Data/MC

  // TH1F*h00=(TH1F*)Hmet_pt_sigQCD->Clone("h00");
  // TH1F*h01=(TH1F*)Hmet_pt_data->Clone("h01");
  // h00->Sumw2(1);
  // h01->Sumw2(1);
  // h01->Divide(h00);
  // h01->Draw("E0X0");
  // TLine *line=new TLine(90,1,1000,1);
  // line->SetLineStyle(2);
  // line->SetLineWidth(2);
  // line->Draw("same");

  h01->GetXaxis()->SetTitleSize(0.15);
  h01->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]");
  h01->GetXaxis()->SetTitleOffset(1.1);
  h01->GetXaxis()->SetLabelSize(0.13);

  h01->GetYaxis()->SetTitleOffset(0.45);
  h01->GetYaxis()->SetTitleSize(0.135);
  h01->GetYaxis()->SetLabelSize(0.11);
  // h01->GetYaxis()->SetTitle("#frac{Data}{MC}");
  h01->GetYaxis()->SetTitle("Sig. Significance");

  TH1F *Hmet_pt_signif_sigQCD_new = (TH1F *)Hmet_pt_signif_sigQCD->Clone("Hmet_pt_signif_sigQCD_new");
  TH1F *Hmet_pt_signif_sigEWK_new = (TH1F *)Hmet_pt_signif_sigEWK->Clone("Hmet_pt_signif_sigEWK_new");

  TCanvas *c2 = new TCanvas("c2", "MET significance", 1400., 600., 600, 600);

  TPad *pad_met_signif = new TPad("pad_met_signif", "This is pad_met_signif", 0.01, 0.30, 1., 1.);
  TPad *pad_met_signif2 = new TPad("pad_met_signif2", "This is pad_met_signif2", 0.01, 0.01, 1., 0.30);

  pad_met_signif->SetBorderSize(0);
  pad_met_signif->SetBottomMargin(0.02);
  pad_met_signif->Draw();
  pad_met_signif2->SetBottomMargin(0.35);
  pad_met_signif2->SetTopMargin(0.0);
  pad_met_signif2->SetBorderSize(0);
  pad_met_signif2->Draw();
  pad_met_signif->cd();

  // Draw contributions in stacked histogtams
  //
  Hmet_pt_signif_sigQCD->SetFillColor(38);
  Hmet_pt_signif_sigEWK->SetFillColor(6);
  Hmet_pt_signif_Zjets->SetFillColor(8);
  Hmet_pt_signif_WZ->SetFillColor(4);
  Hmet_pt_signif_WW->SetFillColor(5);
  Hmet_pt_signif_Wt->SetFillColor(7);
  Hmet_pt_signif_tt->SetFillColor(42);
  Hmet_pt_signif_othr->SetFillColor(30);
  Hmet_pt_signif_trib->SetFillColor(46);

  Hmet_pt_signif_trib->Add(Hmet_pt_signif_othr);
  Hmet_pt_signif_Wt->Add(Hmet_pt_signif_trib);
  Hmet_pt_signif_WW->Add(Hmet_pt_signif_Wt);
  Hmet_pt_signif_WZ->Add(Hmet_pt_signif_WW);
  Hmet_pt_signif_tt->Add(Hmet_pt_signif_WZ);
  Hmet_pt_signif_Zjets->Add(Hmet_pt_signif_tt);
  Hmet_pt_signif_sigEWK->Add(Hmet_pt_signif_Zjets);
  Hmet_pt_signif_sigQCD->Add(Hmet_pt_signif_sigEWK);

  Hmet_pt_signif_sigQCD->Draw("hist");
  Hmet_pt_signif_sigEWK->Draw("histsame");
  Hmet_pt_signif_Zjets->Draw("histsame");
  Hmet_pt_signif_tt->Draw("histsame");
  Hmet_pt_signif_WZ->Draw("histsame");
  Hmet_pt_signif_WW->Draw("histsame");
  Hmet_pt_signif_Wt->Draw("histsame");
  Hmet_pt_signif_trib->Draw("histsame");
  Hmet_pt_signif_othr->Draw("histsame");
  Hmet_pt_signif_data->Draw("sameE0X0");

  Hmet_pt_signif_sigQCD->GetYaxis()->SetTitleSize(0.06);
  Hmet_pt_signif_sigQCD->GetYaxis()->SetTitleOffset(1.11);
  Hmet_pt_signif_sigQCD->GetYaxis()->SetLabelSize(0.05);
  Hmet_pt_signif_sigQCD->GetXaxis()->SetLabelSize(0.00);
  Hmet_pt_signif_sigQCD->GetYaxis()->SetTitle("Events");

  tex = new TLatex(0.2, 0.8, "#intL dt = 138.9 fb^{-1}, #sqrt{s} = 13 TeV, ee+#mu#mu, mc16ade");
  tex->SetNDC();
  tex->SetTextFont(42);
  tex->SetLineWidth(2);
  tex->Draw();

  leg = new TLegend(0.6, 0.3, 0.9, 0.75, NULL, "brNDC");
  // leg_entry;
  leg_entry = leg->AddEntry(Hmet_pt_signif_data, "Data", "lp");
  leg_entry = leg->AddEntry(Hmet_pt_signif_sigQCD, "ZZQCD", "f");
  leg_entry = leg->AddEntry(Hmet_pt_signif_sigEWK, "ZZEWK", "f");
  leg_entry = leg->AddEntry(Hmet_pt_signif_Zjets, "Zjets", "f");
  leg_entry = leg->AddEntry(Hmet_pt_signif_tt, "TOP", "f");
  leg_entry = leg->AddEntry(Hmet_pt_signif_WZ, "WZ", "f");
  leg_entry = leg->AddEntry(Hmet_pt_signif_WW, "WW", "f");
  leg_entry = leg->AddEntry(Hmet_pt_signif_Wt, "Wt", "f");
  leg_entry = leg->AddEntry(Hmet_pt_signif_trib, "VVV", "f");
  leg_entry = leg->AddEntry(Hmet_pt_signif_othr, "other", "f");
  leg->SetLineColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

  pad_met_signif->RedrawAxis();
  pad_met_signif2->cd();

  // to plot significance
  /*
    Double_t Z_bin=0;
    for (int bin=1; bin<Hmet_pt_signif_Zjets->GetSize(); ++bin)
    {
    Double_t B = Hmet_pt_signif_Zjets->GetBinContent(bin);
    Double_t S = Hmet_pt_signif_sigQCD_new->GetBinContent(bin)+Hmet_pt_signif_sigEWK_new->GetBinContent(bin);

    if (B>0 && S>0){
    Double_t Z_bin = sqrt(2*((S+B)*log(1+(S/B))-S));
    //cout << "B=" << B << ", S=" << S << ", Z_bin=" << Z_bin << endl;
    Hmet_pt_signif_signal_sing->SetBinContent(bin,Z_bin);
    }
    }

    h01=(TH1F*)Hmet_pt_signif_signal_sing->Clone("h01");
    h01->SetLineColor(4);
    h01->SetMarkerColor(4);
    h01->Draw("E0X0");
  */
  // do plot Data/MC

  TH1F *h00 = (TH1F *)Hmet_pt_signif_sigQCD->Clone("h00");
  h01 = (TH1F *)Hmet_pt_signif_data->Clone("h01");
  h00->Sumw2(1);
  h01->Sumw2(1);
  h01->Divide(h00);
  h01->Draw("E0X0");

  TLine *line = new TLine(0, 1, 35, 1);
  line->SetLineStyle(2);
  line->SetLineWidth(2);
  line->Draw("same");

  h01->GetXaxis()->SetTitleSize(0.15);
  h01->GetXaxis()->SetTitle("E_{T}^{miss} significance");
  h01->GetXaxis()->SetTitleOffset(1.1);
  h01->GetXaxis()->SetLabelSize(0.13);

  h01->GetYaxis()->SetTitleOffset(0.45);
  h01->GetYaxis()->SetTitleSize(0.135);
  h01->GetYaxis()->SetLabelSize(0.11);
  h01->GetYaxis()->SetTitle("#frac{Data}{MC}");
  // h01->GetYaxis()->SetTitle("Sig. Significance");

  c1->SaveAs("./MET.png");
  c2->SaveAs("./MET_S.png");

  return 0;
}
