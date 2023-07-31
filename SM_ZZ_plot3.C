
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
std::vector<float> DoReco(TTree *tree, TH1F *hist, TH1F *hist_signif, int MC)
{

  TH1::SetDefaultSumw2(kTRUE);

  Int_t nentries = (Int_t)tree->GetEntries();

  Float_t M2Lep = 0.; // Defines the flags needed
  Float_t M2Lep_signif;
  Float_t met_tst = 0.;
  Float_t met_signif = 0.;
  Float_t dMetZPhi = 0.;
  Float_t frac_pT = 0.;
  Float_t MetOHT = 0.;
  Float_t dPhiJ100met = 0.;
  Float_t dLepR = 0.;
  Int_t n_bjets = 0.;
  Float_t n_bjets_signif = 0.;
  Int_t n_jets = 0.;
  Float_t leading_pT_lepton = 0;
  Float_t subleading_pT_lepton = 0;
  Float_t detajj = 0;
  Float_t mjj = 0;
  Float_t leading_jet_pt = 0;
  Float_t second_jet_pt = 0;
  Int_t event_3CR = 0.;
  Int_t event_type = 0.;
  Float_t weight = 1.;

  std::vector<float> events;
  events.clear();

  // For all available variables see https://twiki.cern.ch/twiki/bin/view/AtlasProtected/HZZllvvMinitreeFullRun2
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
  //
  tree->SetBranchAddress("event_3CR", &event_3CR);
  tree->SetBranchAddress("event_type", &event_type);

  double signal = 0.;
  double signaler = 0.;
  double Norm = 1;

  // give each event the appropriate weight (data events should have weight = 1)
  // global_weight=weight_pileup*weight_gen*weight_exp*weight_trig*weight_jets*weight_jvt + normalization to the MC cross section and L=139 fb^-1
  if (MC == 1)
  {
    tree->SetBranchAddress("global_weight", &weight);
  }

  // Loop over events
  for (Int_t i = 0; i < nentries; i++)
  {
    tree->GetEntry(i);

    // ----- Define signal region -----  Change as needed:
    //
    // inclusive selection:
    //---------------------
    // RECO cut-based:

    if (event_3CR == 0 && (event_type == 0 || event_type == 1) &&
        leading_pT_lepton > 30 && subleading_pT_lepton > 20 && M2Lep > 76 && M2Lep < 116 && n_bjets < 1 &&
        dLepR < 1.9 && dMetZPhi > 2.6 && met_tst > 90 && MetOHT > 0.6)

    // if (event_3CR == 0 && (event_type == 0 || event_type == 1) &&
    //       leading_pT_lepton > 30 && subleading_pT_lepton > 20  && M2Lep > 80 && M2Lep < 100 && n_bjets < 1 &&
    //       dLepR < 1.8 && dMetZPhi > 2.7 && met_tst > 110 && MetOHT > 0.65)

    {
      signal = signal + weight;              // signal yield is sum of weights
      signaler = signaler + weight * weight; // keep Sum(weights^2) for calculating the error on the signal yield
      hist->Fill(met_tst, weight);        // Fills histogram but with wieght
      hist_signif->Fill(met_tst, weight);
    }
  }

  cout << "N=" << signal << "+-" << sqrt(signaler) << endl; // signal yield and error on this yield

  events.push_back(signal);
  events.push_back(sqrt(signaler));
  return events;
}

// the main method
//
int SM_ZZ_plot3()
{

  gStyle->SetLegendFont(42);
  gStyle->SetPalette(1);

  // TString sample = "/Results/trial2.root";
  // ---convert TString to string
  // std::string sample_string(sample.Data());
  // ---convert string to TString
  // TString sample_TString(sample_string)

  // Preliminary minitrees for ZZllvvjj analysis
  TFile *file_data = new TFile("../data/data_fullRun2_MET90.root");
  TTree *tree_data = file_data->Get<TTree>("tree_PFLOW");

  TFile *file_signalQCD = new TFile("../data/mc16ade_QCD_EWKcor_ZZllvv_MET90.root"); // QCD signal with EWK corrections
  TTree *tree_signalQCD = file_signalQCD->Get<TTree>("tree_PFLOW");

  TFile *file_signalEWK = new TFile("../data/mc16ade_EWK_ZZllvvjj_MET90.root"); // EW signal
  TTree *tree_signalEWK = file_signalEWK->Get<TTree>("tree_PFLOW");

  TFile *file_Zjets = new TFile("../data/mc16ade_Zjets_MET90.root");
  TTree *tree_Zjets = file_Zjets->Get<TTree>("tree_PFLOW");

  TFile *file_WZ = new TFile("../data/mc16ade_WZ_MET90.root");
  TTree *tree_WZ = file_WZ->Get<TTree>("tree_PFLOW");

  TFile *file_WW = new TFile("../data/mc16ade_WW_Pow_MET90.root");
  TTree *tree_WW = file_WW->Get<TTree>("tree_PFLOW");

  TFile *file_Wt = new TFile("../data/mc16ade_Wt_MET90.root");
  TTree *tree_Wt = file_Wt->Get<TTree>("tree_PFLOW");

  TFile *file_tt = new TFile("../data/mc16ade_TOP_MET90.root");
  TTree *tree_tt = file_tt->Get<TTree>("tree_PFLOW");

  TFile *file_trib = new TFile("../data/mc16ade_VVV_MET90.root");
  TTree *tree_trib = file_trib->Get<TTree>("tree_PFLOW");

  TFile *file_othr = new TFile("../data/mc16ade_ttV_ttVV_MET90.root");
  TTree *tree_othr = file_othr->Get<TTree>("tree_PFLOW");

  TH1::SetDefaultSumw2(kTRUE);

  Float_t xbins[21] = {90, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 345, 375, 410, 450, 500, 580, 700, 1000}; // Unequal bins
  Float_t binmin = 80;
  Float_t binmax = 500;

  TH1F *hist_sigQCD = new TH1F("N_bjets", "N_{bjets}", 20, binmin, binmax);
  TH1F *hist_sigEWK = new TH1F("hist_sigEWK", "hist_sigEWK", 20, binmin, binmax);
  TH1F *hist_data = new TH1F("hist_data ", "hist_data ", 20, binmin, binmax);
  TH1F *hist_Zjets = new TH1F("hist_Zjets ", "hist_Zjets ", 20, binmin, binmax);
  TH1F *hist_WZ = new TH1F("hist_WZ ", "hist_WZ ", 20, binmin, binmax);
  TH1F *hist_WW = new TH1F("hist_WW ", "hist_WW ", 20, binmin, binmax);
  TH1F *hist_Wt = new TH1F("hist_Wt ", "hist_Wt ", 20, binmin, binmax);
  TH1F *hist_tt = new TH1F("hist_tt ", "hist_tt ", 20, binmin, binmax);
  TH1F *hist_trib = new TH1F("hist_trib ", "hist_trib ", 20, binmin, binmax);
  TH1F *hist_othr = new TH1F("hist_othr ", "hist_othr ", 20, binmin, binmax);
  TH1F *hist_signal_sing = new TH1F("hist_signal_sing ", "hist_signal_sing", 20, binmin, binmax);

  TH1F *hist_signif_sigQCD = new TH1F("hist_signif_sigQCD", "hist_signif_sigQCD", 35, 0, 35);
  TH1F *hist_signif_sigEWK = new TH1F("hist_signif_sigEWK", "hist_signif_sigEWK", 35, 0, 35);
  TH1F *hist_signif_data = new TH1F("hist_signif_data ", "hist_signif_data ", 35, 0, 35);
  TH1F *hist_signif_Zjets = new TH1F("hist_signif_Zjets ", "hist_signif_Zjets ", 35, 0, 35);
  TH1F *hist_signif_WZ = new TH1F("hist_signif_WZ ", "hist_signif_WZ ", 35, 0, 35);
  TH1F *hist_signif_WW = new TH1F("hist_signif_WW ", "hist_signif_WW ", 35, 0, 35);
  TH1F *hist_signif_Wt = new TH1F("hist_signif_Wt ", "hist_signif_Wt ", 35, 0, 35);
  TH1F *hist_signif_tt = new TH1F("hist_signif_tt ", "hist_signif_tt ", 35, 0, 35);
  TH1F *hist_signif_trib = new TH1F("hist_signif_trib ", "hist_signif_trib", 35, 0, 35);
  TH1F *hist_signif_othr = new TH1F("hist_signif_othr ", "hist_signif_othr", 35, 0, 35);

  TH1F *hist_signif_signal_sing = new TH1F("hist_signif_signal_sing ", "hist_signif_signal_sing", 35, 0, 35);
  TH1F *hist_signif_signal_pur = new TH1F("hist_signif_signal_pur ", "hist_signif_signal_pur", 35, 0, 35);

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
  events_data = DoReco(tree_data, hist_data, hist_signif_data, 0); // DoReco Fills the vectors
  cout << "===Signal QCD ZZ===" << endl;
  events_NsigQCD = DoReco(tree_signalQCD, hist_sigQCD, hist_signif_sigQCD, 1);
  cout << "===Signal EWK ZZ===" << endl;
  events_NsigEWK = DoReco(tree_signalEWK, hist_sigEWK, hist_signif_sigEWK, 1);
  cout << "===Background===" << endl;
  cout << "Zjets" << endl;
  events_NZjets = DoReco(tree_Zjets, hist_Zjets, hist_signif_Zjets, 1);
  cout << "WZ" << endl;
  events_NWZ = DoReco(tree_WZ, hist_WZ, hist_signif_WZ, 1);
  cout << "tt" << endl;
  events_Ntt = DoReco(tree_tt, hist_tt, hist_signif_tt, 1);
  cout << "WW" << endl;
  events_NWW = DoReco(tree_WW, hist_WW, hist_signif_WW, 1);
  cout << "Wt" << endl;
  events_NWt = DoReco(tree_Wt, hist_Wt, hist_signif_Wt, 1);
  cout << "trib" << endl;
  events_Ntrib = DoReco(tree_trib, hist_trib, hist_signif_trib, 1);
  cout << "othr" << endl;
  events_Nothr = DoReco(tree_othr, hist_othr, hist_signif_othr, 1);

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

  TH1F *hist_sigQCD_new = (TH1F *)hist_sigQCD->Clone("hist_sigQCD_new");
  TH1F *hist_sigEWK_new = (TH1F *)hist_sigEWK->Clone("hist_sigEWK_new");

  TCanvas *c1 = new TCanvas("c1", "MET", 1400., 600., 600, 600);

  TPad *pad_met = new TPad("pad_met", "This is pad_met", 0.01, 0.30, 1., 1.);
  TPad *pad_met2 = new TPad("pad_met2", "This is pad_met2", 0.01, 0.01, 1., 0.30);

  pad_met->SetBorderSize(0);
  pad_met->SetBottomMargin(0.02);
  pad_met->Draw();
  pad_met2->SetBottomMargin(0.35);
  pad_met2->SetGrid();
  pad_met2->SetTopMargin(0.0);
  pad_met2->SetBorderSize(0);
  pad_met2->Draw();
  pad_met->cd();

  // Draw contributions in stacked histogtams
  //
  hist_sigQCD->SetFillColor(38);
  hist_sigEWK->SetFillColor(6);
  hist_Zjets->SetFillColor(8);
  hist_WZ->SetFillColor(4);
  hist_WW->SetFillColor(5);
  hist_Wt->SetFillColor(7);
  hist_tt->SetFillColor(42);
  hist_othr->SetFillColor(30);
  hist_trib->SetFillColor(46);

  hist_trib->Add(hist_othr);
  hist_Wt->Add(hist_trib);
  hist_WW->Add(hist_Wt);
  hist_WZ->Add(hist_WW);
  hist_tt->Add(hist_WZ);
  hist_Zjets->Add(hist_tt);
  hist_sigEWK->Add(hist_Zjets);
  hist_sigQCD->Add(hist_sigEWK);

  hist_sigQCD->Draw("hist");
  hist_sigEWK->Draw("histsame");
  hist_Zjets->Draw("histsame");
  hist_tt->Draw("histsame");
  hist_WZ->Draw("histsame");
  hist_WW->Draw("histsame");
  hist_Wt->Draw("histsame");
  hist_trib->Draw("histsame");
  hist_othr->Draw("histsame");
  hist_data->Draw("sameE0X0");

  hist_sigQCD->GetYaxis()->SetTitleSize(0.06);
  hist_sigQCD->GetYaxis()->SetTitleOffset(0.85);
  hist_sigQCD->GetYaxis()->SetLabelSize(0.04);
  hist_sigQCD->GetXaxis()->SetLabelSize(0.00);
  hist_sigQCD->GetYaxis()->SetTitle("Events");

  TLatex *tex1 = new TLatex(0.15, 0.8, "#intL dt = 138.9 fb^{-1}");
  tex1->SetNDC();
  tex1->SetTextFont(42);
  tex1->SetLineWidth(2);
  tex1->Draw();

  TLatex *tex2 = new TLatex(0.15, 0.7, "#sqrt{s} = 13 TeV");
  tex2->SetNDC();
  tex2->SetTextFont(42);
  tex2->SetLineWidth(2);
  tex2->Draw();

  TLatex *tex3 = new TLatex(0.15, 0.6, "ee+#mu#mu, mc16ade");
  tex3->SetNDC();
  tex3->SetTextFont(42);
  tex3->SetLineWidth(2);
  tex3->Draw();

  TLegend *leg = new TLegend(0.6, 0.3, 0.7, 0.75, NULL, "brNDC");
  TLegendEntry *leg_entry;
  leg_entry = leg->AddEntry(hist_data, "Data", "lp");
  leg_entry = leg->AddEntry(hist_sigQCD, "ZZQCD", "f");
  leg_entry = leg->AddEntry(hist_sigEWK, "ZZEWK", "f");
  leg_entry = leg->AddEntry(hist_Zjets, "Zjets", "f");
  leg_entry = leg->AddEntry(hist_tt, "TOP", "f");
  leg_entry = leg->AddEntry(hist_WZ, "WZ", "f");
  leg_entry = leg->AddEntry(hist_WW, "WW", "f");
  leg_entry = leg->AddEntry(hist_Wt, "Wt", "f");
  leg_entry = leg->AddEntry(hist_trib, "VVV", "f");
  leg_entry = leg->AddEntry(hist_othr, "other", "f");
  leg->SetLineColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

  pad_met->RedrawAxis();
  pad_met2->cd();

  // to plot significance

  Double_t Z_bin = 0;
  for (int bin = 1; bin < hist_Zjets->GetSize(); ++bin)
  {
    Double_t B = hist_Zjets->GetBinContent(bin);
    Double_t S = hist_sigQCD_new->GetBinContent(bin) + hist_sigEWK_new->GetBinContent(bin);

    if (B > 0 && S > 0)
    {
      Double_t Z_bin = sqrt(2 * ((S + B) * log(1 + (S / B)) - S));
      // cout << "B=" << B << ", S=" << S << ", Z_bin=" << Z_bin << endl;
      hist_signal_sing->SetBinContent(bin, Z_bin);
    }
  }

  TH1F *h01 = (TH1F *)hist_signal_sing->Clone("h01");
  h01->SetLineColor(4);
  h01->SetMarkerColor(4);
  h01->Draw("E0X0");

  // do plot Data/MC
  /*
    TH1F*h00=(TH1F*)hist_sigQCD->Clone("h00");
    TH1F*h01=(TH1F*)hist_data->Clone("h01");
    h00->Sumw2(1);
    h01->Sumw2(1);
    h01->Divide(h00);
    h01->Draw("E0X0");

    TLine *line=new TLine(90,1,1000,1);
    line->SetLineStyle(2);
    line->SetLineWidth(2);
    line->Draw("same");
  */

  h01->GetXaxis()->SetTitleSize(0.15);
  h01->GetXaxis()->SetTitle("N_{bjets}");
  h01->GetXaxis()->SetTitleOffset(1.1);
  h01->GetXaxis()->SetLabelSize(0.13);

  h01->GetYaxis()->SetTitleOffset(0.35);
  h01->GetYaxis()->SetTitleSize(0.110);
  h01->GetYaxis()->SetLabelSize(0.08);
  // h01->GetYaxis()->SetTitle("#frac{Data}{MC}");
  h01->GetYaxis()->SetTitle("Sig. Significance");

  TH1F *hist_signif_sigQCD_new = (TH1F *)hist_signif_sigQCD->Clone("hist_signif_sigQCD_new");
  TH1F *hist_signif_sigEWK_new = (TH1F *)hist_signif_sigEWK->Clone("hist_signif_sigEWK_new");

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
  hist_signif_sigQCD->SetFillColor(38);
  hist_signif_sigEWK->SetFillColor(6);
  hist_signif_Zjets->SetFillColor(8);
  hist_signif_WZ->SetFillColor(4);
  hist_signif_WW->SetFillColor(5);
  hist_signif_Wt->SetFillColor(7);
  hist_signif_tt->SetFillColor(42);
  hist_signif_othr->SetFillColor(30);
  hist_signif_trib->SetFillColor(46);

  hist_signif_trib->Add(hist_signif_othr);
  hist_signif_Wt->Add(hist_signif_trib);
  hist_signif_WW->Add(hist_signif_Wt);
  hist_signif_WZ->Add(hist_signif_WW);
  hist_signif_tt->Add(hist_signif_WZ);
  hist_signif_Zjets->Add(hist_signif_tt);
  hist_signif_sigEWK->Add(hist_signif_Zjets);
  hist_signif_sigQCD->Add(hist_signif_sigEWK);

  hist_signif_sigQCD->Draw("hist");
  hist_signif_sigEWK->Draw("histsame");
  hist_signif_Zjets->Draw("histsame");
  hist_signif_tt->Draw("histsame");
  hist_signif_WZ->Draw("histsame");
  hist_signif_WW->Draw("histsame");
  hist_signif_Wt->Draw("histsame");
  hist_signif_trib->Draw("histsame");
  hist_signif_othr->Draw("histsame");
  hist_signif_data->Draw("sameE0X0");

  hist_signif_sigQCD->GetYaxis()->SetTitleSize(0.04);
  hist_signif_sigQCD->GetYaxis()->SetTitleOffset(1.);
  hist_signif_sigQCD->GetYaxis()->SetLabelSize(0.05);
  hist_signif_sigQCD->GetXaxis()->SetLabelSize(0.00);
  hist_signif_sigQCD->GetYaxis()->SetTitle("Events");

  tex1 = new TLatex(0.15, 0.8, "#intL dt = 138.9 fb^{-1}");
  tex1->SetNDC();
  tex1->SetTextFont(42);
  tex1->SetLineWidth(2);
  tex1->Draw();

  tex2 = new TLatex(0.15, 0.7, "#sqrt{s} = 13 TeV, ee+#mu#mu, mc16ade");
  tex2->SetNDC();
  tex2->SetTextFont(42);
  tex2->SetLineWidth(2);
  tex2->Draw();

  leg = new TLegend(0.6, 0.3, 0.7, 0.75, NULL, "brNDC");
  // leg_entry;
  leg_entry = leg->AddEntry(hist_signif_data, "Data", "lp");
  leg_entry = leg->AddEntry(hist_signif_sigQCD, "ZZQCD", "f");
  leg_entry = leg->AddEntry(hist_signif_sigEWK, "ZZEWK", "f");
  leg_entry = leg->AddEntry(hist_signif_Zjets, "Zjets", "f");
  leg_entry = leg->AddEntry(hist_signif_tt, "TOP", "f");
  leg_entry = leg->AddEntry(hist_signif_WZ, "WZ", "f");
  leg_entry = leg->AddEntry(hist_signif_WW, "WW", "f");
  leg_entry = leg->AddEntry(hist_signif_Wt, "Wt", "f");
  leg_entry = leg->AddEntry(hist_signif_trib, "VVV", "f");
  leg_entry = leg->AddEntry(hist_signif_othr, "other", "f");
  leg->SetLineColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

  pad_met_signif->RedrawAxis();
  pad_met_signif2->cd();

  // to plot significance
  /*
    Double_t Z_bin=0;
    for (int bin=1; bin<hist_signif_Zjets->GetSize(); ++bin)
    {
    Double_t B = hist_signif_Zjets->GetBinContent(bin);
    Double_t S = hist_signif_sigQCD_new->GetBinContent(bin)+hist_signif_sigEWK_new->GetBinContent(bin);

    if (B>0 && S>0){
    Double_t Z_bin = sqrt(2*((S+B)*log(1+(S/B))-S));
    //cout << "B=" << B << ", S=" << S << ", Z_bin=" << Z_bin << endl;
    hist_signif_signal_sing->SetBinContent(bin,Z_bin);
    }
    }

    h01=(TH1F*)hist_signif_signal_sing->Clone("h01");
    h01->SetLineColor(4);
    h01->SetMarkerColor(4);
    h01->Draw("E0X0");
  */
  // do plot Data/MC

  TH1F *h00 = (TH1F *)hist_signif_sigQCD->Clone("h00");
  h01 = (TH1F *)hist_signif_data->Clone("h01");
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

  // c1->SaveAs("../figs/MET_TST.pdf");
  // c2->SaveAs("./MET_S.png");

  return 0;
}
