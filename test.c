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
using namespace std;

// gStyle->SetLegendFont(42);
// gStyle->SetPalette(1);

// Method to fill hitograms (MET and MET_significance) with events properly weighted
//

std::vector<float> DoReco(TFile *file, TH1F *Hmet_pt, TH1F *Hmet_pt_signif, int MC, double testDLepR, double testdMetZPhi, double testmet_tst, double testMetOHT)
{

  TH1::SetDefaultSumw2(kTRUE);

  TTree *tree = file->Get<TTree>("tree_PFLOW");

  Int_t nentries = (Int_t)tree->GetEntries();
  // cout << "number of entries" <<nentries <<endl;

  Float_t M2Lep = 0.;
  Float_t met_tst = 0.;
  Float_t met_signif = 0.;
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

    // FIDUCIAL
    if (event_3CR == 0 && (event_type == 0 || event_type == 1) && 
    leading_pT_lepton > 30 && subleading_pT_lepton > 20 && M2Lep > 76 && M2Lep < 116 && n_bjets < 1 && 
    dLepR < testDLepR && dMetZPhi > testdMetZPhi && met_tst > testmet_tst && MetOHT > testMetOHT && met_signif > 0) // met_tst/sqrt

    // FID-BDT:

    // FID-DNN:

    // VBS selection:

    // Nominal:

    // Nominal + METSignif :

    // Inclusive + di-jet variables:

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
  // Hmet_pt->Reset();
  // Hmet_pt_signif->Recreate;
  return events;
}

// the main method
//

int test()
{

  Float_t Z;
  Float_t maxsig = -1;
  Int_t a = -1;
  Int_t b = -1;
  Int_t c = -1;
  Int_t z = -1;
  Int_t metritiri = 0;

  double testDLepR[4] = {1.8, 2.0, 2.2, 2.5};
  double testdMetZPhi[4] = {2.2, 2.5, 2.7, 3.0};
  double testmet_tst[1] = {90};
  double testMetOHT[4] = {0.5, 0.6, 0.65, 0.7};

  // Preliminary minitrees for ZZllvvjj analysis
  TFile *file_data = new TFile("./data_fullRun2_MET90.root");

  TFile *file_signalQCD = new TFile("./mc16ade_QCD_EWKcor_ZZllvv_MET90.root"); // QCD signal with EWK corrections

  TFile *file_signalEWK = new TFile("./mc16ade_EWK_ZZllvvjj_MET90.root"); // EW signal

  TFile *file_Zjets = new TFile("./mc16ade_Zjets_MET90.root");

  TFile *file_WZ = new TFile("./mc16ade_WZ_MET90.root");

  TFile *file_WW = new TFile("./mc16ade_WW_Pow_MET90.root");

  TFile *file_Wt = new TFile("./mc16ade_Wt_MET90.root");

  TFile *file_tt = new TFile("./mc16ade_TOP_MET90.root");

  TFile *file_trib = new TFile("./mc16ade_VVV_MET90.root");

  TFile *file_othr = new TFile("./mc16ade_ttV_ttVV_MET90.root");

  for (Int_t i = 0; i < 4; i++)
  {
    for (Int_t j = 0; j < 4; j++)
    {
      for (Int_t k = 0; k < 1; k++)
      {
        for (Int_t l = 0; l < 4; l++)
        {

          printf("------------ Combination: i = %d, j = %d, k= %d, l = %d-------------------------\n", i, j, k, l);

          gStyle->SetLegendFont(42);
          gStyle->SetPalette(1);

          TH1::SetDefaultSumw2(kTRUE);

          Float_t xbins[21] = {90, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 345, 375, 410, 450, 500, 580, 700, 1000};

          TH1F *Hmet_pt_sigQCD = new TH1F("Hmet_pt_sigQCD", "Hmet_pt_sigQCD", 20, xbins);
          TH1F *Hmet_pt_sigEWK = new TH1F("Hmet_pt_sigEWK", "Hmet_pt_sigEWK", 20, xbins);
          TH1F *Hmet_pt_data = new TH1F("Hmet_pt_data", "Hmet_pt_data", 20, xbins);
          TH1F *Hmet_pt_Zjets = new TH1F("Hmet_pt_Zjets", "Hmet_pt_Zjets", 20, xbins);
          TH1F *Hmet_pt_WZ = new TH1F("Hmet_pt_WZ", "Hmet_pt_WZ", 20, xbins);
          TH1F *Hmet_pt_WW = new TH1F("Hmet_pt_WW", "Hmet_pt_WW", 20, xbins);
          TH1F *Hmet_pt_Wt = new TH1F("Hmet_pt_Wt", "Hmet_pt_Wt", 20, xbins);
          TH1F *Hmet_pt_tt = new TH1F("Hmet_pt_tt", "Hmet_pt_tt", 20, xbins);
          TH1F *Hmet_pt_trib = new TH1F("Hmet_pt_trib", "Hmet_pt_trib", 20, xbins);
          TH1F *Hmet_pt_othr = new TH1F("Hmet_pt_othr", "Hmet_pt_othr", 20, xbins);
          TH1F *Hmet_pt_signal_sing = new TH1F("Hmet_pt_signal_sing", "Hmet_pt_signal_sing", 20, xbins);

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
          events_data = DoReco(file_data, Hmet_pt_data, Hmet_pt_signif_data, 0, testDLepR[i], testdMetZPhi[j], testmet_tst[k], testMetOHT[l]);
          cout << "===Signal QCD ZZ===" << endl;
          events_NsigQCD = DoReco(file_signalQCD, Hmet_pt_sigQCD, Hmet_pt_signif_sigQCD, 1, testDLepR[i], testdMetZPhi[j], testmet_tst[k], testMetOHT[l]);
          cout << "===Signal EWK ZZ===" << endl;
          events_NsigEWK = DoReco(file_signalEWK, Hmet_pt_sigEWK, Hmet_pt_signif_sigEWK, 1, testDLepR[i], testdMetZPhi[j], testmet_tst[k], testMetOHT[l]);
          cout << "===Background===" << endl;
          cout << "Zjets" << endl;
          events_NZjets = DoReco(file_Zjets, Hmet_pt_Zjets, Hmet_pt_signif_Zjets, 1, testDLepR[i], testdMetZPhi[j], testmet_tst[k], testMetOHT[l]);
          cout << "WZ" << endl;
          events_NWZ = DoReco(file_WZ, Hmet_pt_WZ, Hmet_pt_signif_WZ, 1, testDLepR[i], testdMetZPhi[j], testmet_tst[k], testMetOHT[l]);
          cout << "tt" << endl;
          events_Ntt = DoReco(file_tt, Hmet_pt_tt, Hmet_pt_signif_tt, 1, testDLepR[i], testdMetZPhi[j], testmet_tst[k], testMetOHT[l]);
          cout << "WW" << endl;
          events_NWW = DoReco(file_WW, Hmet_pt_WW, Hmet_pt_signif_WW, 1, testDLepR[i], testdMetZPhi[j], testmet_tst[k], testMetOHT[l]);
          cout << "Wt" << endl;
          events_NWt = DoReco(file_Wt, Hmet_pt_Wt, Hmet_pt_signif_Wt, 1, testDLepR[i], testdMetZPhi[j], testmet_tst[k], testMetOHT[l]);
          cout << "trib" << endl;
          events_Ntrib = DoReco(file_trib, Hmet_pt_trib, Hmet_pt_signif_trib, 1, testDLepR[i], testdMetZPhi[j], testmet_tst[k], testMetOHT[l]);
          cout << "othr" << endl;
          events_Nothr = DoReco(file_othr, Hmet_pt_othr, Hmet_pt_signif_othr, 1, testDLepR[i], testdMetZPhi[j], testmet_tst[k], testMetOHT[l]);

          //  file_data->Close();
          //  file_signalQCD->Close();
          //  file_signalEWK->Close();
          //  file_Zjets->Close();
          //  file_WZ->Close();
          //  file_WW->Close();
          //  file_tt->Close();
          //  file_Wt->Close();
          //  file_trib->Close();
          //  file_othr->Close();

          // Some useful outputs
          double totalBKG = events_NZjets.at(0) + events_Ntt.at(0) + events_NWt.at(0) + events_NWW.at(0) + events_NWZ.at(0) + events_Ntrib.at(0) + events_Nothr.at(0);
          double totalBKG_er = sqrt(events_NZjets.at(1) * events_NZjets.at(1) + events_Ntt.at(1) * events_Ntt.at(1) + events_NWt.at(1) * events_NWt.at(1) + events_NWW.at(1) * events_NWW.at(1) +
                                    events_NWZ.at(1) * events_NWZ.at(1) + events_Ntrib.at(1) * events_Ntrib.at(1) + events_Nothr.at(1) * events_Nothr.at(1));
          cout << "total Bkg=" << totalBKG << "+-" << totalBKG_er << endl;
          cout << "S/B=" << (events_NsigQCD.at(0) + events_NsigEWK.at(0)) / totalBKG << endl;
          Double_t S = (events_NsigQCD.at(0) + events_NsigEWK.at(0));
          Double_t B = totalBKG;
          Z = sqrt(2 * ((S + B) * log(1 + (S / B)) - S));
          cout << " Significance = " << Z << endl;

          // return 0; // scaning in overall of histogram. Take all the bins at once.

          // Plotting

          TH1F *Hmet_pt_sigQCD_new = (TH1F *)Hmet_pt_sigQCD->Clone("Hmet_pt_sigQCD_new");
          TH1F *Hmet_pt_sigEWK_new = (TH1F *)Hmet_pt_sigEWK->Clone("Hmet_pt_sigEWK_new");

          gROOT->SetBatch(kTRUE); // gia na min emfanizei tipota apo ta istogrammata
          TCanvas *c1 = new TCanvas("c1", "MET", 0., 0., 600, 600);

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
          TLatex *tex = new TLatex(0.2, 0.8, "#intL dt = 138.9 fb^{-1}, #sqrt{s} = 13 TeV, ee+#mu#mu, mc16ade");
          tex->SetNDC();
          tex->SetTextFont(42);
          tex->SetLineWidth(2);
          tex->Draw();

          TLegend *leg = new TLegend(0.6, 0.3, 0.9, 0.75, NULL, "brNDC");
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

          // Double_t Z_bin=0;
          for (int bin = 1; bin < Hmet_pt_Zjets->GetSize(); ++bin)
          {
            // bin per bin significance
            // Double_t B = Hmet_pt_Zjets->GetBinContent(bin);
            // Double_t S = Hmet_pt_sigQCD_new->GetBinContent(bin)+Hmet_pt_sigEWK_new->GetBinContent(bin);

            // count from this bin and up
            Double_t B = Hmet_pt_Zjets->Integral(bin, Hmet_pt_Zjets->GetSize());
            Double_t S = Hmet_pt_sigQCD_new->Integral(bin, Hmet_pt_sigQCD_new->GetSize()) + Hmet_pt_sigEWK_new->Integral(bin, Hmet_pt_sigEWK_new->GetSize());

            if (B > 0 && S > 0)
            {
              Double_t Z_bin = sqrt(2 * ((S + B) * log(1 + (S / B)) - S));
              // cout << "B=" << B << ", S=" << S << ", Z_bin=" << Z_bin << endl;
              Hmet_pt_signal_sing->SetBinContent(bin, Z_bin);

              if (Z_bin > maxsig)
              {
                a = i;
                b = j;
                c = k;
                z = l;
                maxsig = Z_bin;
                cout << "The maximum Significance is: " << maxsig << endl;
                // printf("Found for bin = %f\n",bin);
                printf("So met_tst was: %f\n", xbins[bin - 1]);
                printf("Combination i = %d, j = %d, k = %d, l = %d\n", i, j, k, l);
              }
            }
          }

          TH1F *h01 = (TH1F *)Hmet_pt_signal_sing->Clone("h01");
          h01->SetLineColor(4);
          h01->SetMarkerColor(4);
          h01->Draw("E0X0");

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

          TCanvas *c2 = new TCanvas("c2", "MET significance", 0., 0., 600, 600);

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

          delete Hmet_pt_sigQCD;
          delete Hmet_pt_sigEWK;
          delete Hmet_pt_data;
          delete Hmet_pt_Zjets;
          delete Hmet_pt_WZ;
          delete Hmet_pt_WW;
          delete Hmet_pt_Wt;
          delete Hmet_pt_tt;
          delete Hmet_pt_trib;
          delete Hmet_pt_othr;
          delete Hmet_pt_signal_sing;

          delete Hmet_pt_signif_sigQCD;
          delete Hmet_pt_signif_sigEWK;
          delete Hmet_pt_signif_data;
          delete Hmet_pt_signif_Zjets;
          delete Hmet_pt_signif_WZ;
          delete Hmet_pt_signif_WW;
          delete Hmet_pt_signif_Wt;
          delete Hmet_pt_signif_tt;
          delete Hmet_pt_signif_trib;
          delete Hmet_pt_signif_othr;

          delete Hmet_pt_signif_signal_sing;
          delete Hmet_pt_signif_signal_pur;
          
          metritiri += 1;
          Printf("TO METRITIRIIIIII # %d #\n", metritiri);
        }
      }
    }
  }

  printf("Maximum significance is %f and found for combination of i = %d, j = %d, k = %d, l = %d \n", maxsig, a, b, c, z);
  printf("The best cuts that can be applied are: DLepR < %f, DMetZPhi > %f, met_tst > %f, MetOHT > %f \n", testDLepR[a], testdMetZPhi[b], testmet_tst[c], testMetOHT[z]);
  return 0;
}