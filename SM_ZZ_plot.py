import ROOT
from ROOT import TRandom3, TChain, TFile, TLorentzVector, TH1, TH2, TProfile, TMath, TStyle, TLatex, TPad, TLegend, TCanvas, TLine
import numpy as np
import math 
from math import sqrt, log
import array

def DoReco(tree, Hmet_pt, Hmet_pt_signif, MC):
    ROOT.TH1.SetDefaultSumw2(True)
    
    nentries = tree.GetEntries()
    
    M2Lep = 0.0 # Defines the flags needed
    met_tst = 0.0
    met_signif = 0.0
    dMetZPhi = 0.0
    frac_pT = 0.0
    MetOHT = 0.0
    dPhiJ100met = 0.0
    dLepR = 0.0
    n_bjets = 0
    n_jets = 0
    leading_pT_lepton = 0.0
    subleading_pT_lepton = 0.0
    detajj = 0.0
    mjj = 0.0
    leading_jet_pt = 0.0
    second_jet_pt = 0.0
    event_3CR = 0
    event_type = 0
    weight = 1.0
    
    
    tree.SetBranchAddress("M2Lep", M2Lep)
    tree.SetBranchAddress("met_tst", met_tst)
    tree.SetBranchAddress("met_signif", met_signif)
    tree.SetBranchAddress("dMetZPhi", dMetZPhi)
    tree.SetBranchAddress("frac_pT", frac_pT)
    tree.SetBranchAddress("MetOHT", MetOHT)
    tree.SetBranchAddress("dPhiJ100met", dPhiJ100met)
    tree.SetBranchAddress("dLepR", dLepR)
    tree.SetBranchAddress("leading_pT_lepton", leading_pT_lepton)
    tree.SetBranchAddress("subleading_pT_lepton", subleading_pT_lepton)
    tree.SetBranchAddress("n_jets", n_jets)
    tree.SetBranchAddress("n_bjets", n_bjets)
    tree.SetBranchAddress("detajj", detajj)
    tree.SetBranchAddress("mjj", mjj)
    tree.SetBranchAddress("leading_jet_pt", leading_jet_pt)
    tree.SetBranchAddress("second_jet_pt", second_jet_pt)
    tree.SetBranchAddress("event_3CR", event_3CR)
    tree.SetBranchAddress("event_type", event_type)
  

    signal = 0.
    signaler = 0.
    Norm = 1.
    if MC == 1:
        tree.SetBranchAddress("global_weight", weight)
        # Loop over events
    for i in range(nentries):
        tree.GetEntry(i)
    
        # Define signal region - Change as needed
        #
        # Inclusive selection:
        # ---------------------
        # RECO cut-based:
        # if (event_3CR==0 && (event_type==0 || event_type==1) && leading_pT_lepton>30 && subleading_pT_lepton>20 && M2Lep>80 && M2Lep<100 && n_bjets<1 && dLepR<1.8 && dMetZPhi>2.7 && met_tst > 110 && MetOHT>0.65)
        # if (event_3CR==0 && (event_type==0 || event_type==1) && leading_pT_lepton>30 && subleading_pT_lepton>20 && M2Lep>80 && M2Lep<100 && n_bjets<1 && dLepR<1.8 && dMetZPhi>2.7 && met_tst > 110 && MetOHT>0.65 && met_signif>10)
        # if (event_3CR==0 && (event_type==0 || event_type==1) && leading_pT_lepton>30 && subleading_pT_lepton>20 && M2Lep>80 && M2Lep<100 && n_bjets<1 && dLepR<1.8 && dMetZPhi>2.7 && met_tst > 110 && MetOHT>0.65 && met_signif>0) # met_tst/sqrt(met_tst/MetOHT)>10)
        # FIDUCIAL
        # if (event_3CR==0 && (event_type==0 || event_type==1) && leading_pT_lepton>30 && subleading_pT_lepton>20 && M2Lep>80 && M2Lep<100 && n_bjets<1 && dLepR<1.8 && dMetZPhi>2.7 && met_tst > 110 && MetOHT>0.65)
        if (event_3CR==0 and (event_type==0 or event_type==1) and leading_pT_lepton>30 and subleading_pT_lepton>20 and M2Lep>76 and M2Lep<116 and n_bjets<1 and dLepR<1.9 and dMetZPhi>2.6 and met_tst > 90 and MetOHT>0.6):
    
        # FID-BDT:
        # if (event_3CR==0 && (event_type==0 || event_type==1) &&  leading_pT_lepton>30 && subleading_pT_lepton>20 && M2Lep>76 && M2Lep<106 && n_bjets<1 && dLepR<2.2 && dMetZPhi>1.3 && met_tst > 70 && MetOHT>0.3 && met_signif>0) # met_tst/sqrt(met_tst/MetOHT)>10)
        # FID-DNN:
        # if (event_3CR==0 && (event_type==0 || event_type==1) &&  leading_pT_lepton>30 && subleading_pT_lepton>20 && M2Lep>76 && M2Lep<106 && n_bjets<1 && dLepR<2.2 && dMetZPhi>1.3 && met_tst > 90 && MetOHT>0.1 && met_signif>0) # met_tst/sqrt(met_tst/MetOHT)>10)
        #
        # VBS selection:
        # --------------
        # if (event_3CR==0 && (event_type==0 || event_type==1) && met_tst>120 &&  met_signif>0 && dLepR<1.8 && dMetZPhi>2.5 && n_bjets<1 && MetOHT>0.6 && n_jets > 1)
        # if (event_3CR==0 && (event_type==0 || event_type==1) && met_tst>90 && met_signif>12 && dLepR<1.5 && dMetZPhi>2.6 && n_bjets<1 && MetOHT>0.6 )
        #(e.g., here: not an event in the 3-lepton Control-Region, event_type==0 or 1, MissingET > 90 GeV)
        # if (event_3CR==0 && (event_type==0 || event_type==1) && met_tst>90)
        # Nominal:
        # if (event_3CR==0 && (event_type==0 || event_type==1) &&  leading_pT_lepton>30  && subleading_pT_lepton>20 && M2Lep>80 && M2Lep<100 && n_jets>1 && n_bjets<1 && leading_jet_pt>30 && second_jet_pt>30 && detajj>1  && mjj>100 && dMetZPhi>2.2 && dLepR<1.8 && met_tst>150 && MetOHT>0.4)
        # Nominal + METSignif :
        # if (event_3CR==0 && (event_type==0 || event_type==1) &&  leading_pT_lepton>30  && subleading_pT_lepton>20 && M2Lep>80 && M2Lep<100 && n_jets>1 && n_bjets<1 && leading_jet_pt>30 && second_jet_pt>30 && detajj>1  && mjj>100 && dMetZPhi>2.2 && dLepR<1.8 && met_tst>150 && MetOHT>0.4 && met_signif>10 )
        # Inclusive + di-jet variables:
        # if (event_3CR==0 && (event_type==0 || event_type==1) &&  leading_pT_lepton>30  && subleading_pT_lepton>20 && M2Lep>80 && M2Lep<100 && n_jets>1 && n_bjets<1 && leading_jet_pt>30 && second_jet_pt>30 && detajj>1  && mjj>100 && dMetZPhi>2.2 && dLepR<1.8 && met_tst>110 && MetOHT>0.65 && met_signif>0 )
          
            signal += weight  # signal yield is sum of weights
            signaler += weight * weight  # keep Sum(weights^2) for calculating the error on the signal yield
            Hmet_pt.Fill(met_tst, weight)  # Fills histogram but with weight
            Hmet_pt_signif.Fill(met_signif, weight)
    print("N =", signal, "+-", math.sqrt(signaler))
    # Store signal yield and error in a list
    events = []
    events.append(signal)
    events.append(math.sqrt(signaler))
    
    return events
ROOT.gStyle.SetLegendFont(42)
ROOT.gStyle.SetPalette(1)
#TString sample = "/Results/trial2.root";
# ---convert TString to string
#std::string sample_string(sample.Data());
# ---convert string to TString
#TString sample_TString(sample_string)

# Preliminary minitrees for ZZllvvjj analysis
file_data = ROOT.TFile("data/data_fullRun2_MET90.root")
tree_data = file_data.Get("tree_PFLOW")
file_signalQCD = ROOT.TFile("data/mc16ade_QCD_EWKcor_ZZllvv_MET90.root") #QCD signal with EWK corrections
tree_signalQCD = file_signalQCD.Get("tree_PFLOW")
file_signalEWK = ROOT.TFile("data/mc16ade_EWK_ZZllvvjj_MET90.root") # EW signal
tree_signalEWK = file_signalEWK.Get("tree_PFLOW")
file_Zjets = ROOT.TFile("data/mc16ade_Zjets_MET90.root")
tree_Zjets = file_Zjets.Get("tree_PFLOW")

file_WZ = ROOT.TFile("data/mc16ade_WZ_MET90.root")
tree_WZ = file_WZ.Get("tree_PFLOW")

file_WW = ROOT.TFile("data/mc16ade_WW_Pow_MET90.root")
tree_WW = file_WW.Get("tree_PFLOW")

file_Wt = ROOT.TFile("data/mc16ade_Wt_MET90.root")
tree_Wt = file_Wt.Get("tree_PFLOW")

file_tt = ROOT.TFile("data/mc16ade_TOP_MET90.root")
tree_tt = file_tt.Get("tree_PFLOW")

file_trib = ROOT.TFile("data/mc16ade_VVV_MET90.root")
tree_trib = file_trib.Get("tree_PFLOW")

file_othr = ROOT.TFile("data/mc16ade_ttV_ttVV_MET90.root")
tree_othr = file_othr.Get("tree_PFLOW")



ROOT.TH1.SetDefaultSumw2(ROOT.kTRUE)
xbins = [90, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 345, 375, 410, 450, 500, 580, 700, 1000]

Hmet_pt_sigQCD = ROOT.TH1F("Hmet_pt_sigQCD", "Hmet_pt_sigQCD", 20, array.array('f', xbins))
Hmet_pt_sigEWK = ROOT.TH1F("Hmet_pt_sigEWK", "Hmet_pt_sigEWK", 20, array.array('f', xbins))
Hmet_pt_data = ROOT.TH1F("Hmet_pt_data", "Hmet_pt_data", 20, array.array('f', xbins))
Hmet_pt_Zjets = ROOT.TH1F("Hmet_pt_Zjets", "Hmet_pt_Zjets", 20, array.array('f', xbins))
Hmet_pt_WZ = ROOT.TH1F("Hmet_pt_WZ", "Hmet_pt_WZ", 20, array.array('f', xbins))
Hmet_pt_WW = ROOT.TH1F("Hmet_pt_WW", "Hmet_pt_WW", 20, array.array('f', xbins))
Hmet_pt_Wt = ROOT.TH1F("Hmet_pt_Wt", "Hmet_pt_Wt", 20, array.array('f', xbins))
Hmet_pt_tt = ROOT.TH1F("Hmet_pt_tt", "Hmet_pt_tt", 20, array.array('f', xbins))
Hmet_pt_trib = ROOT.TH1F("Hmet_pt_trib", "Hmet_pt_trib", 20, array.array('f', xbins))
Hmet_pt_othr = ROOT.TH1F("Hmet_pt_othr", "Hmet_pt_othr", 20, array.array('f', xbins))
Hmet_pt_signal_sing = ROOT.TH1F("Hmet_pt_signal_sing", "Hmet_pt_signal_sing", 20, array.array('f', xbins))
Hmet_pt_signif_sigQCD = ROOT.TH1F("Hmet_pt_signif_sigQCD", "Hmet_pt_signif_sigQCD", 35, 0, 35)
Hmet_pt_signif_sigEWK = ROOT.TH1F("Hmet_pt_signif_sigEWK", "Hmet_pt_signif_sigEWK", 35, 0, 35)
Hmet_pt_signif_data = ROOT.TH1F("Hmet_pt_signif_data", "Hmet_pt_signif_data", 35, 0, 35)
Hmet_pt_signif_Zjets = ROOT.TH1F("Hmet_pt_signif_Zjets", "Hmet_pt_signif_Zjets", 35, 0, 35)
Hmet_pt_signif_WZ = ROOT.TH1F("Hmet_pt_signif_WZ", "Hmet_pt_signif_WZ", 35, 0, 35)
Hmet_pt_signif_WW = ROOT.TH1F("Hmet_pt_signif_WW", "Hmet_pt_signif_WW", 35, 0, 35)
Hmet_pt_signif_Wt = ROOT.TH1F("Hmet_pt_signif_Wt", "Hmet_pt_signif_Wt", 35, 0, 35)
Hmet_pt_signif_tt = ROOT.TH1F("Hmet_pt_signif_tt", "Hmet_pt_signif_tt", 35, 0, 35)
Hmet_pt_signif_trib = ROOT.TH1F("Hmet_pt_signif_trib", "Hmet_pt_signif_trib", 35, 0, 35)
Hmet_pt_signif_othr = ROOT.TH1F("Hmet_pt_signif_othr", "Hmet_pt_signif_othr", 35, 0, 35)

Hmet_pt_signif_signal_sing = ROOT.TH1F("Hmet_pt_signif_signal_sing", "Hmet_pt_signif_signal_sing", 35, 0, 35)
Hmet_pt_signif_signal_pur = ROOT.TH1F("Hmet_pt_signif_signal_pur", "Hmet_pt_signif_signal_pur", 35, 0, 35)
events_data = []
events_NsigQCD = []
events_NsigEWK = []
events_NZjets = []
events_NWZ = []
events_NWW = []
events_NWt = []
events_Ntt = []
events_Ntrib = []
events_Nothr = []

# Event Yields
print("===Data===")
events_data = DoReco(tree_data, Hmet_pt_data, Hmet_pt_signif_data, 0)  # DoReco Fills the vectors
print("===Signal QCD ZZ===")
events_NsigQCD = DoReco(tree_signalQCD, Hmet_pt_sigQCD, Hmet_pt_signif_sigQCD, 1)
print("===Signal EWK ZZ===")
events_NsigEWK = DoReco(tree_signalEWK, Hmet_pt_sigEWK, Hmet_pt_signif_sigEWK, 1)
print("===Background===")
print("Zjets")
events_NZjets = DoReco(tree_Zjets, Hmet_pt_Zjets, Hmet_pt_signif_Zjets, 1)
print("WZ")
events_NWZ = DoReco(tree_WZ, Hmet_pt_WZ, Hmet_pt_signif_WZ, 1)
print("tt")
events_Ntt = DoReco(tree_tt, Hmet_pt_tt, Hmet_pt_signif_tt, 1)
print("WW")
events_NWW = DoReco(tree_WW, Hmet_pt_WW, Hmet_pt_signif_WW, 1)
print("Wt")
events_NWt = DoReco(tree_Wt, Hmet_pt_Wt, Hmet_pt_signif_Wt, 1)
print("trib")
events_Ntrib = DoReco(tree_trib, Hmet_pt_trib, Hmet_pt_signif_trib, 1)
print("othr")
events_Nothr = DoReco(tree_othr, Hmet_pt_othr, Hmet_pt_signif_othr, 1)
totalBKG = (
events_NZjets[0]
+ events_Ntt[0]
+ events_NWt[0]
+ events_NWW[0]
+ events_NWZ[0]
+ events_Ntrib[0]
+ events_Nothr[0]
)
totalBKG_er = sqrt(
    events_NZjets[1] ** 2
    + events_Ntt[1] ** 2
    + events_NWt[1] ** 2
    + events_NWW[1] ** 2
    + events_NWZ[1] ** 2
    + events_Ntrib[1] ** 2
    + events_Nothr[1] ** 2
)
print("total Bkg =", totalBKG, "+-", totalBKG_er)
print("S/B =", (events_NsigQCD[0] + events_NsigEWK[0]) / totalBKG)
S = events_NsigQCD[0] + events_NsigEWK[0]
B = totalBKG
Z = sqrt(2 * ((S + B) * log(1 + (S / B)) - S))
print("Significance =", Z)
Hmet_pt_sigQCD_new = Hmet_pt_sigQCD.Clone("Hmet_pt_sigQCD_new")
Hmet_pt_sigEWK_new = Hmet_pt_sigEWK.Clone("Hmet_pt_sigEWK_new")

c1 = ROOT.TCanvas("c1", "MET", 1400, 600)
pad_met = ROOT.TPad("pad_met", "This is pad_met", 0.01, 0.30, 1.0, 1.0)
pad_met2 = ROOT.TPad("pad_met2", "This is pad_met2", 0.01, 0.01, 1.0, 0.30)

pad_met.SetBorderSize(0)
pad_met.SetBottomMargin(0.02)
pad_met.Draw()
pad_met2.SetBottomMargin(0.35)
pad_met2.SetTopMargin(0.0)
pad_met2.SetBorderSize(0)
pad_met2.Draw()
pad_met.cd()
Hmet_pt_sigQCD.SetFillColor(38)
Hmet_pt_sigEWK.SetFillColor(6)
Hmet_pt_Zjets.SetFillColor(8)
Hmet_pt_WZ.SetFillColor(4)
Hmet_pt_WW.SetFillColor(5)
Hmet_pt_Wt.SetFillColor(7)
Hmet_pt_tt.SetFillColor(42)
Hmet_pt_othr.SetFillColor(30)
Hmet_pt_trib.SetFillColor(46)

Hmet_pt_trib.Add(Hmet_pt_othr)
Hmet_pt_Wt.Add(Hmet_pt_trib)
Hmet_pt_WW.Add(Hmet_pt_Wt)
Hmet_pt_WZ.Add(Hmet_pt_WW)
Hmet_pt_tt.Add(Hmet_pt_WZ)
Hmet_pt_Zjets.Add(Hmet_pt_tt)
Hmet_pt_sigEWK.Add(Hmet_pt_Zjets)
Hmet_pt_sigQCD.Add(Hmet_pt_sigEWK)

Hmet_pt_sigQCD.Draw("hist")
Hmet_pt_sigEWK.Draw("histsame")
Hmet_pt_Zjets.Draw("histsame")
Hmet_pt_tt.Draw("histsame")
Hmet_pt_WZ.Draw("histsame")
Hmet_pt_WW.Draw("histsame")
Hmet_pt_Wt.Draw("histsame")
Hmet_pt_trib.Draw("histsame")
Hmet_pt_othr.Draw("histsame")
Hmet_pt_data.Draw("sameE0X0")

Hmet_pt_sigQCD.GetYaxis().SetTitleSize(0.06)
Hmet_pt_sigQCD.GetYaxis().SetTitleOffset(1.11)
Hmet_pt_sigQCD.GetYaxis().SetLabelSize(0.05)
Hmet_pt_sigQCD.GetXaxis().SetLabelSize(0.0)
Hmet_pt_sigQCD.GetYaxis().SetTitle("Events")
tex = ROOT.TLatex(0.2, 0.8, "#intL dt = 138.9 fb^{-1}, #sqrt{s} = 13 TeV, ee+#mu#mu, mc16ade")
tex.SetNDC()
tex.SetTextFont(42)
tex.SetLineWidth(2)
tex.Draw()

leg = ROOT.TLegend(0.6, 0.3, 0.9, 0.75, None, "brNDC")
leg_entry = leg.AddEntry(Hmet_pt_data, "Data", "lp")
leg_entry = leg.AddEntry(Hmet_pt_sigQCD, "ZZQCD", "f")
leg_entry = leg.AddEntry(Hmet_pt_sigEWK, "ZZEWK", "f")
leg_entry = leg.AddEntry(Hmet_pt_Zjets, "Zjets", "f")
leg_entry = leg.AddEntry(Hmet_pt_tt, "TOP", "f")
leg_entry = leg.AddEntry(Hmet_pt_WZ, "WZ", "f")
leg_entry = leg.AddEntry(Hmet_pt_WW, "WW", "f")
leg_entry = leg.AddEntry(Hmet_pt_Wt, "Wt", "f")
leg_entry = leg.AddEntry(Hmet_pt_trib, "VVV", "f")
leg_entry = leg.AddEntry(Hmet_pt_othr, "other", "f")
leg.SetLineColor(0)
leg.SetBorderSize(0)
leg.Draw()

pad_met.RedrawAxis()
pad_met2.cd()
Z_bin = 0
for bin in range(1, Hmet_pt_Zjets.GetSize()):
    B = Hmet_pt_Zjets.GetBinContent(bin)
    S = Hmet_pt_sigQCD_new.GetBinContent(bin) + Hmet_pt_sigEWK_new.GetBinContent(bin)

    if B > 0 and S > 0:
        Z_bin = sqrt(2 * ((S + B) * log(1 + (S / B)) - S))
        Hmet_pt_signal_sing.SetBinContent(bin, Z_bin)

h01 = ROOT.TH1F(Hmet_pt_signal_sing.Clone("h01"))
h01.SetLineColor(4)
h01.SetMarkerColor(4)
h01.Draw("E0X0")
# Plot Data/MC ratio

# h00 = Hmet_pt_sigQCD.Clone("h00")
# h01 = Hmet_pt_data.Clone("h01")
# h00.Sumw2()
# h01.Sumw2()
# h01.Divide(h00)
# h01.Draw("E0X0")

# line = ROOT.TLine(90, 1, 1000, 1)
# line.SetLineStyle(2)
# line.SetLineWidth(2)
# line.Draw("same")
h01.GetXaxis().SetTitleSize(0.15)
h01.GetXaxis().SetTitle("E_{T}^{miss} [GeV]")
h01.GetXaxis().SetTitleOffset(1.1)
h01.GetXaxis().SetLabelSize(0.13)

h01.GetYaxis().SetTitleOffset(0.45)
h01.GetYaxis().SetTitleSize(0.135)
h01.GetYaxis().SetLabelSize(0.11)
#h01.GetYaxis().SetTitle("#frac{Data}{MC}")
h01.GetYaxis().SetTitle("Sig. Significance")

Hmet_pt_signif_sigEWK_new = Hmet_pt_signif_sigEWK.Clone("Hmet_pt_signif_sigEWK_new")
Hmet_pt_signif_sigQCD_new = Hmet_pt_signif_sigQCD.Clone("Hmet_pt_signif_sigQCD_new")

c2 = ROOT.TCanvas("c2", "MET significance", 1400, 600)
pad_met_signif = ROOT.TPad("pad_met_signif", "This is pad_met_signif", 0.01, 0.30, 1.0, 1.0)
pad_met_signif2 = ROOT.TPad("pad_met_signif2", "This is pad_met_signif2", 0.01, 0.01, 1.0, 0.30)

pad_met_signif.SetBorderSize(0)
pad_met_signif.SetBottomMargin(0.02)
pad_met_signif.Draw()
pad_met_signif2.SetBottomMargin(0.35)
pad_met_signif2.SetTopMargin(0.0)
pad_met_signif2.SetBorderSize(0)
pad_met_signif2.Draw()
pad_met_signif.cd()
Hmet_pt_signif_sigQCD.SetFillColor(38)
Hmet_pt_signif_sigEWK.SetFillColor(6)
Hmet_pt_signif_Zjets.SetFillColor(8)
Hmet_pt_signif_WZ.SetFillColor(4)
Hmet_pt_signif_WW.SetFillColor(5)
Hmet_pt_signif_Wt.SetFillColor(7)
Hmet_pt_signif_tt.SetFillColor(42)
Hmet_pt_signif_othr.SetFillColor(30)
Hmet_pt_signif_trib.SetFillColor(46)

Hmet_pt_signif_trib.Add(Hmet_pt_signif_othr)
Hmet_pt_signif_Wt.Add(Hmet_pt_signif_trib)
Hmet_pt_signif_WW.Add(Hmet_pt_signif_Wt)
Hmet_pt_signif_WZ.Add(Hmet_pt_signif_WW)
Hmet_pt_signif_tt.Add(Hmet_pt_signif_WZ)
Hmet_pt_signif_Zjets.Add(Hmet_pt_signif_tt)
Hmet_pt_signif_sigEWK.Add(Hmet_pt_signif_Zjets)
Hmet_pt_signif_sigQCD.Add(Hmet_pt_signif_sigEWK)

Hmet_pt_signif_sigQCD.Draw("hist")
Hmet_pt_signif_sigEWK.Draw("histsame")
Hmet_pt_signif_Zjets.Draw("histsame")
Hmet_pt_signif_tt.Draw("histsame")
Hmet_pt_signif_WZ.Draw("histsame")
Hmet_pt_signif_WW.Draw("histsame")
Hmet_pt_signif_Wt.Draw("histsame")
Hmet_pt_signif_trib.Draw("histsame")
Hmet_pt_signif_othr.Draw("histsame")
Hmet_pt_signif_data.Draw("sameE0X0")

Hmet_pt_signif_sigQCD.GetYaxis().SetTitleSize(0.04)
Hmet_pt_signif_sigQCD.GetYaxis().SetTitleOffset(1.)
Hmet_pt_signif_sigQCD.GetYaxis().SetLabelSize(0.05)
Hmet_pt_signif_sigQCD.GetXaxis().SetLabelSize(0.00)
Hmet_pt_signif_sigQCD.GetYaxis().SetTitle("Events")

tex = ROOT.TLatex(0.2, 0.8, "#intL dt = 138.9 fb^{-1}, #sqrt{s} = 13 TeV, ee+#mu#mu, mc16ade")
tex.SetNDC()
tex.SetTextFont(42)
tex.SetLineWidth(2)
tex.Draw()

leg = ROOT.TLegend(0.6, 0.3, 0.7, 0.75, None, "brNDC")
# leg_entry
leg_entry = leg.AddEntry(Hmet_pt_signif_data, "Data", "lp")
leg_entry = leg.AddEntry(Hmet_pt_signif_sigQCD, "ZZQCD", "f")
leg_entry = leg.AddEntry(Hmet_pt_signif_sigEWK, "ZZEWK", "f")
leg_entry = leg.AddEntry(Hmet_pt_signif_Zjets, "Zjets", "f")
leg_entry = leg.AddEntry(Hmet_pt_signif_tt, "TOP", "f")
leg_entry = leg.AddEntry(Hmet_pt_signif_WZ, "WZ", "f")
leg_entry = leg.AddEntry(Hmet_pt_signif_WW, "WW", "f")
leg_entry = leg.AddEntry(Hmet_pt_signif_Wt, "Wt", "f")
leg_entry = leg.AddEntry(Hmet_pt_signif_trib, "VVV", "f")
leg_entry = leg.AddEntry(Hmet_pt_signif_othr, "other", "f")
leg.SetLineColor(0)
leg.SetBorderSize(0)
leg.Draw()

pad_met_signif.RedrawAxis()
pad_met_signif2.cd()
# to plot significance
# Z_bin = 0
# for bin in range(1, Hmet_pt_signif_Zjets.GetSize()):
#     B = Hmet_pt_signif_Zjets.GetBinContent(bin)
#     S = Hmet_pt_signif_sigQCD_new.GetBinContent(bin) + Hmet_pt_signif_sigEWK_new.GetBinContent(bin)

#     if B > 0 and S > 0:
#         Z_bin = math.sqrt(2 * ((S + B) * math.log(1 + (S / B)) - S))
#         # print("B =", B, ", S =", S, ", Z_bin =", Z_bin)
#         Hmet_pt_signif_signal_sing.SetBinContent(bin, Z_bin)

# h01 = Hmet_pt_signif_signal_sing.Clone("h01")
# h01.SetLineColor(4)
# h01.SetMarkerColor(4)
# h01.Draw("E0X0")
# do plot Data/MC
h00 = Hmet_pt_signif_sigQCD.Clone("h00")
h01 = Hmet_pt_signif_data.Clone("h01")
h00.Sumw2()
h01.Sumw2()
h01.Divide(h00)
h01.Draw("E0X0")

line = ROOT.TLine(0, 1, 35, 1)
line.SetLineStyle(2)
line.SetLineWidth(2)
line.Draw("same")

h01.GetXaxis().SetTitleSize(0.15)
h01.GetXaxis().SetTitle("E_{T}^{miss} significance")
h01.GetXaxis().SetTitleOffset(1.1)
h01.GetXaxis().SetLabelSize(0.13)

h01.GetYaxis().SetTitleOffset(0.45)
h01.GetYaxis().SetTitleSize(0.135)
h01.GetYaxis().SetLabelSize(0.11)
h01.GetYaxis().SetTitle("#frac{Data}{MC}")
# h01.GetYaxis().SetTitle("Sig. Significance")

c1.Draw()
c1.SaveAs("./figs/METadsfa.png")
c2.SaveAs("./figs/MET_S.png")


