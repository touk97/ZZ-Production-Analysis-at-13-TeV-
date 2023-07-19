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
//
//
//
//
// MAIN
void zjets_splitter()
{
     vector<string> filenames = {"DATA", "WZ", "Z_jets_ee", "Z_jets_mumu",
                                 "top", "ttbarV_ttbarVV", "Wt", "WW",
                                 "llll", "llqq", "VVV", "W_jets", "Ztt",
                                 "lllljj", "llvv", "llvvjj", "llvvjj_WW",
                                 "WZ_jj"};

     string filepath = "../data/SAMPLES/Zjets/";

     for (string &filename : filenames)
     {
          cout << endl
               << "   Creating file " << filename << "..." << endl
               << endl;

          Double_t signal = 0.0;
          Double_t weight, n_jets;
          Double_t M2Lep, met_tst, met_signif, dMetZPhi, MetOHT, dLepR;
          Double_t leading_pT_lepton, subleading_pT_lepton, Z_pT;
          Double_t n_bjets, detajj, mjj, leading_jet_pt, second_jet_pt;
          Double_t event_3CR, event_type;

          // Open the original TFile and TTree
          TFile *file = new TFile((string(filepath) + filename + ".root").c_str(), "READ");
          TTree *tree = (TTree *)file->Get("tree");

          // Create a new TFile and TTree
          TFile *file_output = new TFile(("../data/SAMPLES/Zjets3/" + filename + ".root").c_str(), "RECREATE");
          TTree *tree_output = tree->CloneTree(0); // Clone the structure of the tree

          // Retrieve the list of branches
          TObjArray *branches = tree->GetListOfBranches();
          tree->SetBranchAddress("n_jets", &n_jets);
          tree->SetBranchAddress("global_weight", &weight);
          tree->SetBranchAddress("M2Lep", &M2Lep);
          tree->SetBranchAddress("met_tst", &met_tst);
          tree->SetBranchAddress("met_signif", &met_signif);
          tree->SetBranchAddress("dMetZPhi", &dMetZPhi);
          tree->SetBranchAddress("MetOHT", &MetOHT);
          tree->SetBranchAddress("dLepR", &dLepR);
          tree->SetBranchAddress("leading_pT_lepton", &leading_pT_lepton);
          tree->SetBranchAddress("subleading_pT_lepton", &subleading_pT_lepton);
          tree->SetBranchAddress("Z_pT", &Z_pT);
          tree->SetBranchAddress("n_bjets", &n_bjets);
          tree->SetBranchAddress("detajj", &detajj);
          tree->SetBranchAddress("mjj", &mjj);
          tree->SetBranchAddress("leading_jet_pt", &leading_jet_pt);
          tree->SetBranchAddress("second_jet_pt", &second_jet_pt);
          tree->SetBranchAddress("event_3CR", &event_3CR);
          tree->SetBranchAddress("event_type", &event_type);

          // Create new branches for the desired parameters in the output TTree
          TBranch *branch_M2Lep = tree_output->Branch("M2Lep", &M2Lep, "M2Lep/D");
          TBranch *branch_met_tst = tree_output->Branch("met_tst", &met_tst, "met_tst/D");
          TBranch *branch_met_signif = tree_output->Branch("met_signif", &met_signif, "met_signif/D");
          TBranch *branch_dMetZPhi = tree_output->Branch("dMetZPhi", &dMetZPhi, "dMetZPhi/D");
          TBranch *branch_MetOHT = tree_output->Branch("MetOHT", &MetOHT, "MetOHT/D");
          TBranch *branch_dLepR = tree_output->Branch("dLepR", &dLepR, "dLepR/D");
          TBranch *branch_leading_pT_lepton = tree_output->Branch("leading_pT_lepton", &leading_pT_lepton, "leading_pT_lepton/D");
          TBranch *branch_subleading_pT_lepton = tree_output->Branch("subleading_pT_lepton", &subleading_pT_lepton, "subleading_pT_lepton/D");
          TBranch *branch_Z_pT = tree_output->Branch("Z_pT", &Z_pT, "Z_pT/D");
          TBranch *branch_n_bjets = tree_output->Branch("n_bjets", &n_bjets, "n_bjets/D");
          TBranch *branch_detajj = tree_output->Branch("detajj", &detajj, "detajj/D");
          TBranch *branch_mjj = tree_output->Branch("mjj", &mjj, "mjj/D");
          TBranch *branch_leading_jet_pt = tree_output->Branch("leading_jet_pt", &leading_jet_pt, "leading_jet_pt/D");
          TBranch *branch_second_jet_pt = tree_output->Branch("second_jet_pt", &second_jet_pt, "second_jet_pt/D");
          TBranch *branch_event_3CR = tree_output->Branch("event_3CR", &event_3CR, "event_3CR/D");
          TBranch *branch_event_type = tree_output->Branch("event_type", &event_type, "event_type/D");

          Int_t eventCounter = 0;
          const Int_t batchSize = 1000; // Adjust the batch size as needed

          for (Int_t i = 0; i < tree->GetEntriesFast(); i++)
          {
               tree->GetEntry(i);
               if (n_jets > 2)
               {
                    // Copy the desired event(s) from the original TTree to the new TTree
                    tree->GetEntry(i);
                    tree_output->Fill();
                    signal = signal + weight;

                    eventCounter++;

                    // Save the output TTree periodically
                    if (eventCounter % batchSize == 0)
                    {
                         tree_output->AutoSave();
                         cout << "   Processed " << eventCounter << " events." << endl;
                    }
               }
          }
          cout << "   TOTAL WEIGHT:   " << signal << endl;
          cout << "   TOTAL EVENTS:   " << tree->GetEntries() << endl;

          // Save the final output TTree
          tree_output->AutoSave();

          // Close the TFiles
          file_output->Close();
          file->Close();

          delete file_output;
          delete file;

          cout << "   -------------------------------------" << endl;
     }

     return;
}