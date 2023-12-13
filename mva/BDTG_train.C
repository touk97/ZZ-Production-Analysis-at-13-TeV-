

// To use TMVA GUI
//  root[0] TMVA::TMVAGui("TMVA.root")

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
// #include "../tmva/test/TMVAGui.C"

#if not defined(__CINT__) || defined(__MAKECINT__)
// needs to be included when makecint runs (ACLIC)
// #include "TMVA/Factory.h"
// #include "TMVA/Tools.h"

#include "TROOT.h"
#include "TString.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"

using namespace std;

#endif

//
// CUSTOM STREAM BUFFER CLASS FOR OUTPUT LOG FILE
class DualStreamBuffer : public streambuf
{
public:
   // Constructor
   DualStreamBuffer(streambuf *primaryBuffer, streambuf *secondaryBuffer)
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
   streambuf *primaryBuffer;   // Primary stream buffer (cout)
   streambuf *secondaryBuffer; // Secondary stream buffer (log file)
};
// END OF LOG FILE CLASS
//

void BDTG_train(TString myMethodList = "")
{

   // Output log file
   // ofstream logFile("./train/BDTG_train_0jets.txt");
   ofstream logFile("./train/BDTG_train.txt");

   DualStreamBuffer dualBuffer(cout.rdbuf(), logFile.rdbuf());

   streambuf *oldBuffer = cout.rdbuf(&dualBuffer);

   //---------------------------------------------------------------

   TMVA::Tools::Instance();

   // Default MVA methods to be trained + tested
   map<string, int> Use;

   //
   // --- Boosted Decision Trees
   Use["BDT"] = 0;  // uses Adaptive Boost
   Use["BDTG"] = 1; // uses Gradient Boost
   Use["BDTB"] = 0; // uses Bagging
   Use["BDTD"] = 0; // decorrelation + Adaptive Boost
   Use["BDTF"] = 0; // allow usage of fisher discriminant for node splitting

   // ---------------------------------------------------------------

   cout << endl;
   cout << "==> Start TMVAClassification" << endl;

   if (myMethodList != "")
   {
      for (map<string, int>::iterator it = Use.begin(); it != Use.end(); it++)
         it->second = 0;

      vector<TString> mlist = TMVA::gTools().SplitString(myMethodList, ',');
      for (UInt_t i = 0; i < mlist.size(); i++)
      {
         string regMethod(mlist[i]);

         if (Use.find(regMethod) == Use.end())
         {
            cout << "Method \"" << regMethod << "\" not known in TMVA under this name. Choose among the following:" << endl;
            for (map<string, int>::iterator it = Use.begin(); it != Use.end(); it++)
               cout << it->first << " ";
            cout << endl;
            return;
         }
         Use[regMethod] = 1;
      }
   }

   // --------------------------------------------------------------------------------------------------

   // TString filename( "./train/BDTG_train_0jets.root" );  //Change for jet multiplicity
   TString filename("./train/BDTG_train.root");
   TFile *output = TFile::Open(filename, "RECREATE");

   // Create the factory object.

   TMVA::Factory factory("TMVAClassification", output,
                         "!V:ROC:!Correlations:!Silent:Color:DrawProgressBar:AnalysisType=Classification");
   TMVA::DataLoader loader("dataset");

   loader.AddVariable("met_tst", 'F');
   loader.AddVariable("met_signif", 'F');
   loader.AddVariable("mT_ZZ", 'F');
   loader.AddVariable("dLepR", 'F');
   loader.AddVariable("leading_pT_lepton", 'F');
   loader.AddVariable("subleading_pT_lepton", 'F');
   loader.AddVariable("Z_pT", 'F');
   loader.AddVariable("frac_pT", 'F');
   loader.AddVariable("sumpT_scalar", 'F');
   loader.AddVariable("MetOHT", 'F');
   // loader.AddVariable( "RhoZ := Z_pT/(leading_pT_lepton+subleading_pT_lepton)", 'F' );
   loader.AddVariable("dLepPhi:=fabs(lepplus_phi-lepminus_phi)", 'F');
   loader.AddVariable("dMetZPhi", 'F');
   loader.AddVariable("LepRatio := subleading_pT_lepton/leading_pT_lepton", 'F');
   loader.AddVariable("dLepEta:=fabs(lepplus_eta-lepminus_eta)", 'F');
   loader.AddVariable("M2Lep", 'F');
   // loader.AddVariable( "n_jets", 'I' );

   // loader.AddVariable( "Z_eta", 'F' );

   //   loader.AddVariable( "mTlmet", 'F' );

   //    loader.AddVariable( "n_bjets", 'I' );
   //  loader.AddVariable( "mT_Hinv", 'F' );
   //  loader.AddVariable( "ZpTomT_Hinv", 'F' );
   //  loader.AddVariable( "met_RefEle", 'F' );

   //   loader.AddVariable( "leading_jet_pt", 'F' );
   //   loader.AddVariable( "leading_jet_eta", 'F' );
   //   loader.AddVariable( "max_mjj", 'F' );
   //   loader.AddVariable( "max_detajj", 'F' );

   // loader.AddSpectator("event_3CR",'F');
   // //loader.AddSpectator("event_4CR",'F');
   // //loader.AddSpectator("SR_HM_LM",'I');
   // //loader.AddSpectator("met_tst",'F');
   // //loader.AddSpectator("met_signif",'F');
   // loader.AddSpectator( "n_bjets", 'I' );
   // loader.AddSpectator( "event_type", 'I' );

   

   TString filepath_signal = "/home/touk/Desktop/touk/physics/master/thesis/data/mc16ade_QCD_EWKcor_ZZllvv_MET90.root"; // QCD signal with EWK corrections
   TString filepath_Zjets = " /home/touk/Desktop/touk/physics/master/thesis/data/mc16ade_Zjets_MET90.root";
   TString filepath_WZ = "/home/touk/Desktop/touk/physics/master/thesis/data/mc16ade_WZ_MET90.root";
   TString filepath_WW = "/home/touk/Desktop/touk/physics/master/thesis/data/mc16ade_WW_Pow_MET90.root";
   TString filepath_Wt = "/home/touk/Desktop/touk/physics/master/thesis/data/mc16ade_Wt_MET90.root";
   TString filepath_tt = "/home/touk/Desktop/touk/physics/master/thesis/data/mc16ade_TOP_MET90.root";
   TString filepath_othr = "/home/touk/Desktop/touk/physics/master/thesis/data/mc16ade_ttV_ttVV_MET90.root";
   TString filepath_trib = "/home/touk/Desktop/touk/physics/master/thesis/data/mc16ade_VVV_MET90.root";

   TFile *file_signal = TFile::Open(filepath_signal);
   TFile *file_Zjets = TFile::Open(filepath_Zjets);
   TFile *file_WZ = TFile::Open(filepath_WZ);
   TFile *file_WW = TFile::Open(filepath_WW);
   TFile *file_Wt = TFile::Open(filepath_Wt);
   TFile *file_tt = TFile::Open(filepath_tt);
   TFile *file_othr = TFile::Open(filepath_othr);
   TFile *file_trib = TFile::Open(filepath_trib);

   cout << "--- TMVAClassification       : Using input file: " << file_signal->GetName() << endl;

   
   TTree *tree_signal = (TTree *)file_signal->Get("tree_PFLOW");
   TTree *tree_Zjets = (TTree *)file_Zjets->Get("tree_PFLOW");
   TTree *tree_WZ = (TTree *)file_WZ->Get("tree_PFLOW");
   TTree *tree_WW = (TTree *)file_WW->Get("tree_PFLOW");
   TTree *tree_Wt = (TTree *)file_Wt->Get("tree_PFLOW");
   TTree *tree_tt = (TTree *)file_tt->Get("tree_PFLOW");
   TTree *tree_trib = (TTree *)file_trib->Get("tree_PFLOW");
   TTree *tree_othr = (TTree *)file_othr->Get("tree_PFLOW");

   double Weight = 1; // already all normalized to lumi
   loader.AddSignalTree(tree_signal, Weight);
   loader.AddBackgroundTree(tree_Zjets, Weight);
   loader.AddBackgroundTree(tree_WZ, Weight);
   loader.AddBackgroundTree(tree_WW, Weight);
   loader.AddBackgroundTree(tree_Wt, Weight);
   loader.AddBackgroundTree(tree_tt, Weight);
   loader.AddBackgroundTree(tree_trib, Weight);
   loader.AddBackgroundTree(tree_othr, Weight);

   loader.SetSignalWeightExpression("global_weight");
   loader.SetBackgroundWeightExpression("global_weight");

   // Cuts
   // TCut SR_cut = "event_3CR==0 && (event_type==0 || event_type==1) && met_tst>70 && dLepR<1.8 && n_bjets==0 && dMetZPhi>2.2 && global_weight>0 && M2Lep > 80 && M2Lep <100 && leading_pT_lepton > 40 && subleading_pT_lepton > 20"; //
   TCut SR_cut = "event_3CR==0 && (event_type==0 || event_type==1) && MetOHT > 0.1 &&  met_tst>90 && dLepR<2.2 && n_bjets==0 && dMetZPhi>1.3 && M2Lep > 76 && M2Lep <106 && leading_pT_lepton > 30 && subleading_pT_lepton > 20"; //
   TCut jet_cut = " n_jets > -1";                                                                                                                                                                                                 // Change for jet multiplicity

   // If no numbers of events are given, half of the events in the tree are used
   // for training, and the other half for testing:
   loader.PrepareTrainingAndTestTree(SR_cut && jet_cut, "SplitMode=random:!V");

   // Method
   if (Use["BDTG"])
      factory.BookMethod(&loader, TMVA::Types::kBDT, "BDTG",
                         //   "!H:!V:NTrees=800:MinNodeSize=5%:BoostType=Grad::nCuts=20:MaxDepth=3");
                         // "!H:!V:NTrees=400:MinNodeSize=10%:BoostType=Grad::shrinkage=0.2");
                         "!H:!V:NTrees=400:MinNodeSize=10%:BoostType=Grad::shrinkage=0.2:MaxDepth=3");

   // Train sample
   factory.TrainAllMethods();

   // Test sample
   factory.TestAllMethods();

   // Evaluate Performance
   factory.EvaluateAllMethods();

   // Save the output
   output->Close();

   cout << "==> Wrote root file: " << output->GetName() << endl;
   cout << "==> TMVAClassification is done!" << endl;

   // TMVA Gui
   //   if (!gROOT->IsBatch()) TMVAGui( filename );
   if (!gROOT->IsBatch())
      TMVA::TMVAGui(filename);

   // For the log file
   cout.rdbuf(oldBuffer);
   logFile.close();
}
