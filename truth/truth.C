
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


vector<Float_t> Counter(string category, string index1, string index2)
{

  TH1::SetDefaultSumw2(kTRUE);

  
  Double_t weight;
  vector<Float_t> events;
  events.clear();
  

  Double_t signal = 0.;
  Double_t signaler = 0.;

  string filepath = "/home/touk/Desktop/touk/master/thesis/data/SR_nocut_aggelos/" + category + index1 + index2;
  TFile *file = new TFile((string(filepath)).c_str());
  TTree *tree = file->Get<TTree>("tree");

  Int_t nentries = (Int_t)tree->GetEntries();
  tree->SetBranchAddress("global_weight", &weight);


  for (Int_t i = 0; i < nentries; i++)
  {
      tree->GetEntry(i);
      signal = signal + weight;
      signaler = signaler + weight*weight;
  }

  cout << "    ENTRIES = " << tree->GetEntries() << endl << endl; 
  cout << "          N = " << signal << "+-" << sqrt(signaler) << endl << endl; 


  events.push_back(signal);
  events.push_back(sqrt(signaler));


  return events;


}



// MAIN
void truth()
{

    Double_t weight;
    Float_t signal_qq = 0.;
    Float_t signal_qq_er = 0.;
    Float_t signal_gg = 0.;
    Float_t signal_gg_er = 0.;
    Float_t signal_EWK = 0.;
    Float_t signal_EWK_er = 0.;

    vector<Float_t> signal_qq1;
    vector<Float_t> signal_qq2;
    vector<Float_t> signal_qq3;
    vector<Float_t> signal_gg1;
    vector<Float_t> signal_gg2;
    vector<Float_t> signal_gg3;
    vector<Float_t> signal_EWK1;
    vector<Float_t> signal_EWK2;
    vector<Float_t> signal_EWK3;

    TH1::SetDefaultSumw2(kTRUE);

    vector<string> categories = {"345666.Sherpa_222_NNPDF30NNLO_", "345723.Sherpa_222_NNPDF30NNLO_", "363724.MadGraphPythia8EvtGen_"};
    vector<string> indexes1 = {"llvvZZ.deriv.DAOD_STDM3.e6240_s3126_", "ggllvvZZ.deriv.DAOD_STDM3.e6213_s3126_", "ZZllvv2jEWK.deriv.DAOD_STDM3.e5712_s3126_"};
    vector<string> indexes2 = {"r9364_p4252.root.root", "r10201_p4252.root.root", "r10724_p4252.root.root"};

    for (string &category : categories)
    {
      for (string &index1 : indexes1)
      {
          for (string &index2 : indexes2)
          {
              if (category == "345666.Sherpa_222_NNPDF30NNLO_" && index1 == "llvvZZ.deriv.DAOD_STDM3.e6240_s3126_")
              {
                  if (index2 == "r9364_p4252.root.root")
                  {
                    cout << "   QCD qq1" << endl << endl;
                    signal_qq1 = Counter(category, index1, index2);
                  }
                  else if (index2 == "r10201_p4252.root.root")
                  {
                    cout << "   QCD qq2" << endl << endl;
                    signal_qq2 = Counter(category, index1, index2);
                  }
                  else if (index2 == "r10724_p4252.root.root")
                  {
                    cout << "   QCD qq3" << endl << endl;
                    signal_qq3 = Counter(category, index1, index2);
                  }
              }
              else if (category == "345723.Sherpa_222_NNPDF30NNLO_" && index1 == "ggllvvZZ.deriv.DAOD_STDM3.e6213_s3126_")
              {
                  if (index2 == "r9364_p4252.root.root")
                  {
                    cout << "   QCD gg1" << endl << endl;
                    signal_gg1 = Counter(category, index1, index2);
                  }
                  else if (index2 == "r10201_p4252.root.root")
                  {
                    cout << "   QCD gg2" << endl << endl;
                    signal_gg2 = Counter(category, index1, index2);
                  }
                  else if (index2 == "r10724_p4252.root.root")
                  {
                    cout << "   QCD gg3" << endl << endl;
                    signal_gg3 = Counter(category, index1, index2);
                  }
              }
              else if (category == "363724.MadGraphPythia8EvtGen_" && index1 == "ZZllvv2jEWK.deriv.DAOD_STDM3.e5712_s3126_")
              {
                  if (index2 == "r9364_p4252.root.root")
                  {
                    cout << "   QCD EWK2" << endl << endl;
                    signal_EWK1 = Counter(category, index1, index2);
                  }
                  else if (index2 == "r10201_p4252.root.root")
                  {
                    cout << "   QCD EWK2" << endl << endl;
                    signal_EWK2 = Counter(category, index1, index2);
                  }
                  else if (index2 == "r10724_p4252.root.root")
                  {
                    cout << "   QCD EWK3" << endl << endl;
                    signal_EWK3 = Counter(category, index1, index2);
                  }
              }
          }
      }
    }

    signal_qq = signal_qq1[0] + signal_qq2[0] +  signal_qq3[0];
    signal_qq_er = signal_qq1[1] + signal_qq2[1] +  signal_qq3[1];
    signal_gg = signal_gg1[0] + signal_gg2[0] +  signal_gg3[0];
    signal_gg_er = signal_gg1[1] + signal_gg2[1] +  signal_gg3[1];
    signal_EWK = signal_EWK1[0] + signal_EWK2[0] +  signal_EWK3[0];
    signal_EWK_er = signal_EWK1[1] + signal_EWK2[1] +  signal_EWK3[1];

    cout << "Signal QCD qq:  " << signal_qq << " +- " << signal_qq_er << endl << endl;
    cout << "Signal QCD gg:  " << signal_gg << " +- " << signal_gg_er << endl << endl;
    cout << "Signal EWK:     " << signal_EWK << " +- " << signal_EWK_er << endl << endl;
}
