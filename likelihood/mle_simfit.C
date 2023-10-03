// ROOT headers
#include <TROOT.h>
#include <TRandom3.h>
#include <TGraph.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TH1.h>

// cpp headers
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <random>
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

Double_t calc_neglnL(Double_t n_obs, Double_t n_signal, Double_t mu_signal, Double_t n_WZ, Double_t mu_WZ, Double_t n_top, Double_t mu_top, Double_t n_Zjets0, Double_t mu_Zjets0, Double_t n_Zjets1, Double_t mu_Zjets1, Double_t n_Zjets2, Double_t mu_Zjets2, Double_t n_WW, Double_t mu_WW, Double_t n_other)
{
  Double_t lnLi;
  Double_t n_exp;
  n_exp = mu_signal * n_signal + mu_top * n_top + mu_WZ * n_WZ + mu_Zjets2 * n_Zjets2 + mu_Zjets1 * n_Zjets1 + mu_Zjets0 * n_Zjets0 + mu_WW * n_WW + n_other;

  lnLi = +n_obs * log(n_exp) - n_exp;

  return -lnLi;
}

Double_t calc_total_neglnL(Double_t mu_signal, Double_t mu_WZ, Double_t mu_top, Double_t mu_Zjets0, Double_t mu_Zjets1, Double_t mu_Zjets2, Double_t mu_WW)
{
  Double_t lnL_SR;
  Double_t lnL_3l;
  Double_t lnL_Zjets0;
  Double_t lnL_Zjets1;
  Double_t lnL_Zjets2;
  Double_t lnL_emuA;
  Double_t lnL_emuB;
  Double_t lnL_WZ;
  Double_t lnL_total;

  //                   obs   sig   wz   t   z0  z1  z2  ww  ot  sf
  // lnL_SR = calc_neglnL(2891, 1577, mu_signal, 669, mu_WZ, 91, mu_top, 102, mu_Zjets0, 32, mu_Zjets1, 1, mu_Zjets2, 37, mu_WW, 45); //DATA UNBLINDED

  // lnL_SR = calc_neglnL(2636, 1577, mu_signal, 669, mu_WZ, 91, mu_top, 102, mu_Zjets0, 32, mu_Zjets1, 1, mu_Zjets2, 37, mu_WW, 45); //DATA = MC TOTAL

  lnL_3l = calc_neglnL(2409, 2.3, mu_signal, 2099, mu_WZ, 57, mu_top, 35.7, mu_Zjets0, 34.5, mu_Zjets1, 15.7, mu_Zjets2, 0.6, mu_WW, 114.1);

  lnL_emuA = calc_neglnL(2903, 0.4, mu_signal, 22.3, mu_WZ, 1834.8, mu_top, 0.6, mu_Zjets0, 2.2, mu_Zjets1, 0.4, mu_Zjets2, 632.1, mu_WW, 96.9);

  lnL_emuB = calc_neglnL(5736, 0, mu_signal, 0.9, mu_WZ, 5653.7, mu_top, 0, mu_Zjets0, 0.1, mu_Zjets1, 0.1, mu_Zjets2, 13, mu_WW, 4.6);

  lnL_Zjets0 = calc_neglnL(2691, 97.7, mu_signal, 79.4, mu_WZ, 91, mu_top, 1730.1, mu_Zjets0, 0, mu_Zjets1, 0, mu_Zjets2, 47.7, mu_WW, 12.4);

  lnL_Zjets1 = calc_neglnL(5733, 145, mu_signal, 248.2, mu_WZ, 196.5, mu_top, 0, mu_Zjets0, 3798.5, mu_Zjets1, 0, mu_Zjets2, 46, mu_WW, 46.1);

  lnL_Zjets2 = calc_neglnL(5769, 161.9, mu_signal, 308.7, mu_WZ, 244.1, mu_top, 0, mu_Zjets0, 0, mu_Zjets1, 4659.2, mu_Zjets2, 18.4, mu_WW, 57.9);

  lnL_SR = calc_neglnL(2636, 1577, mu_signal, 669.2, mu_WZ, 90.9, mu_top, 101.9, mu_Zjets0, 31.6, mu_Zjets1, 1.2, mu_Zjets2, 37.3, mu_WW, 44.9); // DATA = MC TOTAL

  lnL_total = lnL_SR + lnL_3l + lnL_emuA + lnL_emuB + lnL_Zjets0 + lnL_Zjets1 + lnL_Zjets2;
  // lnL_total = lnL_3l + lnL_emuA + lnL_emuB + lnL_Zjets0 + lnL_Zjets1 + lnL_Zjets2;

  return lnL_total;
}

void mle_simfit()
{

  gROOT->SetBatch(kTRUE);

  auto start = std::chrono::high_resolution_clock::now();

  // Output log file
  ofstream logFile("./sim_fit_new.txt");

  DualStreamBuffer dualBuffer(cout.rdbuf(), logFile.rdbuf());

  std::streambuf *oldBuffer = cout.rdbuf(&dualBuffer);

  Double_t lnL_SR;
  Double_t lnL_3l;
  Double_t lnL_Zjets0;
  Double_t lnL_Zjets1;
  Double_t lnL_Zjets2;
  Double_t lnL_emuA;
  Double_t lnL_emuB;
  Double_t lnL_WZ;
  Double_t lnL_total;
  vector<Double_t> lnL_total_vec;

  Double_t lnL_min;
  int mod = 10000000;
  int min_i;

  // Scaling Factors

  Double_t mu_signal = 1;
  Double_t mu_top = 1;
  Double_t mu_WZ = 1;
  Double_t mu_Zjets0 = 1;
  Double_t mu_Zjets1 = 1;
  Double_t mu_Zjets2 = 1;
  Double_t mu_WW = 1;
  vector<Double_t> sf_vec;

  Double_t min_mu_signal = 1;
  Double_t min_mu_top = 1;
  Double_t min_mu_WZ = 1;
  Double_t min_mu_Zjets0 = 1;
  Double_t min_mu_Zjets1 = 1;
  Double_t min_mu_Zjets2 = 1;
  Double_t min_mu_WW = 1;

  //----------------------------------

  // // Min
  // Double_t low_signal = 1.168 - value;
  // Double_t low_WZ = 1.023 - value;
  // Double_t low_WW = 1.466 - value;
  // Double_t low_top = 1.011 - value;
  // Double_t low_Zjets0 = 1.364 - value;
  // Double_t low_Zjets1 = 1.328 - value;
  // Double_t low_Zjets2 = 1.066 - value;

  // // Max
  // Double_t up_signal = 1.168 + value;
  // Double_t up_WZ = 1.023 + value;
  // Double_t up_WW = 1.466 + value;
  // Double_t up_top = 1.011 + value;
  // Double_t up_Zjets0 = 1.364 + value;
  // Double_t up_Zjets1 = 1.328 + value;
  // Double_t up_Zjets2 = 1.066 + value;

  Double_t value = 0.01;

  // Min
  Double_t low_signal = 1.01 - value;
  Double_t low_WZ = 1.013 - value;
  Double_t low_WW = 1.461 - value;
  Double_t low_top = 1.011 - value;
  Double_t low_Zjets0 = 1.354 - value;
  Double_t low_Zjets1 = 1.323 - value;
  Double_t low_Zjets2 = 1.061 - value;

  // Max
  Double_t up_signal = 1.01 + value;
  Double_t up_WZ = 1.013 + value;
  Double_t up_WW = 1.461 + value;
  Double_t up_top = 1.011 + value;
  Double_t up_Zjets0 = 1.354 + value;
  Double_t up_Zjets1 = 1.323 + value;
  Double_t up_Zjets2 = 1.061 + value;

  Double_t step = 0.1;

  int only_CRs = 0;
  long int iter = 0;
  long int iter_best;
  long int iter_er_plus;
  long int iter_er_minus;

  for (mu_signal = low_signal; mu_signal <= up_signal; mu_signal += step)
  {
    for (mu_WW = low_WW; mu_WW <= up_WW; mu_WW += step)
    {
      for (mu_Zjets0 = low_Zjets0; mu_Zjets0 <= up_Zjets0; mu_Zjets0 += step)
      {
        for (mu_Zjets1 = low_Zjets1; mu_Zjets1 <= up_Zjets1; mu_Zjets1 += step)
        {
          for (mu_Zjets2 = low_Zjets2; mu_Zjets2 <= up_Zjets2; mu_Zjets2 += step)
          {
            for (mu_WZ = low_WZ; mu_WZ <= up_WZ; mu_WZ += step)
            {
              for (mu_top = low_top; mu_top <= up_top; mu_top += step)
              {

                lnL_total = calc_total_neglnL(mu_signal, mu_WZ, mu_top, mu_Zjets0, mu_Zjets1, mu_Zjets2, mu_WW);

                if (iter % mod == 0)
                {
                  cout << " " << iter / mod << "0 mil:  "
                       << " mu_signal = " << mu_signal << "|  mu_WZ = " << mu_WZ << "|  mu_WW = " << mu_WW << "|  mu_top = " << mu_top << "|  mu_Zjets0 = " << mu_Zjets0 << "|  mu_Zjets1 =  " << mu_Zjets1 << "|  mu_Zjets2 = " << mu_Zjets2 << "  -> lnL = " << lnL_total << endl
                       << endl;
                }

                if (lnL_total <= lnL_min or iter == 0)
                {
                  if (lnL_total == lnL_min)
                  {
                    cout << "    SAME VALUE!!   " << endl << endl;
                  }

                  lnL_min = lnL_total;
                  min_mu_signal = mu_signal;
                  min_mu_WZ = mu_WZ;
                  min_mu_WW = mu_WW;
                  min_mu_top = mu_top;
                  min_mu_Zjets0 = mu_Zjets0;
                  min_mu_Zjets1 = mu_Zjets1;
                  min_mu_Zjets2 = mu_Zjets2;
                }

                iter++;
              }
            }
          }
        }
      }
    }
  }

  if (only_CRs == 0)
  {
    cout << "   ----------------------------------------------------------------------------" << endl << endl;
    cout << "   SIMULTANEOUS FIT:        " << endl;
    cout << "   ____________________________        " << endl << endl;
    cout << "   lnL =          " << setprecision(10) << lnL_min << endl << endl;
    cout << "   mu_signal =    " << min_mu_signal << endl << endl;
    cout << "   mu_WZ =        " << min_mu_WZ << endl << endl;
    cout << "   mu_WW =        " << min_mu_WW << endl << endl;
    cout << "   mu_top =       " << min_mu_top << endl << endl;
    cout << "   mu_Zjets0 =    " << min_mu_Zjets0 << endl << endl;
    cout << "   mu_Zjets1 =    " << min_mu_Zjets1 << endl << endl;
    cout << "   mu_Zjets2 =    " << min_mu_Zjets2 << endl << endl;
    cout << "   ----------------------------------------------------------------------------" << endl << endl;
  }
  else if (only_CRs == 1)
  {
    cout << "   ----------------------------------------------------------------------------" << endl << endl;
    cout << "   SIMULTANEOUS FIT - CRs ONLY:        " << endl;
    cout << "   ____________________________        " << endl << endl;
    cout << "   lnL =          " << lnL_min << setprecision(10) << endl << endl;
    cout << "   mu_signal =    " << min_mu_signal << endl << endl;
    cout << "   mu_WZ =        " << min_mu_WZ << endl << endl;
    cout << "   mu_WW =        " << min_mu_WW << endl << endl;
    cout << "   mu_top =       " << min_mu_top << endl << endl;
    cout << "   mu_Zjets0 =    " << min_mu_Zjets0 << endl << endl;
    cout << "   mu_Zjets1 =    " << min_mu_Zjets1 << endl << endl;
    cout << "   mu_Zjets2 =    " << min_mu_Zjets2 << endl << endl;
    cout << "   ----------------------------------------------------------------------------" << endl << endl;
  }

  //Statistical uncertainty

  float sf_low = 1 - 1.;
  float sf_up  = 1 + 1.;
  float sf;
  float sf_min;
  step = 0.0001;
  iter = 0;

  for (sf = sf_low; sf <= sf_up; sf += step)
  {

    //mu_signal, mu_WZ, mu_top, mu_Zjets0, mu_Zjets1, mu_Zjets2, mu_WW

    // lnL_total = calc_total_neglnL(1, 1.011, 1.010, sf, 1.314, 1.050, 1.466);         // SIM CR ONLY
    // lnL_total = calc_total_neglnL(1.005, 1.011, 1.010, 1.384, 1.314, 1.049, 1.466);     // SIM ALL REGIONS
    lnL_total = calc_total_neglnL(sf, 1.012, 1.010, 1.351, 1.323, 1.065, 1.465);     // SIM ALL REGIONS2
    // lnL_total = calc_total_neglnL(0.989, 1.024, 1.011, 1.401, 1.321, 1.051, sf);        // PER REGION FIT

    lnL_total_vec.push_back(lnL_total);
    sf_vec.push_back(sf);

    if (lnL_total <= lnL_min or iter == 0)
    {
      if (lnL_total == lnL_min && iter != 0)
      {
        cout << "    SAME VALUE:  "  << sf << "  :    " << lnL_total << "   |    " << lnL_min <<  endl << endl;
      }

      lnL_min = lnL_total;
      sf_min = sf;
      iter_best = iter;
    }

    iter++;
  }

  TCanvas *c1 = new TCanvas("c1", "-lnL vs mu", 800, 600);
  TCanvas *c2 = new TCanvas("c2", "mu distribution", 800, 600);
  TCanvas *c3 = new TCanvas("c3", "rnd values", 800, 600);

 

  // for (int i = 0; i < lnL_total_vec.size(); i++)
  // {
  //   gr1->SetPoint(i, sf_vec.at(i), lnL_total_vec.at(i));
  // }

  for (int i = 0; i < lnL_total_vec.size(); i++)
  {
    if (lnL_total_vec.at(i) <= lnL_total_vec.at(iter_best) + 0.5)
    {
      iter_er_minus = i;
      break;
    }

  }

  for (int i = iter_best + 1; i < lnL_total_vec.size(); i++)
  {
    if (lnL_total_vec.at(i) >= lnL_total_vec.at(iter_best) + 0.5)
    {
      iter_er_plus = i;
      break;
    }
  }

  // // cout << "   ----------------------------------------------------------------------------" << endl << endl;
  // // cout << "   SIMULTANEOUS FIT - CRs ONLY:        " << endl;
  // // cout << "   ____________________________        " << endl << endl;
  // // cout << "   lnL =          " << lnL_min << setprecision(10) << endl << endl;
  // // cout << "   mu_signal =    " << mu_signal << endl << endl;
  // // cout << "   mu_WZ =        " << mu_WZ << endl << endl;
  // // cout << "   mu_WW =        " << mu_WW << endl << endl;
  // // cout << "   mu_top =       " << mu_top << endl << endl;
  // // cout << "   mu_Zjets0 =    " << mu_Zjets0 << endl << endl;
  // // cout << "   mu_Zjets1 =    " << mu_Zjets1 << endl << endl;
  // // cout << "   mu_Zjets2 =    " << mu_Zjets2 << endl << endl;
  // // cout << "   ----------------------------------------------------------------------------" << endl << endl;

  cout << "   Iter minus:  " << iter_best << endl << endl;
  cout << "   Iter minus:  " << iter_er_minus << "  ->  -" << sf_vec.at(iter_best) - sf_vec.at(iter_er_minus) << endl << endl;
  cout << "   Iter plus:   " << iter_er_plus  << "  ->  +" << sf_vec.at(iter_er_plus) - sf_vec.at(iter_best) << endl << endl;

  // gr1->Draw("ALP");
  // c1->SaveAs("./MLE_plot.png");



  // Timer stop
  auto end = chrono::high_resolution_clock::now();

  chrono::duration<float> duration = end - start;

  cout << endl << "   Script executed in: " << int(duration.count() / 60.0) << " minutes"
       << " and " << int((duration.count() / 60.0 - int(duration.count() / 60.0)) * 60) << " s" << endl
       << endl;

  // For the log file
  std::cout.rdbuf(oldBuffer); // Restore the original stream buffer for std::cout

  logFile.close(); // Close the log file

  return;
}
