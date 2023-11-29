// ROOT headers
#include <TROOT.h>
#include <TRandom1.h>
#include <TRandom2.h>
#include <TRandom3.h>
#include <TGraph.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TH1.h>

// cpp headers
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <random>
#include <chrono> //Timer
using namespace std;

Double_t calc_neglnL(Double_t n_obs, Double_t n_signal, Double_t mu_signal, Double_t n_WZ, Double_t mu_WZ, Double_t n_top, Double_t mu_top, Double_t n_Zjets0, Double_t mu_Zjets0, Double_t n_Zjets1, Double_t mu_Zjets1, Double_t n_Zjets2, Double_t mu_Zjets2, Double_t n_WW, Double_t mu_WW, Double_t n_other)
{
  Double_t lnLi;
  Double_t n_exp;
  n_exp = mu_signal * n_signal + mu_top * n_top + mu_WZ * n_WZ + mu_Zjets2 * n_Zjets2 + mu_Zjets1 * n_Zjets1 + mu_Zjets0 * n_Zjets0 + mu_WW * n_WW + n_other;

  lnLi = +n_obs * log(n_exp) - n_exp;

  return -lnLi;
}




void simfit_statistic_new()
{
  // Timer start
  auto start = std::chrono::high_resolution_clock::now();


  gROOT->SetBatch(kTRUE);


  Double_t xbins[19] = {-20, -4, -3.5, -3, -2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 20};
  
  TH1F *hist_mu = new TH1F("hist_mu", "mu hist", 100, 0.8, 1.2);
  TH1F *hist_pull = new TH1F("hist_pull", "pull hist", 100, -6, 6);
  // TH1F *hist_pull = new TH1F("hist_pull", "pull hist", sizeof(xbins) / sizeof(xbins[0]) - 1, xbins );
  TGraph *gr1 = new TGraph();


  Double_t lnL_SR1;
  Double_t lnL_SR2;
  Double_t lnL_SR;

  Double_t lnL_3l1;
  Double_t lnL_3l2;
  Double_t lnL_3l;

  Double_t lnL_Zjets01;
  Double_t lnL_Zjets02;
  Double_t lnL_Zjets0;

  Double_t lnL_Zjets11;
  Double_t lnL_Zjets12;
  Double_t lnL_Zjets1;

  Double_t lnL_Zjets21;
  Double_t lnL_Zjets22;
  Double_t lnL_Zjets2;

  Double_t lnL_emuA1;
  Double_t lnL_emuA2;
  Double_t lnL_emuA;

  Double_t lnL_emuB1;
  Double_t lnL_emuB2;
  Double_t lnL_emuB;

  Double_t lnL_WZ1;
  Double_t lnL_WZ2;
  Double_t lnL_WZ;

  Double_t lnL_total;


  Double_t lnL_min;
  int mod = 10;

  // Scaling Factors

  Double_t mu_signal = 1.006;
  Double_t mu_top = 1.013;
  Double_t mu_WZ = 1.012;
  Double_t mu_Zjets0 = 1.352;
  Double_t mu_Zjets1 = 1.322;
  Double_t mu_Zjets2 = 1.065;
  Double_t mu_WW = 1.435;

  Double_t min_mu_signal = 1;
  Double_t min_mu_top = 1;
  Double_t min_mu_WZ = 1;
  Double_t min_mu_Zjets0 = 1;
  Double_t min_mu_Zjets1 = 1;
  Double_t min_mu_Zjets2 = 1;
  Double_t min_mu_WW = 1;

  //----------------------------------



  long int iter = 0;
  long int ibest;
  long int iter_er_plus;
  long int iter_er_minus;




  // Data Statistic uncertainty 

  double rnd_values1[7];
  double rnd_values2[7];

  //DATA
  //                       emub  emua    3l   zj2   zj1   zj0    SR  
  double data_values1[7] = {2867, 1445, 1229, 2871, 2836, 1355, 1314};   //DATA1
  double data_values2[7] = {2869, 1458, 1180, 2898, 2897, 1336, 1320};   //DATA2
          
  //                                emub  emua    3l   zj2   zj1   zj0    SR     
  double data_uncertainties1[7] = { 53.5, 38.0, 35.1,  53.6, 53.3, 36.8, sqrt(1314)}; // DATA1
  double data_uncertainties2[7] = { 53.6, 38.2, 34.4 , 53.8, 53.8, 36.6, sqrt(1320)}; // DATA2


  // mu_signal = 1.006;
  // mu_WZ = 1.012;
  // mu_top = 1.013;
  // mu_Zjets0 = 1.352;
  // mu_Zjets1 = 1.323;
  // mu_Zjets2 = 1.065;
  // mu_WW = 1.435;
  
  
  vector<Double_t> min_mu_signal_vec;
  vector<Double_t> pull_param_vec;
  vector<Double_t> error_minus_vec;
  vector<Double_t> error_plus_vec;

  TRandom3 rnd_value;

  Double_t mu_sum = 0;
  Double_t mu_mean = 0;
  Double_t mu_var = 0;
  Double_t mu_var_sum = 0;
  Double_t pull_param = 0;
  Double_t pull_sum = 0;
  Double_t pull_mean = 0;
  Double_t pull_var_sum = 0;
  Double_t pull_var = 0;


 Double_t value = 0.02;

  // Min
  Double_t low_signal = 1.006 - value;
  Double_t low_WZ = 1.012 - value;
  Double_t low_WW = 1.435 - value;
  Double_t low_top = 1.013 - value;
  Double_t low_Zjets0 = 1.352 - value;
  Double_t low_Zjets1 = 1.322 - value;
  Double_t low_Zjets2 = 1.065 - value;

  // Max
  Double_t up_signal = 1.006 + value;
  Double_t up_WZ = 1.012 + value;
  Double_t up_WW = 1.435 + value;
  Double_t up_top = 1.013 + value;
  Double_t up_Zjets0 = 1.352 + value;
  Double_t up_Zjets1 = 1.322 + value;
  Double_t up_Zjets2 = 1.065 + value;
  
  
  int pseudoexp = 1000;
  float step = 0.01;
  float scan = 0.001;

  float error_plus_max = 0;
  float error_minus_max = 0;

  TCanvas *c1 = new TCanvas("c1", "mu distribution", 800, 600);
  TCanvas *c2 = new TCanvas("c2", "pull parameter", 800, 600);


  for (int k = 0; k < pseudoexp; k++)
  {
    vector<Double_t> mu_signal_vec;
    vector<Double_t> lnL_total_vec;

    iter = 0;

    // if (k % mod == 0)
    // { 
      cout << endl;
      cout << "   Pseudoexperiment " << k+1 << endl;
      cout << "   ---------------- " << endl << endl;
    // }

    for (int i = 0; i < 7; i++)
    {
      // uniform_real_distribution<double> gaus_dist(central_values[i], uncertainties[i]);
      // rnd_values[i] = (int) gaus_dist(central_values[i], uncertainties[i]);
      rnd_values1[i] = (int) rnd_value.Gaus(data_values1[i], data_uncertainties1[i]);
      rnd_values2[i] = (int) rnd_value.Gaus(data_values2[i], data_uncertainties2[i]);
    }


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

                  //                       obs   sig   wz   t   z0  z1  z2  ww  ot  sf
                  lnL_SR = calc_neglnL(2891, 1577, mu_signal, 669, mu_WZ, 91, mu_top, 102, mu_Zjets0, 32, mu_Zjets1, 1, mu_Zjets2, 37, mu_WW, 45); //DATA UNBLINDED


                  lnL_3l1 = calc_neglnL(rnd_values1[0], 0.86, mu_signal, 1080.97, mu_WZ, 34.0, mu_top, 2.45, mu_Zjets0, 20.52, mu_Zjets1, 9.18, mu_Zjets2, 0.37, mu_WW, 51.30);
                  lnL_3l2 = calc_neglnL(rnd_values2[0], 1.42, mu_signal, 1018.03, mu_WZ, 22.97, mu_top, 33.28, mu_Zjets0, 13.96, mu_Zjets1, 6.53, mu_Zjets2, 0.27, mu_WW, 62.80);
                  lnL_3l = lnL_3l1 + lnL_3l2;

                  lnL_emuA1 = calc_neglnL(rnd_values1[1], -0.8, mu_signal, 8.99, mu_WZ, 1067.54, mu_top, 0.28, mu_Zjets0, 2.29, mu_Zjets1, 0.19, mu_Zjets2, 206.19, mu_WW, 28.21);
                  lnL_emuA2 = calc_neglnL(rnd_values2[1], 0.49, mu_signal, 13.30, mu_WZ, 767.31, mu_top, 0.35, mu_Zjets0, 0, mu_Zjets1, 0.21, mu_Zjets2, 425.87, mu_WW, 68.67);
                  lnL_emuA = lnL_emuA1 + lnL_emuA2;

                  lnL_emuB1 = calc_neglnL(rnd_values1[2], 0, mu_signal, 0.38, mu_WZ, 2860.72, mu_top, 0, mu_Zjets0, 0., mu_Zjets1, 0.06, mu_Zjets2, 6.07, mu_WW, 0.92);
                  lnL_emuB2 = calc_neglnL(rnd_values2[2], 0, mu_signal, 0.48, mu_WZ, 2792.99, mu_top, 0, mu_Zjets0, 0.12, mu_Zjets1, 0.1, mu_Zjets2, 6.96, mu_WW, 0.68);
                  lnL_emuB = lnL_emuB1 + lnL_emuB2;

                  lnL_Zjets01 = calc_neglnL(rnd_values1[3], 36.96, mu_signal, 40.47, mu_WZ, 63.07, mu_top, 855.74, mu_Zjets0, 0, mu_Zjets1, 0, mu_Zjets2, 19.71, mu_WW, 7.57);
                  lnL_Zjets02 = calc_neglnL(rnd_values2[3], 60.77, mu_signal, 38.90, mu_WZ, 27.94, mu_top, 874.37, mu_Zjets0, 0, mu_Zjets1, 0, mu_Zjets2, 27.96, mu_WW, 4.84);
                  lnL_Zjets0 = lnL_Zjets01 + lnL_Zjets02;

                  lnL_Zjets11 = calc_neglnL(rnd_values1[4], 86.51, mu_signal, 144.37, mu_WZ, 143.99, mu_top, 0, mu_Zjets0, 1790.57, mu_Zjets1, 0, mu_Zjets2, 34.78, mu_WW, 18.44);
                  lnL_Zjets12 = calc_neglnL(rnd_values2[4], 58.53, mu_signal, 103.87, mu_WZ, 52.53, mu_top, 0, mu_Zjets0, 2007.9, mu_Zjets1, 0, mu_Zjets2, 11.20, mu_WW, 27.66);
                  lnL_Zjets1 = lnL_Zjets11 + lnL_Zjets12;

                  lnL_Zjets21 = calc_neglnL(rnd_values1[5], 83.64, mu_signal, 158.64, mu_WZ, 152.34, mu_top, 0, mu_Zjets0, 0, mu_Zjets1, 2296.7, mu_Zjets2, 12.44, mu_WW, 24.86);
                  lnL_Zjets22 = calc_neglnL(rnd_values2[5], 78.29, mu_signal, 150.05, mu_WZ, 91.75, mu_top, 0, mu_Zjets0, 0, mu_Zjets1, 2632.56, mu_Zjets2, 5.95, mu_WW, 33.03);
                  lnL_Zjets2 = lnL_Zjets21 + lnL_Zjets22;

                  lnL_SR1 = calc_neglnL(rnd_values1[6], 713.47, mu_signal, 356.24, mu_WZ, 71.79, mu_top, 51.78, mu_Zjets0, 31.18, mu_Zjets1, 0.94, mu_Zjets2, 19.70, mu_WW, 22.82); // DATA = MC TOTAL
                  lnL_SR2 = calc_neglnL(rnd_values2[6], 863.56, mu_signal, 312.95, mu_WZ, 19.14, mu_top, 50.12, mu_Zjets0, 0.47, mu_Zjets1, 0.22, mu_Zjets2, 17.66, mu_WW, 22.12); // DATA = MC TOTAL
                  lnL_SR = lnL_SR1 + lnL_SR2;

                  lnL_total = lnL_SR + lnL_3l + lnL_emuA + lnL_emuB + lnL_Zjets0 + lnL_Zjets1 + lnL_Zjets2;
                  lnL_total = lnL_3l + lnL_emuA + lnL_emuB + lnL_Zjets0 + lnL_Zjets1 + lnL_Zjets2;

                  if (lnL_total <= lnL_min || iter == 0)
                  {
                    lnL_min = lnL_total;
                    min_mu_signal = mu_signal;
                    min_mu_top = mu_top;
                    min_mu_WZ = mu_WZ;
                    min_mu_WW = mu_WW;
                    min_mu_Zjets0 = mu_Zjets0;
                    min_mu_Zjets1 = mu_Zjets1;
                    min_mu_Zjets2 = mu_Zjets2;
                    if (lnL_total == lnL_min && iter != 0)
                    {
                      cout << "    SAME VALUE:  " << mu_signal << "  :    " << lnL_total << "   |    " << lnL_min << endl << endl;
                    }
                  }
                }
              }
            }
          }
        }
      }
    }


    min_mu_signal = 1.006;
    min_mu_WZ = 1.012;
    min_mu_top = 1.013;
    min_mu_Zjets0 = 1.352;
    min_mu_Zjets1 = 1.322;
    min_mu_Zjets2 = 1.065;
    min_mu_WW = 1.435;

    float value2 = 0.3;

    for (mu_signal = 1 - value2; mu_signal < 1 + value2; mu_signal += scan)
    // for (mu_signal = 1 - value; mu_signal < 1 + value; mu_signal += scan)
    {
      //                       obs   sig   wz   t   z0  z1  z2  ww  ot  sf
      lnL_3l1 = calc_neglnL(rnd_values1[0], 0.86, mu_signal, 1080.97, min_mu_WZ, 34.0, min_mu_top, 2.45, min_mu_Zjets0, 20.52, min_mu_Zjets1, 9.18, min_mu_Zjets2, 0.37, min_mu_WW, 51.30);
      lnL_3l2 = calc_neglnL(rnd_values2[0], 1.42, mu_signal, 1018.03, min_mu_WZ, 22.97, min_mu_top, 33.28, min_mu_Zjets0, 13.96, min_mu_Zjets1, 6.53, min_mu_Zjets2, 0.27, min_mu_WW, 62.80);
      lnL_3l = lnL_3l1 + lnL_3l2;

      lnL_emuA1 = calc_neglnL(rnd_values1[1], -0.8, mu_signal, 8.99, min_mu_WZ, 1067.54, min_mu_top, 0.28, min_mu_Zjets0, 2.29, min_mu_Zjets1, 0.19, min_mu_Zjets2, 206.19, min_mu_WW, 28.21);
      lnL_emuA2 = calc_neglnL(rnd_values2[1], 0.49, mu_signal, 13.30, min_mu_WZ, 767.31, min_mu_top, 0.35, min_mu_Zjets0, 0, min_mu_Zjets1, 0.21, min_mu_Zjets2, 425.87, min_mu_WW, 68.67);
      lnL_emuA = lnL_emuA1 + lnL_emuA2;

      lnL_emuB1 = calc_neglnL(rnd_values1[2], 0, mu_signal, 0.38, min_mu_WZ, 2860.72, min_mu_top, 0, min_mu_Zjets0, 0., min_mu_Zjets1, 0.06, min_mu_Zjets2, 6.07, min_mu_WW, 0.92);
      lnL_emuB2 = calc_neglnL(rnd_values2[2], 0, mu_signal, 0.48, min_mu_WZ, 2792.99, min_mu_top, 0, min_mu_Zjets0, 0.12, min_mu_Zjets1, 0.1, min_mu_Zjets2, 6.96, min_mu_WW, 0.68);
      lnL_emuB = lnL_emuB1 + lnL_emuB2;

      lnL_Zjets01 = calc_neglnL(rnd_values1[3], 36.96, mu_signal, 40.47, min_mu_WZ, 63.07, min_mu_top, 855.74, min_mu_Zjets0, 0, min_mu_Zjets1, 0, min_mu_Zjets2, 19.71, min_mu_WW, 7.57);
      lnL_Zjets02 = calc_neglnL(rnd_values2[3], 60.77, mu_signal, 38.90, min_mu_WZ, 27.94, min_mu_top, 874.37, min_mu_Zjets0, 0, min_mu_Zjets1, 0, min_mu_Zjets2, 27.96, min_mu_WW, 4.84);
      lnL_Zjets0 = lnL_Zjets01 + lnL_Zjets02;

      lnL_Zjets11 = calc_neglnL(rnd_values1[4], 86.51, mu_signal, 144.37, min_mu_WZ, 143.99, min_mu_top, 0, min_mu_Zjets0, 1790.57, min_mu_Zjets1, 0, min_mu_Zjets2, 34.78, min_mu_WW, 18.44);
      lnL_Zjets12 = calc_neglnL(rnd_values2[4], 58.53, mu_signal, 103.87, min_mu_WZ, 52.53, min_mu_top, 0, min_mu_Zjets0, 2007.9, min_mu_Zjets1, 0, min_mu_Zjets2, 11.20, min_mu_WW, 27.66);
      lnL_Zjets1 = lnL_Zjets11 + lnL_Zjets12;

      lnL_Zjets21 = calc_neglnL(rnd_values1[5], 83.64, mu_signal, 158.64, min_mu_WZ, 152.34, min_mu_top, 0, min_mu_Zjets0, 0, min_mu_Zjets1, 2296.7, min_mu_Zjets2, 12.44, min_mu_WW, 24.86);
      lnL_Zjets22 = calc_neglnL(rnd_values2[5], 78.29, mu_signal, 150.05, min_mu_WZ, 91.75, min_mu_top, 0, min_mu_Zjets0, 0, min_mu_Zjets1, 2632.56, min_mu_Zjets2, 5.95, min_mu_WW, 33.03);
      lnL_Zjets2 = lnL_Zjets21 + lnL_Zjets22;

      lnL_SR1 = calc_neglnL(rnd_values1[6], 713.47, mu_signal, 356.24, min_mu_WZ, 71.79, min_mu_top, 51.78, min_mu_Zjets0, 31.18, min_mu_Zjets1, 0.94, min_mu_Zjets2, 19.70, min_mu_WW, 22.82); // DATA = MC TOTAL
      lnL_SR2 = calc_neglnL(rnd_values2[6], 863.56, mu_signal, 312.95, min_mu_WZ, 19.14, min_mu_top, 50.12, min_mu_Zjets0, 0.47, min_mu_Zjets1, 0.22, min_mu_Zjets2, 17.66, min_mu_WW, 22.12); // DATA = MC TOTAL
      lnL_SR = lnL_SR1 + lnL_SR2;

      lnL_total = lnL_SR + lnL_3l + lnL_emuA + lnL_emuB + lnL_Zjets0 + lnL_Zjets1 + lnL_Zjets2;

      if (lnL_total <= lnL_min || iter == 0)
      {
        lnL_min = lnL_total;
        min_mu_signal = mu_signal;
        ibest = iter;
        // if (lnL_total == lnL_min && iter != 0)
        // {
        //   cout << "    SAME VALUE:  " << mu_signal << "  :    " << lnL_total << "   |    " << lnL_min << endl << endl;
        // }
      }

      lnL_total_vec.push_back(lnL_total);
      mu_signal_vec.push_back(mu_signal);

      iter++;
    }

    for (int i = 0; i < lnL_total_vec.size(); i++)
    {
      if (lnL_total_vec.at(i) <= lnL_total_vec.at(ibest) + 0.5)
      {
        iter_er_minus = i;
        break;
      }
    }

    for (int i = ibest + 1; i < lnL_total_vec.size(); i++)
    {
      if (lnL_total_vec.at(i) >= lnL_total_vec.at(ibest) + 0.5)
      {
        iter_er_plus = i;
        break;
      }
    }

    if (min_mu_signal - mu_signal_vec.at(iter_er_minus) > error_minus_max)
    {
      error_minus_max = min_mu_signal - mu_signal_vec.at(iter_er_minus);
    }

    if (mu_signal_vec.at(iter_er_plus) - min_mu_signal > error_plus_max)
    {
      error_plus_max = mu_signal_vec.at(iter_er_plus) - min_mu_signal;
    }

    min_mu_signal_vec.push_back(min_mu_signal);
    error_minus_vec.push_back(min_mu_signal - mu_signal_vec.at(iter_er_minus));
    error_plus_vec.push_back(mu_signal_vec.at(iter_er_plus) - min_mu_signal);


    cout << "   mu_S   =  " << min_mu_signal_vec.at(k) << endl << endl;
    cout << "   error+ = +" << mu_signal_vec.at(iter_er_plus) - min_mu_signal << endl << endl;
    cout << "   error- = -" << min_mu_signal - mu_signal_vec.at(iter_er_minus) << endl << endl;
    

    // for (int i = 0; i < lnL_total_vec.size(); i++)
    // {
    //   gr1->SetPoint(i, mu_signal_vec.at(i), lnL_total_vec.at(i));
    // }

    pull_param = (min_mu_signal - 1.006) / ((error_plus_vec.at(k) + error_minus_vec.at(k)) / 2);
    pull_param_vec.push_back(pull_param);


    hist_mu->Fill(min_mu_signal);
    hist_mu->Draw();



    hist_pull->Fill(pull_param);
    hist_pull->Draw();

    if ((k+1) % mod == 0)
    {
      float mean_mu = hist_mu->GetMean();
      float mean_pull = hist_pull->GetMean();
      float rms_mu = hist_mu->GetRMS();
      float rms_pull = hist_pull->GetRMS();
      TF1 *gauss_mu = new TF1("gaussianFit", "gaus", mean_mu, rms_mu);
      TF1 *gauss_pull = new TF1("gaussianFitPull", "gaus", mean_pull, rms_pull);

      gauss_mu->SetLineColor(kBlack);
      gauss_pull->SetLineColor(kBlack);

      hist_mu->Fit(gauss_mu);
      hist_pull->Fit(gauss_pull);

      TLatex *tex1;
      TLatex *tex2;
      TLatex *tex3;
      TLatex *tex4;
      TLegend *leg;
      string pseudoexpStr = Form("%.2d", pseudoexp);

      tex1 = new TLatex(0.15, 0.8, "#mu_{S} Histogram");
      tex1->SetNDC();
      tex1->SetTextSize(0.04);

      tex2 = new TLatex(0.15, 0.72, "Data Statistics");
      tex2->SetNDC();
      tex2->SetTextSize(0.04);

      tex3 = new TLatex(0.15, 0.64, ("Experiments: " + pseudoexpStr).c_str());
      tex3->SetNDC();
      tex3->SetTextSize(0.04);

      TLatex *tex_mu_mean = new TLatex(0.7, 0.8, Form("mean = %.3f", gauss_mu->GetParameter(1)));
      tex_mu_mean->SetNDC();
      tex_mu_mean->SetTextSize(0.035);

      TLatex *tex_mu_rms = new TLatex(0.7, 0.7, Form("rms = %.3f", gauss_mu->GetParameter(2)));
      tex_mu_rms->SetNDC();
      tex_mu_rms->SetTextSize(0.035);

      TLatex *tex_pull_mean = new TLatex(0.7, 0.8, Form("mean = %.3f", gauss_pull->GetParameter(1)));
      tex_pull_mean->SetNDC();
      tex_pull_mean->SetTextSize(0.035);

      TLatex *tex_pull_rms = new TLatex(0.7, 0.7, Form("rms = %.3f", gauss_pull->GetParameter(2)));
      tex_pull_rms->SetNDC();
      tex_pull_rms->SetTextSize(0.035);

      c1->cd();
      c1->SetTicks();

      hist_mu->SetLineColor(kBlack);
      hist_mu->SetLineWidth(1);
      hist_mu->SetStats(0);
      hist_mu->SetFillColorAlpha(TColor::GetColor("#48c9b0"), 1.);
      hist_mu->SetTitle("  ");
      hist_mu->GetXaxis()->SetTitleSize(0.04);
      hist_mu->GetXaxis()->SetTitleFont(62);
      hist_mu->GetYaxis()->SetTitleFont(62);
      hist_mu->GetXaxis()->SetTitle("#hat{#mu}_{S}");
      hist_mu->GetYaxis()->SetTitle("Entries");
      hist_mu->GetXaxis()->SetTitleOffset(1.1);
      hist_mu->Draw();

      tex1->Draw("same");
      tex2->Draw("same");
      tex3->Draw("same");
      tex_mu_mean->Draw("same");
      tex_mu_rms->Draw("same");

      c1->SaveAs("./errors/stat_mu.png");

      tex1 = new TLatex(0.15, 0.8, "Pull Histogram");
      tex1->SetNDC();
      tex1->SetTextSize(0.04);

      tex2 = new TLatex(0.15, 0.72, "Data Statistics");
      tex2->SetNDC();
      tex2->SetTextSize(0.04);

      tex3 = new TLatex(0.15, 0.64, ("Experiments: " + pseudoexpStr).c_str());
      tex3->SetNDC();
      tex3->SetTextSize(0.04);

      c2->cd();
      c2->SetTicks();
      hist_pull->SetLineColor(kBlack);
      hist_pull->SetLineWidth(1);
      hist_pull->SetStats(0);
      hist_pull->SetFillColorAlpha(TColor::GetColor("#48c9b0"), 1.);
      hist_pull->SetTitle("  ");
      hist_pull->GetXaxis()->SetTitleSize(0.04);
      hist_pull->GetXaxis()->SetTitleFont(62);
      hist_pull->GetYaxis()->SetTitleFont(62);
      hist_pull->GetXaxis()->SetTitle("(#hat{#mu}_{S} - #mu_{S}) / #hat{#sigma}_{S}");
      hist_pull->GetYaxis()->SetTitle("Entries");
      hist_pull->GetXaxis()->SetTitleOffset(1.1);
      hist_pull->Draw();
      // gr1->Draw("ALP");

      tex1->Draw("same");
      tex2->Draw("same");
      tex3->Draw("same");
      tex_pull_mean->Draw("same");
      tex_pull_rms->Draw("same");

      c2->SaveAs("./errors/stat_pull.png");
    }
  }

  for (int i = 0; i < min_mu_signal_vec.size(); i++)
  {
    mu_sum += min_mu_signal_vec.at(i);
    pull_sum += pull_param_vec.at(i);
  }

  mu_mean = mu_sum / min_mu_signal_vec.size();
  pull_mean = pull_sum / pull_param_vec.size();


  for (int i = 0; i < min_mu_signal_vec.size(); i++)
  {
    mu_var_sum += pow(min_mu_signal_vec.at(i) - mu_mean, 2);
    pull_var_sum += pow((min_mu_signal_vec.at(i) - 1.005) / ((error_plus_vec.at(i) + error_minus_vec.at(i))/2) - pull_mean, 2);
  }


  mu_var = mu_var_sum / min_mu_signal_vec.size();
  pull_var = pull_var_sum / min_mu_signal_vec.size();

  cout << endl << endl;
  cout << "   ----------------------------------------------"  << endl << endl;
  cout << "   min lnL = " <<  lnL_min << endl << endl;    
  cout << "   mean(mu) = " <<  mu_mean << endl << endl;    
  cout << "   std(mu) = " <<  sqrt(mu_var) << endl << endl;    
  cout << "   mean(pull) = " <<  pull_mean << endl << endl;    
  cout << "   std(pull) = " <<  sqrt(pull_var) << endl << endl; 
  cout << "   mu_plus_max =     +" << error_plus_max << endl << endl;
  cout << "   mu_minus_max =    -" << error_minus_max << endl << endl;
  cout << "   ----------------------------------------------"  << endl << endl;

  // Timer stop
  auto end = chrono::high_resolution_clock::now();
  chrono::duration<float> duration = end - start;
  cout << endl << "   Script executed in: " << int(duration.count() / 60.0) << " minutes" << " and " <<
  int((duration.count() / 60.0 - int(duration.count() / 60.0))*60) << " s" <<  endl << endl;

  return;
}
