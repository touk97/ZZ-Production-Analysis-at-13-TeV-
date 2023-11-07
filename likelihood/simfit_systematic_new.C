// ROOT headers
#include <TROOT.h>
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


void simfit_systematic_new()
{
  // Timer start
  auto start = std::chrono::high_resolution_clock::now();


  gROOT->SetBatch(kTRUE);

  TH1F *hist_mu = new TH1F("hist_mu", "mu hist", 80, 0.8, 1.2);
  TH1F *hist_pull = new TH1F("hist_pull", "pull hist", 80, -6, 6);
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
  int mod = 100;

  // Scaling Factors

  Double_t mu_signal = 1.006;
  Double_t mu_top = 1.013;
  Double_t mu_WZ = 1.012;
  Double_t mu_Zjets0 = 1.352;
  Double_t mu_Zjets1 = 1.322;
  Double_t mu_Zjets2 = 1.065;
  Double_t mu_WW = 1.435;

  Double_t min_mu_signal = 1.006;
  Double_t min_mu_WZ = 1.012;
  Double_t min_mu_top = 1.013;
  Double_t min_mu_Zjets0 = 1.352;
  Double_t min_mu_Zjets1 = 1.323;
  Double_t min_mu_Zjets2 = 1.065;
  Double_t min_mu_WW = 1.435;

  //----------------------------------



  long int iter = 0;
  long int ibest;
  long int iter_er_plus;
  long int iter_er_minus;




  // MC Systematic Uncertainty


  double rnd_values[14][8];

  double central_values[14][8] = {
      // sig, WZ, top, zj0, zj1, zj2, WW, othr
      {0.86, 1080.97, 34.0, 2.45, 20.52, 9.18, 0.37, 51.30},   // 3l1
      {1.42, 1018.03, 22.97, 33.28, 13.96, 6.53, 0.27, 62.80}, // 3l2

      {0., 8.99, 1067.54, 0.28, 2.29, 0.19, 206.19, 28.21}, // emuA1
      {0.49, 13.30, 767.31, 0.35, 0, 0.21, 425.87, 68.67}, // emuA2

      {0, 0.38, 2860.72, 0, 0, 0.06, 6.07, 0.92},   // emuB1
      {0, 0.48, 2792.99, 0, 0.12, 0.1, 6.96, 0.68}, // emuB2

      {36.96, 40.47, 63.07, 855.74, 0, 0, 19.71, 7.57}, // Zjets01
      {60.77, 38.90, 27.94, 874.37, 0, 0, 27.96, 4.84}, // Zjets02

      {86.51, 144.37, 143.99, 0, 1790.57, 0, 34.78, 18.44}, // Zjets11
      {58.53, 103.87, 52.53, 0, 2007.9, 0, 11.20, 27.66},   // Zjets12

      {83.64, 158.64, 152.34, 0, 0, 2296.7, 12.44, 24.86}, // Zjets21
      {78.29, 150.05, 91.75, 0, 0, 2362.56, 5.95, 33.03},  // Zjets22

      {713.47, 356.24, 71.79, 51.78, 31.18, 0.94, 19.70, 22.82}, // SR11
      {863.56, 312.95, 19.14, 50.12, 0.47, 0.22, 17.66, 22.12}  // SR12
  };

  double uncertainties[14][8] = {
      // sig, WZ, top, zj0, zj1, zj2, WW, othr

      {0.23, 7.11, 1.35, 3.59, 5.14, 3.28, 0.10, 0.59}, // 3l1
      {0.49, 7.36, 1.22, 7.47, 2.72, 2.44, 0.09, 0.72}, // 3l2

      {0.08, 0.59, 8.25, 0.29, 1.27, 0.31, 2.62, 8.09},  // emuA1
      {0.18, 0.83, 7.48, 0.32, 0.40, 0.18, 3.71, 14.27}, // emuA2

      {0, 0.11, 12.46, 0, 0, 0.03, 0.45, 0.30},    // emuB1
      {0, 0.13, 12.45, 0, 0.05, 0.04, 0.54, 1.95}, // emuB1

      {2.49, 1.44, 2.12, 56.10, 0, 0, 0.78, 2.41}, // Zjets01
      {2.58, 1.40, 1.63, 82.34, 0, 0, 0.94, 0.32}, // Zjets02

      {3.08, 2.61, 3.14, 0, 64.13, 0, 1.03, 1.47}, // Zjets11
      {2.55, 2.29, 1.88, 0, 62.17, 0, 0.60, 2.04}, // Zjets12
      
      {1.64, 2.07, 2.95, 0, 0, 79.29, 0.64, 1.56}, // Zjets21
      {1.86, 2.21, 2.35, 0, 0, 77.61, 0.46, 2.04}, // Zjets22
      
      {9.48, 4.45, 2.21, 21.54, 8.37, 0.7, 0.78, 1.45},   // SR1
      {10.07, 4.23, 1.18, 10.51, 2.35, 0.16, 0.74, 1.47}, // SR2
  };

  // double central_values[7][8] = {
  //     {2.285, 2098.99, 56.962, 35.737, 34.480, 15.704, 0.646, 114.097},              // 3l
  //     {0.414, 22.295, 1834.840, 0.631, 2.222, 0.398, 632.067, 96.879},            // emuA
  //     {0.004, 0.869, 5653.710, 0, 0.119, 0.153, 13.033, 4.564},            // emuB
  //     {97.732, 79.368, 91.010, 1730.100, 0, 0, 47.677, 12.409},                      // Zjets0
  //     {145.039, 248.238, 196.525, 0, 3798.460, 0, 45.983, 46.104},                // Zjets1
  //     {161.935, 308.685, 244.092, 0, 0, 4659.240, 18.390, 57.887},            // Zjets2
  //     {1577.03, 669.197, 90.937, 101.899, 31.651, 1.161, 37.355, 44.945}           // SR
  // };

  // double uncertainties[7][8] = {
  //     {0.542, 10.232, 1.818, 8.284, 5.812, 4.090, 0.138, 0.930},  // 3l
  //     {0.202, 1.018, 11.134, 0.431, 1.332, 0.362, 4.548, 16.403}, // emuA
  //     {0.004, 0.167, 17.611, 0, 0.119, 0.153, 0.706, 1.977},      // emuB
  //     {3.587, 2.012, 2.676, 99.639, 0, 0, 1.218, 2.427},          // Zjets0
  //     {3.999, 3.471, 3.660, 0, 89.313, 0, 1.191, 2.517},          // Zjets1
  //     {2.478, 3.025, 3.771, 0, 0, 110.953, 0.788, 2.567},         // Zjets2
  //     {13.836, 6.144, 2.503, 23.970, 8.693, 0.718, 1.078, 2.060}, // SR
  // };

  
  
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

  Double_t value = 0.002;

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
  
  int pseudoexp = 2000;
  float step = 0.001;
  float scan = 0.0001;

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

    for (int i = 0; i < 14; i++)
    {
      for (int j = 0; j < 8; j++)
      {
        rnd_values[i][j] = rnd_value.Gaus(central_values[i][j], uncertainties[i][j]);
      }
    }

    iter = 0;

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

                  //                   obs   sig   wz   t   z0  z1  z2  ww  ot  sf
                  // lnL_SR = calc_neglnL(2891, 1577, mu_signal, 669, mu_WZ, 91, mu_top, 102, mu_Zjets0, 32, mu_Zjets1, 1, mu_Zjets2, 37, mu_WW, 45); //DATA UNBLINDED

                  lnL_3l1 = calc_neglnL(1229, rnd_values[0][0], mu_signal, rnd_values[0][1], mu_WZ, rnd_values[0][2], mu_top, rnd_values[0][3], mu_Zjets0, rnd_values[0][4], mu_Zjets1, rnd_values[0][5], mu_Zjets2, rnd_values[0][6], mu_WW, rnd_values[0][7]);
                  lnL_3l2 = calc_neglnL(1180, rnd_values[1][0], mu_signal, rnd_values[1][1], mu_WZ, rnd_values[1][2], mu_top, rnd_values[1][3], mu_Zjets0, rnd_values[1][4], mu_Zjets1, rnd_values[1][5], mu_Zjets2, rnd_values[1][6], mu_WW, rnd_values[1][7]);
                  lnL_3l = lnL_3l1 + lnL_3l2;

                  lnL_emuA1 = calc_neglnL(1445, rnd_values[2][0], mu_signal, rnd_values[2][1], mu_WZ, rnd_values[2][2], mu_top, rnd_values[2][3], mu_Zjets0, rnd_values[2][4], mu_Zjets1, rnd_values[2][5], mu_Zjets2, rnd_values[2][6], mu_WW, rnd_values[2][7]);
                  lnL_emuA2 = calc_neglnL(1458, rnd_values[3][0], mu_signal, rnd_values[3][1], mu_WZ, rnd_values[3][2], mu_top, rnd_values[3][3], mu_Zjets0, rnd_values[3][4], mu_Zjets1, rnd_values[3][5], mu_Zjets2, rnd_values[3][6], mu_WW, rnd_values[3][7]);
                  lnL_emuA = lnL_emuA1 + lnL_emuA2;

                  lnL_emuB1 = calc_neglnL(2867, rnd_values[4][0], mu_signal, rnd_values[4][1], mu_WZ, rnd_values[4][2], mu_top, rnd_values[4][3], mu_Zjets0, rnd_values[4][4], mu_Zjets1, rnd_values[4][5], mu_Zjets2, rnd_values[4][6], mu_WW, rnd_values[4][7]);
                  lnL_emuB2 = calc_neglnL(2869, rnd_values[5][0], mu_signal, rnd_values[5][1], mu_WZ, rnd_values[5][2], mu_top, rnd_values[5][3], mu_Zjets0, rnd_values[5][4], mu_Zjets1, rnd_values[5][5], mu_Zjets2, rnd_values[5][6], mu_WW, rnd_values[5][7]);
                  lnL_emuB = lnL_emuB1 + lnL_emuB2;

                  lnL_Zjets01 = calc_neglnL(1355, rnd_values[6][0], mu_signal, rnd_values[6][1], mu_WZ, rnd_values[6][2], mu_top, rnd_values[6][3], mu_Zjets0, rnd_values[6][4], mu_Zjets1, rnd_values[6][5], mu_Zjets2, rnd_values[6][6], mu_WW, rnd_values[6][7]);
                  lnL_Zjets02 = calc_neglnL(1336, rnd_values[7][0], mu_signal, rnd_values[7][1], mu_WZ, rnd_values[7][2], mu_top, rnd_values[7][3], mu_Zjets0, rnd_values[7][4], mu_Zjets1, rnd_values[7][5], mu_Zjets2, rnd_values[7][6], mu_WW, rnd_values[7][7]);
                  lnL_Zjets0 = lnL_Zjets01 + lnL_Zjets02;

                  lnL_Zjets11 = calc_neglnL(2836, rnd_values[8][0], mu_signal, rnd_values[8][1], mu_WZ, rnd_values[8][2], mu_top, rnd_values[8][3], mu_Zjets0, rnd_values[8][4], mu_Zjets1, rnd_values[8][5], mu_Zjets2, rnd_values[8][6], mu_WW, rnd_values[8][7]);
                  lnL_Zjets12 = calc_neglnL(2897, rnd_values[9][0], mu_signal, rnd_values[9][1], mu_WZ, rnd_values[9][2], mu_top, rnd_values[9][3], mu_Zjets0, rnd_values[9][4], mu_Zjets1, rnd_values[9][5], mu_Zjets2, rnd_values[9][6], mu_WW, rnd_values[9][7]);
                  lnL_Zjets1 = lnL_Zjets11 + lnL_Zjets12;

                  lnL_Zjets21 = calc_neglnL(2871, rnd_values[10][0], mu_signal, rnd_values[10][1], mu_WZ, rnd_values[10][2], mu_top, rnd_values[10][3], mu_Zjets0, rnd_values[10][4], mu_Zjets1, rnd_values[10][5], mu_Zjets2, rnd_values[10][6], mu_WW, rnd_values[10][7]);
                  lnL_Zjets22 = calc_neglnL(2898, rnd_values[11][0], mu_signal, rnd_values[11][1], mu_WZ, rnd_values[11][2], mu_top, rnd_values[11][3], mu_Zjets0, rnd_values[11][4], mu_Zjets1, rnd_values[11][5], mu_Zjets2, rnd_values[11][6], mu_WW, rnd_values[11][7]);
                  lnL_Zjets2 = lnL_Zjets21 + lnL_Zjets22;

                  // lnL_SR = calc_neglnL(2554, rnd_values[6][0], mu_signal, rnd_values[6][1], mu_WZ, rnd_values[6][2], mu_top, rnd_values[6][3], mu_Zjets0, rnd_values[6][4], mu_Zjets1, rnd_values[6][5], mu_Zjets2, rnd_values[6][6], mu_WW, rnd_values[6][7]); //DATA = MC TOTAL

                  lnL_SR1 = calc_neglnL(1314.01, rnd_values[12][0], mu_signal, rnd_values[12][1], mu_WZ, rnd_values[12][2], mu_top, rnd_values[12][3], mu_Zjets0, rnd_values[12][4], mu_Zjets1, rnd_values[12][5], mu_Zjets2, rnd_values[12][6], mu_WW, rnd_values[12][7]); // DATA = MC TOTAL SCALED
                  lnL_SR2 = calc_neglnL(1320.53, rnd_values[13][0], mu_signal, rnd_values[13][1], mu_WZ, rnd_values[13][2], mu_top, rnd_values[13][3], mu_Zjets0, rnd_values[13][4], mu_Zjets1, rnd_values[13][5], mu_Zjets2, rnd_values[13][6], mu_WW, rnd_values[13][7]); // DATA = MC TOTAL SCALED
                  lnL_SR = lnL_SR1 + lnL_SR2;

                  lnL_total = lnL_SR + lnL_3l + lnL_emuA + lnL_emuB + lnL_Zjets0 + lnL_Zjets1 + lnL_Zjets2;
                  // lnL_total = lnL_3l + lnL_emuA + lnL_emuB + lnL_Zjets0 + lnL_Zjets1 + lnL_Zjets2;

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



    float value2 = 0.3;

    for (mu_signal = 1 - value2; mu_signal < 1 + value2; mu_signal += scan)
    {
      //                   obs   sig   wz   t   z0  z1  z2  ww  ot  sf
      lnL_3l1 = calc_neglnL(1229, rnd_values[0][0], mu_signal, rnd_values[0][1], min_mu_WZ, rnd_values[0][2], min_mu_top, rnd_values[0][3], min_mu_Zjets0, rnd_values[0][4], min_mu_Zjets1, rnd_values[0][5], min_mu_Zjets2, rnd_values[0][6], min_mu_WW, rnd_values[0][7]);
      lnL_3l2 = calc_neglnL(1180, rnd_values[1][0], mu_signal, rnd_values[1][1], min_mu_WZ, rnd_values[1][2], min_mu_top, rnd_values[1][3], min_mu_Zjets0, rnd_values[1][4], min_mu_Zjets1, rnd_values[1][5], min_mu_Zjets2, rnd_values[1][6], min_mu_WW, rnd_values[1][7]);
      lnL_3l = lnL_3l1 + lnL_3l2;

      lnL_emuA1 = calc_neglnL(1445, rnd_values[2][0], mu_signal, rnd_values[2][1], min_mu_WZ, rnd_values[2][2], min_mu_top, rnd_values[2][3], min_mu_Zjets0, rnd_values[2][4], min_mu_Zjets1, rnd_values[2][5], min_mu_Zjets2, rnd_values[2][6], min_mu_WW, rnd_values[2][7]);
      lnL_emuA2 = calc_neglnL(1458, rnd_values[3][0], mu_signal, rnd_values[3][1], min_mu_WZ, rnd_values[3][2], min_mu_top, rnd_values[3][3], min_mu_Zjets0, rnd_values[3][4], min_mu_Zjets1, rnd_values[3][5], min_mu_Zjets2, rnd_values[3][6], min_mu_WW, rnd_values[3][7]);
      lnL_emuA = lnL_emuA1 + lnL_emuA2;

      lnL_emuB1 = calc_neglnL(2867, rnd_values[4][0], mu_signal, rnd_values[4][1], min_mu_WZ, rnd_values[4][2], min_mu_top, rnd_values[4][3], min_mu_Zjets0, rnd_values[4][4], min_mu_Zjets1, rnd_values[4][5], min_mu_Zjets2, rnd_values[4][6], min_mu_WW, rnd_values[4][7]);
      lnL_emuB2 = calc_neglnL(2869, rnd_values[5][0], mu_signal, rnd_values[5][1], min_mu_WZ, rnd_values[5][2], min_mu_top, rnd_values[5][3], min_mu_Zjets0, rnd_values[5][4], min_mu_Zjets1, rnd_values[5][5], min_mu_Zjets2, rnd_values[5][6], min_mu_WW, rnd_values[5][7]);
      lnL_emuB = lnL_emuB1 + lnL_emuB2;

      lnL_Zjets01 = calc_neglnL(1355, rnd_values[6][0], mu_signal, rnd_values[6][1], min_mu_WZ, rnd_values[6][2], min_mu_top, rnd_values[6][3], min_mu_Zjets0, rnd_values[6][4], min_mu_Zjets1, rnd_values[6][5], min_mu_Zjets2, rnd_values[6][6], min_mu_WW, rnd_values[6][7]);
      lnL_Zjets02 = calc_neglnL(1336, rnd_values[7][0], mu_signal, rnd_values[7][1], min_mu_WZ, rnd_values[7][2], min_mu_top, rnd_values[7][3], min_mu_Zjets0, rnd_values[7][4], min_mu_Zjets1, rnd_values[7][5], min_mu_Zjets2, rnd_values[7][6], min_mu_WW, rnd_values[7][7]);
      lnL_Zjets0 = lnL_Zjets01 + lnL_Zjets02;

      lnL_Zjets11 = calc_neglnL(2836, rnd_values[8][0], mu_signal, rnd_values[8][1], min_mu_WZ, rnd_values[8][2], min_mu_top, rnd_values[8][3], min_mu_Zjets0, rnd_values[8][4], min_mu_Zjets1, rnd_values[8][5], min_mu_Zjets2, rnd_values[8][6], min_mu_WW, rnd_values[8][7]);
      lnL_Zjets12 = calc_neglnL(2897, rnd_values[9][0], mu_signal, rnd_values[9][1], min_mu_WZ, rnd_values[9][2], min_mu_top, rnd_values[9][3], min_mu_Zjets0, rnd_values[9][4], min_mu_Zjets1, rnd_values[9][5], min_mu_Zjets2, rnd_values[9][6], min_mu_WW, rnd_values[9][7]);
      lnL_Zjets1 = lnL_Zjets11 + lnL_Zjets12;

      lnL_Zjets21 = calc_neglnL(2871, rnd_values[10][0], mu_signal, rnd_values[10][1], min_mu_WZ, rnd_values[10][2], min_mu_top, rnd_values[10][3], min_mu_Zjets0, rnd_values[10][4], min_mu_Zjets1, rnd_values[10][5], min_mu_Zjets2, rnd_values[10][6], min_mu_WW, rnd_values[10][7]);
      lnL_Zjets22 = calc_neglnL(2898, rnd_values[11][0], mu_signal, rnd_values[11][1], min_mu_WZ, rnd_values[11][2], min_mu_top, rnd_values[11][3], min_mu_Zjets0, rnd_values[11][4], min_mu_Zjets1, rnd_values[11][5], min_mu_Zjets2, rnd_values[11][6], min_mu_WW, rnd_values[11][7]);
      lnL_Zjets2 = lnL_Zjets21 + lnL_Zjets22;

      // lnL_SR = calc_neglnL(2554, rnd_values[6][0], mu_signal, rnd_values[6][1], min_mu_WZ, rnd_values[6][2], min_mu_top, rnd_values[6][3], min_mu_Zjets0, rnd_values[6][4], min_mu_Zjets1, rnd_values[6][5], min_mu_Zjets2, rnd_values[6][6], min_mu_WW, rnd_values[6][7]); //DATA = MC TOTAL

      lnL_SR1 = calc_neglnL(1314.01, rnd_values[12][0], mu_signal, rnd_values[12][1], min_mu_WZ, rnd_values[12][2], min_mu_top, rnd_values[12][3], min_mu_Zjets0, rnd_values[12][4], min_mu_Zjets1, rnd_values[12][5], min_mu_Zjets2, rnd_values[12][6], min_mu_WW, rnd_values[12][7]); // DATA = MC TOTAL SCALED
      lnL_SR2 = calc_neglnL(1320.53, rnd_values[13][0], mu_signal, rnd_values[13][1], min_mu_WZ, rnd_values[13][2], min_mu_top, rnd_values[13][3], min_mu_Zjets0, rnd_values[13][4], min_mu_Zjets1, rnd_values[13][5], min_mu_Zjets2, rnd_values[13][6], min_mu_WW, rnd_values[13][7]); // DATA = MC TOTAL SCALED
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


    // gr1->Draw();
    // c2->SaveAs("./errors/lnL_Plot.png");

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

      tex2 = new TLatex(0.15, 0.72, "MC Statistics");
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
      hist_mu->SetFillColorAlpha(TColor::GetColor("#f4d03f"), 1.);
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

      c1->SaveAs("./errors/syst_mu.png");

      tex1 = new TLatex(0.15, 0.8, "Pull Histogram");
      tex1->SetNDC();
      tex1->SetTextSize(0.04);

      tex2 = new TLatex(0.15, 0.72, "MC Statistics");
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
      hist_pull->SetFillColorAlpha(TColor::GetColor("#f4d03f"), 1.);
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

      c2->SaveAs("./errors/syst_pull.png");
    }



  } // pseudoexperiments loop

  

    // for (int i = 0; i < lnL_total_vec.size(); i++)
    // {
    //   gr1->SetPoint(i, mu_signal_vec.at(i), lnL_total_vec.at(i));
    // }


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
    pull_var_sum += pow(pull_param_vec.at(i) - pull_mean, 2);
  }

  mu_var = mu_var_sum / min_mu_signal_vec.size();
  pull_var = pull_var_sum / pull_param_vec.size();

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
