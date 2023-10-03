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



void simfit_systematic()
{

  gROOT->SetBatch(kTRUE);

  auto start = std::chrono::high_resolution_clock::now();

  // Output log file
  ofstream logFile("./simfit_systematic.txt");

  DualStreamBuffer dualBuffer(cout.rdbuf(), logFile.rdbuf());

  std::streambuf *oldBuffer = cout.rdbuf(&dualBuffer);



  TCanvas *c1 = new TCanvas("c1", "-lnL vs mu", 800, 600);
  TCanvas *c2 = new TCanvas("c2", "mu distribution", 800, 600);
  TCanvas *c3 = new TCanvas("c3", "pull parameter ", 800, 600);
  
  
  TH1F *hist_mu = new TH1F("hist_mu", "mu hist", 200, 0.9,1.1 );
  TH1F *hist_pull = new TH1F("hist_pull", "pull hist", 200, -4, 4 );
  TGraph *gr1 = new TGraph();


  Double_t lnL_SR;
  Double_t lnL_3l;
  Double_t lnL_Zjets0;
  Double_t lnL_Zjets1;
  Double_t lnL_Zjets2;
  Double_t lnL_emuA;
  Double_t lnL_emuB;
  Double_t lnL_WZ;
  Double_t lnL_total;


  Double_t lnL_min;
  int mod = 1000;

  // Scaling Factors

  Double_t mu_signal = 1;
  Double_t mu_top = 1;
  Double_t mu_WZ = 1;
  Double_t mu_Zjets0 = 1;
  Double_t mu_Zjets1 = 1;
  Double_t mu_Zjets2 = 1;
  Double_t mu_WW = 1;

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




  // MC Systematic Uncertainty

  
  double rnd_values[7][8];

  double central_values[7][8] = {
      {2.3, 2099, 57, 35.7, 34.5, 15.7, 0.6, 114.1},    // 3l
      {0.4, 22.3, 1834.8, 0.6, 2.2, 0.4, 632.1, 96.9},    // emuA
      {0, 0.9, 5653.7, 0, 0.1, 0.1, 13, 4.6},           // emuB
      {97.7, 79.4, 91, 1730.1, 0, 0, 47.7, 12.4},       // Zjets0
      {145, 248.2, 196.5, 0, 3798.5, 0, 46, 46.1},      // Zjets1
      {161.9, 308.7, 244.1, 0, 0, 4659.2, 18.4, 57.9},  // Zjets2
      {1577, 669.2, 90.9, 101.9, 31.6, 1.2, 37.3, 44.9} // SR
  };


  double uncertainties[7][8] = {
      {0.5, 10.2, 1.8, 8.3, 5.8, 4.1, 0.1, 0.9}, // 3l
      {0.2, 1, 11.1, 0.4, 1.3, 0.4, 4.5, 16.4},  // emuA
      {0, 0.2, 17.6, 0, 0.1, 0.1, 0.7, 2},       // emuB
      {3.6, 2, 2.7, 99.6, 0, 0, 1.2, 2.4},       // Zjets0
      {4, 3.5, 3.7, 0, 89.3, 0, 1.2, 2.5},       // Zjets1
      {2.5, 3., 3.8, 0, 0, 110.9, 0.8, 2.6},     // Zjets2
      {0.2, 6.1, 2.5, 23.4, 8.7, 0.7, 1.1, 2},  // SR
  };

  mu_signal = 1.006;
  mu_WZ = 1.012;
  mu_top = 1.010;
  mu_Zjets0 = 1.351;
  mu_Zjets1 = 1.323;
  mu_Zjets2 = 1.065;
  mu_WW = 1.465;
  
  vector<Double_t> min_mu_signal_vec;
  vector<Double_t> error_minus_vec;
  vector<Double_t> error_plus_vec;

  Double_t step = 0.0001;


  TRandom3 rnd_value;



  Double_t mu_sum = 0;
  Double_t mu_sum_sq = 0;
  Double_t mu_mean = 0;
  Double_t mu_var = 0;
  Double_t mu_var_sum = 0;
  Double_t pull_param = 0;

  int experiments = 20000;

  for (int k = 0; k < experiments; k++)
  {
    vector<Double_t> mu_signal_vec;
    vector<Double_t> lnL_total_vec;

 

    if (k % mod == 0)
    {
      cout << "   iteration " << k / mod << "000"  << endl << endl;
    }

    for (int i = 0; i < 7; i++)
    {
      for (int j = 0; j < 8; j++)
      {
        rnd_values[i][j] = rnd_value.Gaus(central_values[i][j], uncertainties[i][j]);
      }
    }

    iter = 0;

    for (mu_signal = 0.9; mu_signal < 1.1; mu_signal += step)
    {
      

      //                   obs   sig   wz   t   z0  z1  z2  ww  ot  sf
      // lnL_SR = calc_neglnL(2891, 1577, mu_signal, 669, mu_WZ, 91, mu_top, 102, mu_Zjets0, 32, mu_Zjets1, 1, mu_Zjets2, 37, mu_WW, 45); //DATA UNBLINDED

      lnL_3l = calc_neglnL(2409, rnd_values[0][0], mu_signal, rnd_values[0][1], mu_WZ, rnd_values[0][2], mu_top, rnd_values[0][3], mu_Zjets0, rnd_values[0][4], mu_Zjets1, rnd_values[0][5], mu_Zjets2, rnd_values[0][6], mu_WW, rnd_values[0][7]);

      lnL_emuA = calc_neglnL(2903, rnd_values[1][0], mu_signal, rnd_values[1][1], mu_WZ, rnd_values[1][2], mu_top, rnd_values[1][3], mu_Zjets0, rnd_values[1][4], mu_Zjets1, rnd_values[1][5], mu_Zjets2, rnd_values[1][6], mu_WW, rnd_values[1][7]);

      lnL_emuB = calc_neglnL(5736, rnd_values[2][0], mu_signal, rnd_values[2][1], mu_WZ, rnd_values[2][2], mu_top, rnd_values[2][3], mu_Zjets0, rnd_values[2][4], mu_Zjets1, rnd_values[2][5], mu_Zjets2, rnd_values[2][6], mu_WW, rnd_values[2][7]);

      lnL_Zjets0 = calc_neglnL(2691, rnd_values[3][0], mu_signal, rnd_values[3][1], mu_WZ, rnd_values[3][2], mu_top, rnd_values[3][3], mu_Zjets0, rnd_values[3][4], mu_Zjets1, rnd_values[3][5], mu_Zjets2, rnd_values[3][6], mu_WW, rnd_values[3][7]);

      lnL_Zjets1 = calc_neglnL(5733, rnd_values[4][0], mu_signal, rnd_values[4][1], mu_WZ, rnd_values[4][2], mu_top, rnd_values[4][3], mu_Zjets0, rnd_values[4][4], mu_Zjets1, rnd_values[4][5], mu_Zjets2, rnd_values[4][6], mu_WW, rnd_values[4][7]);

      lnL_Zjets2 = calc_neglnL(5769, rnd_values[5][0], mu_signal, rnd_values[5][1], mu_WZ, rnd_values[5][2], mu_top, rnd_values[5][3], mu_Zjets0, rnd_values[5][4], mu_Zjets1, rnd_values[5][5], mu_Zjets2, rnd_values[5][6], mu_WW, rnd_values[5][7]);

      lnL_SR = calc_neglnL(2636, rnd_values[6][0], mu_signal, rnd_values[6][1], mu_WZ, rnd_values[6][2], mu_top, rnd_values[6][3], mu_Zjets0, rnd_values[6][4], mu_Zjets1, rnd_values[6][5], mu_Zjets2, rnd_values[6][6], mu_WW, rnd_values[6][7]); //DATA = MC TOTAL


      lnL_total = lnL_SR + lnL_3l + lnL_emuA + lnL_emuB + lnL_Zjets0 + lnL_Zjets1 + lnL_Zjets2;
      // lnL_total = lnL_3l + lnL_emuA + lnL_emuB + lnL_Zjets0 + lnL_Zjets1 + lnL_Zjets2;

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



    min_mu_signal_vec.push_back(min_mu_signal);

    for (int i = 0; i < lnL_total_vec.size(); i++)
    {
      gr1->SetPoint(i, mu_signal_vec.at(i), lnL_total_vec.at(i));
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
    
    //CDF pulls - p.8
    // pull_param = (min_mu_signal - 1.006) / (1.006/sqrt(experiments));
    // pull_param = (min_mu_signal - 1.006) / ((0.0326 - (-0.0321))/2);
    // pull_param = (min_mu_signal - 1.006) / (0.023);
    pull_param = (min_mu_signal - 1.006) / ((mu_signal_vec.at(iter_er_plus) - mu_signal_vec.at(iter_er_minus))/2);
    // pull_param = (min_mu_signal - 1.005) / ((0.09 - (-0.089)) / 2);

    hist_mu->Fill(min_mu_signal);
    hist_pull->Fill(pull_param);

    // cout << "   Iter minus:  " << iter_er_minus << "  ->  -" << mu_signal_vec.at(ibest) - mu_signal_vec.at(iter_er_minus) << endl
    //      << endl;
    // cout << "   Iter plus:   " << iter_er_plus << "  ->  +" << mu_signal_vec.at(iter_er_plus) - mu_signal_vec.at(ibest) << endl
    //      << endl;

  }

  cout << "RMS  "  << hist_pull->GetRMS() << endl;




  c1->cd();
  hist_mu->Draw();

  c1->SaveAs("./mu_plot.png");


  c2->cd();
  hist_pull->Draw();
  // gr1->Draw("ALP");

  c2->SaveAs("./pull_plot.png");


  for (int i = 0; i < min_mu_signal_vec.size(); i++)
  {
    mu_sum += min_mu_signal_vec.at(i);
    mu_sum_sq += pow(min_mu_signal_vec.at(i), 2);
  }

  mu_mean = mu_sum / min_mu_signal_vec.size();


  for (int i = 0; i < min_mu_signal_vec.size(); i++)
  {
    mu_var_sum += pow(min_mu_signal_vec.at(i) - mu_mean, 2);
  }

  mu_var = mu_var_sum / min_mu_signal_vec.size();


  cout << "   SYSTEMATIC " <<  sqrt(mu_var) << endl << endl;    //SYSTEMATIC DUE TO MC






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
