// ROOT headers
#include <TROOT.h>
#include <TRandom3.h>
#include <TGraph.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TH1.h>
#include <TLine.h>

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
  Double_t lnL_SR1;
  Double_t lnL_3l1;
  Double_t lnL_Zjets01;
  Double_t lnL_Zjets11;
  Double_t lnL_Zjets21;
  Double_t lnL_emuA1;
  Double_t lnL_emuB1; 

  Double_t lnL_SR2;
  Double_t lnL_3l2;
  Double_t lnL_Zjets02;
  Double_t lnL_Zjets12;
  Double_t lnL_Zjets22;
  Double_t lnL_emuA2;
  Double_t lnL_emuB2; 

  Double_t lnL_total1;
  Double_t lnL_total2;
  Double_t lnL_total;
  

  

  //                    obs   sig   wz   t   z0  z1  z2  ww  ot  sf
  lnL_3l1 = calc_neglnL(1229, 0.86, mu_signal, 1081, mu_WZ, 34, mu_top, 2.4, mu_Zjets0, 20.5, mu_Zjets1, 9.2, mu_Zjets2, 0.4, mu_WW, 51.3);
  lnL_3l2 = calc_neglnL(1180, 1.42, mu_signal, 1018, mu_WZ, 23, mu_top, 33.3, mu_Zjets0, 14.0, mu_Zjets1, 6.5, mu_Zjets2, 0.3, mu_WW, 62.8);

  lnL_emuA1 = calc_neglnL(1445, 0., mu_signal, 9.0, mu_WZ, 1067.5, mu_top, 0.3, mu_Zjets0, 2.3, mu_Zjets1, 0.2, mu_Zjets2, 206.2, mu_WW, 28.2);
  lnL_emuA2 = calc_neglnL(1458, 0.4, mu_signal, 13.3, mu_WZ, 767.3, mu_top, 0.3, mu_Zjets0, 0.0, mu_Zjets1, 0.2, mu_Zjets2, 425.9, mu_WW, 68.7);

  lnL_emuB1 = calc_neglnL(2867, 0, mu_signal, 0.4, mu_WZ, 2860.7, mu_top, 0, mu_Zjets0, 0.0, mu_Zjets1, 0.1, mu_Zjets2, 6.1, mu_WW, 0.9);
  lnL_emuB2 = calc_neglnL(2869, 0, mu_signal, 0.5, mu_WZ, 2793, mu_top, 0, mu_Zjets0, 0.1, mu_Zjets1, 0.1, mu_Zjets2, 7.0, mu_WW, 3.6);

  lnL_Zjets01 = calc_neglnL(1355, 37, mu_signal, 40.5, mu_WZ, 63.1, mu_top, 855.7, mu_Zjets0, 0, mu_Zjets1, 0, mu_Zjets2, 19.7, mu_WW, 7.6);
  lnL_Zjets02 = calc_neglnL(1336, 60.8, mu_signal, 38.9, mu_WZ, 27.9, mu_top, 874.4, mu_Zjets0, 0, mu_Zjets1, 0, mu_Zjets2, 28.0, mu_WW, 4.8);

  lnL_Zjets11 = calc_neglnL(2836, 86.5, mu_signal, 144.4, mu_WZ, 144.0, mu_top, 0, mu_Zjets0, 1790.6, mu_Zjets1, 0, mu_Zjets2, 34.8, mu_WW, 18.4);
  lnL_Zjets12 = calc_neglnL(2897, 58.5, mu_signal, 103.9, mu_WZ, 52.5, mu_top, 0, mu_Zjets0, 2007.9, mu_Zjets1, 0, mu_Zjets2, 11.2, mu_WW, 27.7);

  lnL_Zjets21 = calc_neglnL(2871, 83.6, mu_signal, 158.6, mu_WZ, 152.3, mu_top, 0, mu_Zjets0, 0, mu_Zjets1, 2296.7, mu_Zjets2, 12.4, mu_WW, 24.8);
  lnL_Zjets22 = calc_neglnL(2898, 78.3, mu_signal, 150.0, mu_WZ, 91.7, mu_top, 0, mu_Zjets0, 0, mu_Zjets1, 2362.6, mu_Zjets2, 5.9, mu_WW, 33.0);

  lnL_SR1 = calc_neglnL(1314.01, 713.5, mu_signal, 356.2, mu_WZ, 71.8, mu_top, 51.8, mu_Zjets0, 31.2, mu_Zjets1, 0.9, mu_Zjets2, 19.7, mu_WW, 22.8); // DATA = MC SCALED
  lnL_SR2 = calc_neglnL(1320.53, 863.6, mu_signal, 312.9, mu_WZ, 19.1, mu_top, 50.1, mu_Zjets0, 0.5, mu_Zjets1, 0.2, mu_Zjets2, 17.7, mu_WW, 22.1); // DATA = MC SCALED


  lnL_total1 = lnL_SR1 + lnL_3l1 + lnL_emuA1 + lnL_emuB1 + lnL_Zjets01 + lnL_Zjets11 + lnL_Zjets21;
  lnL_total2 = lnL_SR2 + lnL_3l2 + lnL_emuA2 + lnL_emuB2 + lnL_Zjets02 + lnL_Zjets12 + lnL_Zjets22;

  lnL_total = lnL_total1 + lnL_total2;

  return lnL_total;
}

void mle_simfit()
{

  gROOT->SetBatch(kTRUE);

  auto start = std::chrono::high_resolution_clock::now();

  // Output log file
  ofstream logFile("./mle_simfit.txt");

  DualStreamBuffer dualBuffer(cout.rdbuf(), logFile.rdbuf());

  std::streambuf *oldBuffer = cout.rdbuf(&dualBuffer);

  Double_t lnL_SR1;
  Double_t lnL_3l1;
  Double_t lnL_Zjets01;
  Double_t lnL_Zjets11;
  Double_t lnL_Zjets21;
  Double_t lnL_emuA1;
  Double_t lnL_emuB1;

  Double_t lnL_SR2;
  Double_t lnL_3l2;
  Double_t lnL_Zjets02;
  Double_t lnL_Zjets12;
  Double_t lnL_Zjets22;
  Double_t lnL_emuA2;
  Double_t lnL_emuB2; 
  
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

  Double_t step;

  //----------------------------------

  //BEST FIT CENTRAL VALUES AT 0.001 PRECISION -
  Double_t value = 0.003;

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

  step = 0.001;


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

  //Statistical uncertainty from likelihood shape

  float sf_low = 1 - 0.8;
  float sf_up  = 1 + 0.8;
  float sf;
  float sf_min;
  step = 0.00001;
  iter = 0;
  
  for (sf = sf_low; sf <= sf_up; sf += step)
  {

    //mu_signal,  mu_WZ,  mu_top,  mu_Zjets0,  mu_Zjets1,  mu_Zjets2,  mu_WW
    lnL_total = calc_total_neglnL(sf, 1.012, 1.013, 1.352, 1.322, 1.065, 1.435);     // SIM ALL REGIONS
    
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
  TGraph *gr1 = new TGraph();


  for (int i = 0; i < lnL_total_vec.size(); i++)
  {
    gr1->SetPoint(i, sf_vec.at(i), lnL_total_vec.at(i));
  }

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


  cout << "   lnL:           " << lnL_min << setprecision(10) << endl << endl;
  cout << "   sf_min :       " << sf_min << endl << endl;
  cout << "   Iter best:     " << iter_best << endl << endl;
  cout << "   Iter plus:     " << iter_er_plus  << "  ->  +" << sf_vec.at(iter_er_plus) - sf_vec.at(iter_best) << endl << endl;
  cout << "   Iter minus:    " << iter_er_minus << "  ->  -" << sf_vec.at(iter_best) - sf_vec.at(iter_er_minus) << endl << endl;
  
  c1->SetLeftMargin(0.15);
  
  gr1->GetXaxis()->SetTitleFont(62);
  gr1->GetXaxis()->SetTitleSize(0.04);
  gr1->GetXaxis()->SetTitle("#mu_{S}");
  
  gr1->GetYaxis()->SetTitle("-lnL");
  gr1->GetYaxis()->SetTitleOffset(1.5);
  gr1->GetYaxis()->SetTitleSize(0.05);
  // gr1->GetYaxis()->SetNdivisions(0);

  gr1->SetLineColor(kRed);
  gr1->SetLineWidth(3);
  gr1->Draw("AL");


  // TLine *line1 = new TLine(xbins[0], 1, xbins[sizeof(xbins) / sizeof(xbins[0]) - 1], 1);
  // line1->SetLineStyle(2);
  // line1->SetLineWidth(2);
  // line1->Draw("same");

  c1->SetTicks();
  c1->SetGrid();
  c1->SaveAs("./MLE_plot.png");



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
