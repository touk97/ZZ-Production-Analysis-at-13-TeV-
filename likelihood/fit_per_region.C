// ROOT headers
#include <TROOT.h>
#include <TRandom3.h>
#include <TGraph.h>
#include <TMath.h>
#include <TCanvas.h>

// cpp headers
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
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

struct scaling_factors_min
{
  Double_t signal;
  Double_t top;
  Double_t WZ;
  Double_t Zjets0;
  Double_t Zjets1;
  Double_t Zjets2;
  Double_t WW;
};

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

  // lnL_total = lnL_SR + lnL_3l + lnL_emuA + lnL_emuB + lnL_Zjets0 + lnL_Zjets1 + lnL_Zjets2;
  lnL_total = lnL_3l + lnL_emuA + lnL_emuB + lnL_Zjets0 + lnL_Zjets1 + lnL_Zjets2;

  return lnL_total;
}



void print_results(Double_t min_mu_signal, Double_t min_mu_WZ, Double_t min_mu_top, Double_t min_mu_Zjets0, Double_t min_mu_Zjets1, Double_t min_mu_Zjets2, Double_t min_mu_WW, Double_t lnL_min)
{

  cout << "   ----------------------------------------------------------------------------" << endl << endl;
  cout << "   FIT PER REGION:                    " << endl;
  cout << "   ________________________      " << endl << endl;
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


void fit_per_region()
{

  gROOT->SetBatch(kTRUE);

  auto start = std::chrono::high_resolution_clock::now();

  // Output log file
  ofstream logFile("./likelihood_fit.txt");

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

  Double_t lnL_min;
  int mod = 10000000;
  int min_i;

  int only_CRs;

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

  // Double_t value = 0.2;

  // // Min
  // Double_t low_signal = 1.17 - value;
  // Double_t low_top = 1.01 - value;
  // Double_t low_WZ = 1.01 - value;
  // Double_t low_Zjets0 = 1.35 - value;
  // Double_t low_Zjets1 = 1.3 - value;
  // Double_t low_Zjets2 = 1.04 - value;
  // Double_t low_WW = 1.45 - value;

  // // Max
  // Double_t up_signal = 1.17 + value;
  // Double_t up_top = 1.01 + value;
  // Double_t up_WZ = 1.01 + value;
  // Double_t up_Zjets0 = 1.35 + value;
  // Double_t up_Zjets1 = 1.3 + alue;
  // Double_t up_Zjets2 = 1.04 + value;
  // Double_t up_WW = 1.45 + value;

  Double_t value = 1;

  // Min
  Double_t low_signal = 1.1 - value;
  Double_t low_WZ = 1.2 - value;
  Double_t low_WW = 1.2 - value;
  Double_t low_top = 1.1 - value;
  Double_t low_Zjets0 = 1.4 - value;
  Double_t low_Zjets1 = 1.3 - value;
  Double_t low_Zjets2 = 1. - value;

  // Max
  Double_t up_signal = 1.1 + value;
  Double_t up_WZ = 1.2 + value;
  Double_t up_WW = 1.2 + value;
  Double_t up_top = 1.1 + value;
  Double_t up_Zjets0 = 1.4 + value;
  Double_t up_Zjets1 = 1.3 + value;
  Double_t up_Zjets2 = 1. + value;

  long int iter = 0;
  Double_t step = 0.00001;

  vector<string> sfnames = {
      "top",
      "WZ",
      "Zjets2",
      "Zjets1",
      "Zjets0",
      "WW",
      "signal"};



  for (mu_top = low_top; mu_top <= up_top; mu_top += step)
  {

   lnL_emuB = calc_neglnL(5736, 0, mu_signal, 0.9, mu_WZ, 5653.7, mu_top, 0, mu_Zjets0, 0.1, mu_Zjets1, 0.1, mu_Zjets2, 13, mu_WW, 4.6);

    if (lnL_emuB <= lnL_min or iter == 0)
    {
      if (lnL_emuB == lnL_min)
      {
        cout << "    SAME VALUE!!   " << endl
             << endl;
      }

      lnL_min = lnL_emuB;
      min_mu_top = mu_top;
    }

    iter++;
  }

  iter = 0;
  mu_top = min_mu_top;


  

  for (mu_WZ = low_WZ; mu_WZ <= up_WZ; mu_WZ += step)
  {
    lnL_3l = calc_neglnL(2409, 2.3, mu_signal, 2099, mu_WZ, 57, mu_top, 35.7, mu_Zjets0, 34.5, mu_Zjets1, 15.7, mu_Zjets2, 0.6, mu_WW, 114.1);

    if (lnL_3l <= lnL_min or iter == 0)
    {
      if (lnL_3l == lnL_min)
      {
        cout << "    SAME VALUE!!   " << endl
             << endl;
      }

      lnL_min = lnL_3l;
      min_mu_WZ = mu_WZ;
    }

    iter++;
  }

  iter = 0;
  mu_WZ = min_mu_WZ;

  for (mu_Zjets2 = low_Zjets2; mu_Zjets2 <= up_Zjets2; mu_Zjets2 += step)
  {
    lnL_Zjets2 = calc_neglnL(5769, 161.9, mu_signal, 308.7, mu_WZ, 244.1, mu_top, 0, mu_Zjets0, 0, mu_Zjets1, 4659.2, mu_Zjets2, 18.4, mu_WW, 57.9);
    
    if (lnL_Zjets2 <= lnL_min or iter == 0)
    {
      if (lnL_Zjets2 == lnL_min)
      {
        cout << "    SAME VALUE!!   " << endl
             << endl;
      }

      lnL_min = lnL_Zjets2;
      min_mu_Zjets2 = mu_Zjets2;
    }

    iter++;
  }

  iter = 0;
  mu_Zjets2 = min_mu_Zjets2;

  for (mu_Zjets1 = low_Zjets1; mu_Zjets1 <= up_Zjets1; mu_Zjets1 += step)
  {

    lnL_Zjets1 = calc_neglnL(5733, 145, mu_signal, 248.2, mu_WZ, 196.5, mu_top, 0, mu_Zjets0, 3798.5, mu_Zjets1, 0, mu_Zjets2, 46, mu_WW, 46.1);


    if (lnL_Zjets1 <= lnL_min or iter == 0)
    {
      if (lnL_Zjets1 == lnL_min)
      {
        cout << "    SAME VALUE!!   " << endl
             << endl;
      }

      lnL_min = lnL_Zjets1;
      min_mu_Zjets1 = mu_Zjets1;
    }

    iter++;
  }

  iter = 0;
  mu_Zjets1 = min_mu_Zjets1;

  for (mu_Zjets0 = low_Zjets0; mu_Zjets0 <= up_Zjets0; mu_Zjets0 += step)
  {

    lnL_Zjets0 = calc_neglnL(2691, 97.7, mu_signal, 79.4, mu_WZ, 91, mu_top, 1730.1, mu_Zjets0, 0, mu_Zjets1, 0, mu_Zjets2, 47.7, mu_WW, 12.4);

    if (lnL_Zjets0 <= lnL_min or iter == 0)
    {
      if (lnL_Zjets0 == lnL_min)
      {
        cout << "    SAME VALUE!!   " << endl
             << endl;
      }

      lnL_min = lnL_Zjets0;
      min_mu_Zjets0 = mu_Zjets0;
    }

    iter++;
  }

  iter = 0;
  mu_Zjets0 = min_mu_Zjets0;

  for (mu_WW = low_WW; mu_WW <= up_WW; mu_WW += step)
  {

    lnL_emuA = calc_neglnL(2903, 0.4, mu_signal, 22.3, mu_WZ, 1834.8, mu_top, 0.6, mu_Zjets0, 2.2, mu_Zjets1, 0.4, mu_Zjets2, 632.1, mu_WW, 96.9);  

    if (lnL_emuA <= lnL_min or iter == 0)
    {
      if (lnL_emuA == lnL_min)
      {
        cout << "    SAME VALUE!!   " << endl
             << endl;
      }

      lnL_min = lnL_emuA;
      min_mu_WW = mu_WW;
    }

    iter++;
  }

  iter = 0;
  mu_WW = min_mu_WW;

  for (mu_signal = low_signal; mu_signal <= up_signal; mu_signal += step)
  {

    lnL_SR = calc_neglnL(2636, 1577, mu_signal, 669.2, mu_WZ, 90.9, mu_top, 101.9, mu_Zjets0, 31.6, mu_Zjets1, 1.2, mu_Zjets2, 37.3, mu_WW, 44.9); 

    cout << lnL_SR << endl << endl;

    if (lnL_SR <= lnL_min or iter == 0)
    {
      if (lnL_SR == lnL_min)
      {
        cout << "    SAME VALUE!!   " << endl
             << endl;
      }

      lnL_min = lnL_SR;
      min_mu_signal = mu_signal;
    }

    iter++;
  }

  iter = 0;
  mu_signal = min_mu_signal;


print_results(min_mu_signal, min_mu_WZ, min_mu_top, min_mu_Zjets0, min_mu_Zjets1, min_mu_Zjets2, min_mu_WW, lnL_min);



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
