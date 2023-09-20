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


struct scaling_factors
{
  float signal;
  float top;
  float WZ;
  float Zjets0;
  float Zjets1;
  float Zjets2;
  float WW;
};

struct scaling_factors_min
{
  float signal;
  float top;
  float WZ;
  float Zjets0;
  float Zjets1;
  float Zjets2;
  float WW;
};

struct bounds_min
{
  float value = 0.1;
  float signal = 1.06 - value;
  float top = 1.01 - value;
  float WZ = 1.02 - value;
  float Zjets0 = 1.40 - value;
  float Zjets1 = 1.32 - value;
  float Zjets2 = 1.05 - value;
  float WW = 1.46 - value;
};

struct bounds_max
{
  float value = 0.1;
  float signal = 1.06 + value;
  float top = 1.01 + value;
  float WZ = 1.02 + value;
  float Zjets0 = 1.40 + value;
  float Zjets1 = 1.32 + value;
  float Zjets2 = 1.05 + value;
  float WW = 1.46 + value;
};

float calc_neglnL(float n_obs, float n_signal, float n_WZ, float n_top, float n_Zjets0, float n_Zjets1, float n_Zjets2, float n_WW, float n_other, struct scaling_factors sf)
{
  float lnLi;
  float n_exp;
  n_exp = sf.signal * n_signal + sf.top * n_top + sf.WZ * n_WZ + sf.Zjets2 * n_Zjets2 + sf.Zjets1 * n_Zjets1 + sf.Zjets0 * n_Zjets0 + sf.WW * n_WW + n_other;

  lnLi = +n_obs * log(n_exp) - n_exp;

  return -lnLi;
}

float calc_neglnL2(float n_obs, float n_signal, float n_WZ, float n_top, float n_Zjets0, float n_Zjets1, float n_Zjets2, float n_WW, float n_other, float mu_signal_only, struct scaling_factors_min sf_min)
{
  float lnLi;
  float n_exp;
  n_exp = mu_signal_only * n_signal + sf_min.top * n_top + sf_min.WZ * n_WZ + sf_min.Zjets2 * n_Zjets2 + sf_min.Zjets1 * n_Zjets1 + sf_min.Zjets0 * n_Zjets0 + sf_min.WW * n_WW + n_other;

  lnLi = +n_obs * log(n_exp) - n_exp;

  return -lnLi;
}



void likelihood_fit()
{

  gROOT->SetBatch(kTRUE);

  auto start = std::chrono::high_resolution_clock::now();

  //Output log file
  ofstream logFile("./likelihood_fit.txt");

  DualStreamBuffer dualBuffer(std::cout.rdbuf(), logFile.rdbuf());

  std::streambuf *oldBuffer = std::cout.rdbuf(&dualBuffer);

  scaling_factors sf;
  scaling_factors_min sf_min;
  bounds_min bmin;
  bounds_max bmax;

  float lnL_SR;
  float lnL_3l;
  float lnL_Zjets0;
  float lnL_Zjets1;
  float lnL_Zjets2;
  float lnL_emuA;
  float lnL_emuB;
  float lnL_WZ;
  float lnL_total;
 


  float step = 0.001;
  float lnL_min;
  int mod = 100000000;
  int min_i;
  int i = 0;

  sf.WW = 1;
  sf.WZ = 1;
  sf.Zjets0 = 1;
  sf.Zjets1 = 1;
  sf.Zjets2 = 1;
  sf.signal = 1;


  for (sf.top = bmin.top; sf.top <= bmax.top; sf.top += step)
  {
    lnL_emuB = calc_neglnL(5736, 0, 1, 5654, 0, 0, 0, 13, 5, sf);

    if (lnL_emuB < lnL_min or i == 0)
    {
      lnL_min = lnL_emuB;
      sf_min.top = sf.top;
    }
    i++;
  }

  i = 0;
  sf.top = sf_min.top;

  for (sf.WZ = bmin.WZ; sf.WZ <= bmax.WZ; sf.WZ += step)
  {
    lnL_3l = calc_neglnL(2409, 2, 2099, 57, 36, 34, 16, 1, 114, sf);

    if (lnL_3l < lnL_min or i == 0)
    {
      lnL_min = lnL_3l;
      sf_min.WZ = sf.WZ;
    }
    i++;
  }
  
  i = 0;
  sf.WZ = sf_min.WZ;

  for (sf.Zjets2 = bmin.Zjets2; sf.Zjets2 <= bmax.Zjets2; sf.Zjets2 += step)
  {
    lnL_Zjets2 = calc_neglnL(2409, 2, 2099, 57, 36, 34, 16, 1, 114, sf);

    if (lnL_Zjets2 < lnL_min or i == 0)
    {
      lnL_min = lnL_3l;
      sf_min.Zjets2 = sf.Zjets2;
    }
    i++;
  }
  
  i = 0;
  sf.Zjets2 = sf_min.Zjets2;

  for (sf.Zjets1 = bmin.Zjets1; sf.Zjets1 <= bmax.Zjets1; sf.Zjets1 += step)
  {
    lnL_Zjets1 = calc_neglnL(2409, 2, 2099, 57, 36, 34, 16, 1, 114, sf);

    if (lnL_Zjets1 < lnL_min or i == 0)
    {
      lnL_min = lnL_3l;
      sf_min.Zjets1 = sf.Zjets1;
    }
    i++;
  }
  
  i = 0;
  sf.Zjets1 = sf_min.Zjets1;

  for (sf.Zjets0 = bmin.Zjets0; sf.Zjets0 <= bmax.Zjets0; sf.Zjets0 += step)
  {
    lnL_Zjets0 = calc_neglnL(2409, 2, 2099, 57, 36, 34, 16, 1, 114, sf);

    if (lnL_Zjets0 < lnL_min or i == 0)
    {
      lnL_min = lnL_3l;
      sf_min.Zjets0 = sf.Zjets0;
    }
    i++;
  }
  
  i = 0;
  sf.Zjets0 = sf_min.Zjets0;

  for (sf.WW = bmin.WW; sf.WW <= bmax.WW; sf.WW += step)
  {
    lnL_emuA = calc_neglnL(2409, 2, 2099, 57, 36, 34, 16, 1, 114, sf);

    if (lnL_emuA < lnL_min or i == 0)
    {
      lnL_min = lnL_3l;
      sf_min.WW = sf.WW;
    }
    i++;
  }
  
  i = 0;
  sf.WW = sf_min.WW;

  for (sf.signal = bmin.signal; sf.signal <= bmax.signal; sf.signal += step)
  {
    lnL_SR = calc_neglnL(2409, 2, 2099, 57, 36, 34, 16, 1, 114, sf);

    if (lnL_SR < lnL_min or i == 0)
    {
      lnL_min = lnL_3l;
      sf_min.signal = sf.signal;
    }
    i++;
  }
  
  i = 0;
  sf.signal = sf_min.signal;


  cout << "   ----------------------------------------------------------------------------" << endl << endl;
  cout << "   SUCCESSIVE FIT:        " << endl;
  cout << "   _______________________  " << endl << endl;
  cout << "   lnL =          " << lnL_min << endl << endl;
  cout << "   mu_signal =    " << sf_min.signal << endl << endl;
  cout << "   mu_WZ =        " << sf_min.WZ << endl << endl;
  cout << "   mu_WW =        " << sf_min.WW << endl << endl;
  cout << "   mu_top =       " << sf_min.top << endl << endl;
  cout << "   mu_Zjets0 =    " << sf_min.Zjets0 << endl << endl;
  cout << "   mu_Zjets1 =    " << sf_min.Zjets1 << endl << endl;
  cout << "   mu_Zjets2 =    " << sf_min.Zjets2 << endl << endl;
  cout << "   ----------------------------------------------------------------------------" << endl << endl;


  step = 0.01;
  

  for (sf.signal = bmin.signal; sf.signal <= bmax.signal; sf.signal += step)
  {
    for (sf.WW = bmin.WW; sf.WW <= bmax.WW; sf.WW += step)
    {
      for (sf.Zjets0 = bmin.Zjets0; sf.Zjets0 <= bmax.Zjets0; sf.Zjets0 += step)
      {
        for (sf.Zjets1 = bmin.Zjets1; sf.Zjets1 <= bmax.Zjets1; sf.Zjets1 += step)
        {
          for (sf.Zjets2 = bmin.Zjets2; sf.Zjets2 <= bmax.Zjets2; sf.Zjets2 += step)
          {
            for (sf.WZ = bmin.WZ; sf.WZ <= bmax.WZ; sf.WZ += step)
            {
              for (sf.top = bmin.top; sf.top <= bmax.top; sf.top += step)
              {

                //                       obs   sig   wz   t   z0  z1  z2  ww  ot  sf
                lnL_SR = calc_neglnL(2748, 1554, 692, 88, 78, 63, 10, 32, 50, sf);

                lnL_3l = calc_neglnL(2409, 2, 2099, 57, 36, 34, 16, 1, 114, sf);

                lnL_emuA = calc_neglnL(2903, 0, 22, 1835, 1, 2, 0, 632, 97, sf);

                lnL_emuB = calc_neglnL(5736, 0, 1, 5654, 0, 0, 0, 13, 5, sf);

                lnL_Zjets0 = calc_neglnL(6202, 723, 442, 289, 3254, 0, 0, 124, 50, sf);

                lnL_Zjets1 = calc_neglnL(8175, 448, 593, 425, 0, 4940, 0, 87, 78, sf);

                lnL_Zjets2 = calc_neglnL(6710, 336, 499, 343, 0, 0, 5151, 27, 77, sf);

                lnL_total = lnL_SR + lnL_3l + lnL_emuA + lnL_emuB + lnL_Zjets0 + lnL_Zjets1 + lnL_Zjets2;
                // lnL_total = lnL_3l + lnL_emuA + lnL_emuB + lnL_Zjets0 + lnL_Zjets1 + lnL_Zjets2;

                if (i % mod == 0)
                {
                  cout << "   " << i/mod << "00 mil)  " << 
                  // "Signal:" << lnL_SR << " WZ: " << lnL_WZ << " WW: " << lnL_emuA << " Top: " << lnL_emuB << " Zjets0:  " << lnL_Zjets0 << 
                  // " Zjets1:  " << lnL_Zjets1 << " Zjets2:  " << lnL_Zjets2 << " Total: " << lnL_total << 
                  " mu_signal = " << sf.signal << "|  mu_WZ = " << sf.WZ << "|  mu_WW = " << sf.WW << "|  mu_top = " << sf.top << 
                  " mu_Zjets0 = " << sf.Zjets0 << "|  mu_Zjets1 =  " << sf.Zjets1 << "|  mu_Zjets2 = " << sf.Zjets2 << endl << endl;
                }

                if (lnL_total < lnL_min or i == 0)
                {
                  lnL_min = lnL_total;
                  min_i = i;

                  sf_min.signal = sf.signal;
                  sf_min.WZ = sf.WZ;
                  sf_min.WW = sf.WW;
                  sf_min.top = sf.top;
                  sf_min.Zjets0 = sf.Zjets0;
                  sf_min.Zjets1 = sf.Zjets1;
                  sf_min.Zjets2 = sf.Zjets2;
                }
  
                i++;
              }
            }
          }
        }
      }
    }
  }

  

  cout << "   ----------------------------------------------------------------------------" << endl << endl;
  cout << "   SIMULTANEOUS FIT:        " << endl;
  cout << "   _______________________  " << endl << endl;
  cout << "   lnL =          " << lnL_min << endl << endl;
  cout << "   mu_signal =    " << sf_min.signal << endl << endl;
  cout << "   mu_WZ =        " << sf_min.WZ << endl << endl;
  cout << "   mu_WW =        " << sf_min.WW << endl << endl;
  cout << "   mu_top =       " << sf_min.top << endl << endl;
  cout << "   mu_Zjets0 =    " << sf_min.Zjets0 << endl << endl;
  cout << "   mu_Zjets1 =    " << sf_min.Zjets1 << endl << endl;
  cout << "   mu_Zjets2 =    " << sf_min.Zjets2 << endl << endl;
  cout << "   ----------------------------------------------------------------------------" << endl << endl;

  //For successive fit

  step = 0.001;
  float mu_signal_only;
  lnL_min = 9999;
  i = 0;
  


  for (mu_signal_only = 0.7; mu_signal_only <= 1.5; mu_signal_only += step)
  {

    lnL_SR = calc_neglnL2(2748, 1554, 692, 88, 78, 63, 10, 32, 50, mu_signal_only, sf_min);

    lnL_3l = calc_neglnL2(2409, 2, 2099, 57, 36, 34, 16, 1, 114, mu_signal_only, sf_min);

    lnL_emuA = calc_neglnL2(2903, 0, 22, 1835, 1, 2, 0, 632, 97, mu_signal_only, sf_min);

    lnL_emuB = calc_neglnL2(5736, 0, 1, 5654, 0, 0, 0, 13, 5, mu_signal_only, sf_min);

    lnL_Zjets0 = calc_neglnL2(6202, 723, 442, 289, 3254, 0, 0, 124, 50, mu_signal_only, sf_min);

    lnL_Zjets1 = calc_neglnL2(8175, 448, 593, 425, 0, 4940, 0, 87, 78, mu_signal_only, sf_min);

    lnL_Zjets2 = calc_neglnL2(6710, 336, 499, 343, 0, 0, 5151, 27, 77, mu_signal_only, sf_min);

    lnL_total = lnL_SR + lnL_3l + lnL_emuA + lnL_emuB + lnL_Zjets0 + lnL_Zjets1 + lnL_Zjets2;

    if (lnL_total < lnL_min or i == 0)
    {
      lnL_min = lnL_total;
      min_i = i;

      sf_min.signal = sf.signal;
    }
    i++;
  }

  

  cout << "   ----------------------------------------------------------------------------" << endl << endl;
  cout << "   SIGNAL ONLY FIT:        " << endl;
  cout << "   _______________________  " << endl << endl;
  cout << "   lnL =          " << lnL_min << endl << endl;
  cout << "   mu_signal =    " << sf_min.signal << endl << endl;
  cout << "   mu_WZ =        " << sf_min.WZ << endl << endl;
  cout << "   mu_WW =        " << sf_min.WW << endl << endl;
  cout << "   mu_top =       " << sf_min.top << endl << endl;
  cout << "   mu_Zjets0 =    " << sf_min.Zjets0 << endl << endl;
  cout << "   mu_Zjets1 =    " << sf_min.Zjets1 << endl << endl;
  cout << "   mu_Zjets2 =    " << sf_min.Zjets2 << endl << endl;
  cout << "   ----------------------------------------------------------------------------" << endl << endl;


  // Timer stop
  auto end = chrono::high_resolution_clock::now();

  chrono::duration<float> duration = end - start;

  cout << endl
       << "   Script executed in: " << int(duration.count() / 60.0) << " minutes"
       << " and " << int((duration.count() / 60.0 - int(duration.count() / 60.0)) * 60) << " s" << endl
       << endl;

  // For the log file
  std::cout.rdbuf(oldBuffer); // Restore the original stream buffer for std::cout

  logFile.close(); // Close the log file

  return;

  }
