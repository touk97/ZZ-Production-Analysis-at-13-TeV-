

#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TLine.h>
#include <TAxis.h>
#include <TBox.h>
#include <TMarker.h>
#include <TLatex.h>
#include <TROOT.h>
#include "TStyle.h"

#include <iostream>
#include <fstream>

using namespace std;


void xsec_plot()
{
    gROOT->SetBatch(kTRUE);
    gStyle->SetTextFont(52);

    TCanvas *c1 = new TCanvas("c1", "Scaling Factors", 800, 600);


    //xsec    
    double xsec = 16.35;
    double xsec_error = 0.68;
   const Int_t n = 1;
   Double_t x[n]  = {16.35};
   Double_t y[n]  = {6};
   Double_t ex[n] = {0.68};
   Double_t ey[n] = {0};
   auto g1 = new TGraphErrors(n, x, y, ex, ey);
   // TGraphAsymmErrors *g1 = new TGraphAsymmErrors(1);
   //     cout << xsec - xsec_error << endl;
   // g1->SetPoint(1, xsec,1);
   // g1->SetPointError(1, xsec_error, xsec_error, 0, 0);

   g1->SetMarkerStyle(20);
   g1->SetMarkerSize(1.1);
   g1->SetMarkerColor(kBlack);
   g1->GetXaxis()->SetTitle("#sigma^{fid}_{ZZ} [fb]");
   g1->GetXaxis()->SetTitleFont(62);
   g1->GetXaxis()->SetTitleSize(0.04);

   // g1->GetXaxis()->SetNdivisions(505);
   g1->GetYaxis()->SetTitle("  ");
   g1->SetTitle("  ");
   g1->GetYaxis()->SetNdivisions(0);
   g1->SetName("Scaling Factors");

   float xmin = xsec - 4.8 * xsec_error;
   float xmax = xsec + 4.8 * xsec_error;

   float ymin = 4;
   float ymax = 8;

//    g1->GetXaxis()->SetRangeUser(xsec - 4 * xsec_error, xsec + 4 * xsec_error);
   g1->GetYaxis()->SetRangeUser(ymin, ymax);

   TAxis *axis = g1->GetXaxis();
   axis->SetLimits(xmin, xmax);
//    TAxis *axisy = g1->GetYaxis();
//    axisy->SetLimits(5, 8);
   g1->Draw("AP same");


   

   TBox *sigma2 = new TBox(xsec - 2 * xsec_error, ymin, xsec + 2 * xsec_error, ymax);
   // sigma2->SetFillColor(kBlack);
   sigma2->SetFillColorAlpha(kYellow, 0.6);
   sigma2->Draw("same");

   TBox *sigma1 = new TBox(xsec - xsec_error, ymin, xsec + xsec_error, ymax);
   sigma1->SetFillColorAlpha(kGreen, 0.6);
   sigma1->Draw("same");

   TMarker *marker = new TMarker(16.35, 6, 20); // 20 is the marker style

   TLine *hline = new TLine(xsec - xsec_error, 6, xsec + xsec_error, 6);
   hline->SetLineStyle(1);
   hline->SetLineWidth(2);

   TLine *v1line = new TLine(xsec - xsec_error, 6 - 0.06, xsec - xsec_error, 6 + 0.06);
   v1line->SetLineStyle(1);
   v1line->SetLineWidth(2);

   TLine *v2line = new TLine(xsec + xsec_error, 6 - 0.06, xsec + xsec_error, 6 + 0.06);
   v2line->SetLineStyle(1);
   v2line->SetLineWidth(2);

   hline->Draw("same");
   v1line->Draw("same");
   v2line->Draw("same");
   marker->Draw("same");
   // xsec
   const char *sigma_text = "#sigma^{fid}_{ZZ} = 16.35 #pm 0.68 fb";

   // TLatex *tex1 = new TLatex(1.03 * xmin, 0.96 * ymax, "ATLAS");
   TLatex *tex2 = new TLatex(1.02 * xmin, 0.96 * ymax, "ZZ \\rightarrow \\ell \\ell \\nu \\nu");
   TLatex *tex3 = new TLatex(1.02 * xmin, 0.84 * ymax, "SR Data Blinded");
   TLatex *tex4 = new TLatex(1.02 * xmin, 0.88 * ymax, sigma_text);
   TLatex *tex5 = new TLatex(1.02 * xmin, 0.92 * ymax, "#sqrt{s}  = 13 TeV, 140.1 fb^{-1}");

   // tex1->SetTextSize(0.04);
   // tex1->SetTextFont(72);
   // tex1->SetLineWidth(2);
   // tex1->Draw();

   tex2->SetTextSize(0.03);
   tex2->SetTextFont(82);
   tex2->SetLineWidth(2);
   tex2->Draw();

   tex3->SetTextSize(0.03);
   tex3->SetTextFont(72);
   tex3->SetLineWidth(2);
   tex3->Draw();

   tex4->SetTextSize(0.03);
   tex4->SetTextFont(72);
   tex4->SetLineWidth(2);
   tex4->Draw();
   
   tex5->SetTextSize(0.03);
   tex5->SetTextFont(72);
   tex5->SetLineWidth(2);
   tex5->Draw();

   c1->RedrawAxis();
   c1->SetGridx();
   c1->SetTicks();
   c1->Update();

   c1->SaveAs("./xsec_plot.png");
}
