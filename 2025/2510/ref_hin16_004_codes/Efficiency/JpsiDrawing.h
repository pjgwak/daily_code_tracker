#include <TGraph.h>
#include <TH1.h>
#include <TLegend.h>

//int colorArr[] = {kBlack, kRed+1, kGreen+2, kBlue, kViolet-4, kMagenta, kMagenta+2};
int colorArr[] = {kBlack, kRed+1, kOrange+7, kGreen+3, kBlue+1, kViolet-4, kMagenta, kMagenta+2};
//int colorArr[] = {kBlack, kRed+1, kOrange+7, kSpring+4, kGreen+3, kAzure+1, kBlue+2, kViolet+5, kViolet-4, kMagenta, kMagenta+2};
//int colorArr[] = {kRed+1, kSpring+4, kAzure+1, kBlue+2, kMagenta, kMagenta+2};
//int colorArr[] = {kRed+1, kSpring+4, kBlue+2, kMagenta, kMagenta+2};
int markerArr[] = {kOpenCircle, kOpenSquare, kOpenStar, kOpenTriangleUp, kOpenDiamond, kOpenCross, kFullCircle, kFullSquare, kFullStar};
int ncolor = sizeof(colorArr)/sizeof(int);
int nmarker = sizeof(markerArr)/sizeof(int);

void SetHistStyleDefault(TH1 *h, int i, int j) {
//  h->SetMarkerSize(1.800);
//  if (j == 2 || j ==4) h->SetMarkerSize(1.800);
//  if (j == 5) h->SetMarkerSize(1.500);

  if (ncolor>i) {
    h->SetMarkerColor(colorArr[i]);
    h->SetLineColor(colorArr[i]);
  } else {
    h->SetMarkerColor(colorArr[i%ncolor]);
    h->SetLineColor(colorArr[i%ncolor]);
  }
  if (nmarker>j) {
    h->SetMarkerStyle(markerArr[j]);
  } else {
    h->SetMarkerStyle(markerArr[j%nmarker]);
  }

  h->GetXaxis()->SetTitleSize(0.048);
  h->GetYaxis()->SetTitleSize(0.048);
  h->GetXaxis()->SetLabelSize(0.048);
  h->GetYaxis()->SetLabelSize(0.048);
//  h->GetYaxis()->SetTitleOffset(1.9);
}

void SetHistStyleDefault(TGraph *h, int i, int j) {
  h->SetMarkerSize(1.200);
//  if (j == 2 || j ==4) h->SetMarkerSize(1.800);
//  if (j == 5) h->SetMarkerSize(1.500);

  if (ncolor>i) {
    h->SetMarkerColor(colorArr[i]);
    h->SetLineColor(colorArr[i]);
  } else {
    h->SetMarkerColor(colorArr[i%ncolor]);
    h->SetLineColor(colorArr[i%ncolor]);
  }
  if (nmarker>j) {
    h->SetMarkerStyle(markerArr[j]);
  } else {
    h->SetMarkerStyle(markerArr[j%nmarker]);
  }

  h->GetXaxis()->SetTitleSize(0.048);
  h->GetYaxis()->SetTitleSize(0.048);
  h->GetXaxis()->SetLabelSize(0.048);
  h->GetYaxis()->SetLabelSize(0.048);
}

void SetHistStyle(TGraph *h, int i, int j, double rmin, double rmax){
  h->SetMinimum(rmin);
  h->SetMaximum(rmax);
  SetHistStyleDefault(h, i, j);
}

void SetHistStyle(TH1 *h, int i, int j, double rmin, double rmax){
  h->GetYaxis()->SetRangeUser(rmin,rmax);
  SetHistStyleDefault(h, i, j);
}

void SetLegendStyle(TLegend* l) {
  l->SetFillColor(0);
  l->SetFillStyle(4000);
  l->SetBorderSize(0);
  l->SetMargin(0.15);
}

