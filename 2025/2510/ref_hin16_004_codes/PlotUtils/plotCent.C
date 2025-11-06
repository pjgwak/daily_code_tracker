#include "TFile.h"
#include "TH1.h"
#include "TArrow.h"
#include "TLine.h"
#include "TLatex.h"
#include "TCanvas.h"
#include <iostream>
#include "../Fitter/Macros/Utilities/EVENTUTILS.h"

using namespace std;
using namespace HI;

float getmean(TH1F* hist, double xmin, double xmax);
float getmeanB(TH1F* hist, int imin, int imax);

void plotCent(const char* filename="histos.root") {
   TFile *f = TFile::Open(filename);
   if (!f || !f->IsOpen()) return;

   TH1F *hCentrality_mid = (TH1F*) f->Get("hCentrality_mid");
   TH1F *hCentrality_fwd = (TH1F*) f->Get("hCentrality_fwd");
   TH1F *hNpart_mid = (TH1F*) f->Get("hdNpart_mid");
   TH1F *hNpart_fwd = (TH1F*) f->Get("hdNpart_fwd");

   TCanvas *c1 = new TCanvas();
   hCentrality_mid->GetXaxis()->SetTitle("Centrality percentile (%)");
   hCentrality_mid->GetYaxis()->SetTitle("Entries");
   hCentrality_mid->Draw();
   double ymin = hCentrality_mid->GetMinimum();
   double ymax = 1.*hCentrality_mid->GetMaximum();
   cout << ymin << " " << ymax << endl;
   // bin limits
   TLine *tl = new TLine(); 
   tl->SetLineStyle(2); tl->SetLineWidth(2);
   tl->DrawLine(10,ymin,10,ymax);
   tl->DrawLine(20,ymin,20,ymax);
   tl->DrawLine(30,ymin,30,ymax);
   tl->DrawLine(40,ymin,40,ymax);
   tl->DrawLine(50,ymin,50,ymax);
   // average centrality in each bin
   TArrow *ta = new TArrow(); ta->SetLineStyle(1); ta->SetLineWidth(2);
   float x = getmean(hCentrality_mid,0,10);
   ta->DrawArrow(x,0.3*ymax,x,ymin+0.0*(ymax-ymin));
   x = getmean(hCentrality_mid,10,20);
   ta->DrawArrow(x,0.3*ymax,x,ymin+0.0*(ymax-ymin));
   x = getmean(hCentrality_mid,20,30);
   ta->DrawArrow(x,0.3*ymax,x,ymin+0.0*(ymax-ymin));
   x = getmean(hCentrality_mid,30,40);
   ta->DrawArrow(x,0.3*ymax,x,ymin+0.0*(ymax-ymin));
   x = getmean(hCentrality_mid,40,50);
   ta->DrawArrow(x,0.3*ymax,x,ymin+0.0*(ymax-ymin));
   x = getmean(hCentrality_mid,50,100);
   ta->DrawArrow(x,0.3*ymax,x,ymin+0.0*(ymax-ymin));

   TLatex *tlx = new TLatex();
   tlx->DrawLatex(60,0.8*ymax,"PbPb data, |y|<1.6");
   c1->SaveAs("centrality_mid.pdf");

   // FWD
   hCentrality_fwd->GetXaxis()->SetTitle("Centrality percentile (%)");
   hCentrality_fwd->GetYaxis()->SetTitle("Entries");
   hCentrality_fwd->Draw();
   ymin = hCentrality_fwd->GetMinimum();
   ymax = 1.*hCentrality_fwd->GetMaximum();
   cout << ymin << " " << ymax << endl;
   // bin limits
   tl->SetLineStyle(2); tl->SetLineWidth(2);
   tl->DrawLine(20,ymin,20,ymax);
   tl->DrawLine(40,ymin,40,ymax);
   // average centrality in each bin
   x = getmean(hCentrality_fwd,0,20);
   ta->DrawArrow(x,0.3*ymax,x,ymin+0.0*(ymax-ymin));
   x = getmean(hCentrality_fwd,20,40);
   ta->DrawArrow(x,0.3*ymax,x,ymin+0.0*(ymax-ymin));
   x = getmean(hCentrality_fwd,40,100);
   ta->DrawArrow(x,0.3*ymax,x,ymin+0.0*(ymax-ymin));

   tlx->DrawLatex(50,0.8*ymax,"PbPb data, 1.6<|y|<2.4");
   c1->SaveAs("centrality_fwd.pdf");

   // NPART
   // MID
   hNpart_mid->GetXaxis()->SetTitle("N_{part}");
   hNpart_mid->GetYaxis()->SetTitle("Entries");
   hNpart_mid->Draw();
   ymin = hNpart_mid->GetMinimum();
   ymax = 1.*hNpart_mid->GetMaximum();
   cout << ymin << " " << ymax << endl;
   // bin limits
   tl->SetLineStyle(2); tl->SetLineWidth(2);
   tl->DrawLine(findNpart(20),ymin,findNpart(20),ymax);
   tl->DrawLine(findNpart(40),ymin,findNpart(40),ymax);
   tl->DrawLine(findNpart(60),ymin,findNpart(60),ymax);
   tl->DrawLine(findNpart(80),ymin,findNpart(80),ymax);
   tl->DrawLine(findNpart(100),ymin,findNpart(100),ymax);
   // average centrality in each bin
   x = getmeanB(hNpart_mid,181,201);
   ta->DrawArrow(x,0.3*ymax,x,ymin+0.0*(ymax-ymin));
   x = getmeanB(hNpart_mid,161,181);
   ta->DrawArrow(x,0.3*ymax,x,ymin+0.0*(ymax-ymin));
   x = getmeanB(hNpart_mid,141,161);
   ta->DrawArrow(x,0.3*ymax,x,ymin+0.0*(ymax-ymin));
   x = getmeanB(hNpart_mid,121,141);
   ta->DrawArrow(x,0.3*ymax,x,ymin+0.0*(ymax-ymin));
   x = getmeanB(hNpart_mid,101,121);
   ta->DrawArrow(x,0.3*ymax,x,ymin+0.0*(ymax-ymin));
   x = getmeanB(hNpart_mid,1,101);
   ta->DrawArrow(x,0.3*ymax,x,ymin+0.0*(ymax-ymin));

   tlx->DrawLatex(50,0.8*ymax,"PbPb data, |y|<1.6");
   c1->SaveAs("Npart_mid.pdf");

   // FWD
   hNpart_fwd->GetXaxis()->SetTitle("N_{part}");
   hNpart_fwd->GetYaxis()->SetTitle("Entries");
   hNpart_fwd->Draw();
   ymin = hNpart_fwd->GetMinimum();
   ymax = 1.*hNpart_fwd->GetMaximum();
   cout << ymin << " " << ymax << endl;
   // bin limits
   tl->SetLineStyle(2); tl->SetLineWidth(2);
   tl->DrawLine(findNpart(40),ymin,findNpart(40),ymax);
   tl->DrawLine(findNpart(80),ymin,findNpart(80),ymax);
   // average centrality in each bin
   x = getmeanB(hNpart_fwd,161,201);
   ta->DrawArrow(x,0.3*ymax,x,ymin+0.0*(ymax-ymin));
   x = getmeanB(hNpart_fwd,121,161);
   ta->DrawArrow(x,0.3*ymax,x,ymin+0.0*(ymax-ymin));
   x = getmeanB(hNpart_fwd,1,121);
   ta->DrawArrow(x,0.3*ymax,x,ymin+0.0*(ymax-ymin));

   tlx->DrawLatex(50,0.8*ymax,"PbPb data, 1.6<|y|<2.4");
   c1->SaveAs("Npart_fwd.pdf");
}

float getmean(TH1F* hist, double xmin, double xmax) {
   int imin = hist->FindBin(xmin);
   int imax = hist->FindBin(xmax);
   cout << imin << " " << imax << endl;
   return getmeanB(hist, imin, imax);
}

float getmeanB(TH1F* hist, int imin, int imax) {
   double tot1=0, tot2=0;
   for (int i=imin; i<imax; i++) {
      tot1 += hist->GetBinContent(i)*hist->GetBinCenter(i);
      tot2 += hist->GetBinContent(i);
   }
   return (tot2!=0) ? tot1/tot2 : 0;
}
