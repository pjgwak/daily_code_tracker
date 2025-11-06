// pt-dependent ctau cuts
const double ctaucut_a_mid_pp = 0.010;
const double ctaucut_a_fwd_pp = 0.013;
const double ctaucut_b_mid_pp = 0.25;
const double ctaucut_b_fwd_pp = 0.29;

const bool datamc = false;

const bool integrate = true;
// if integrate = false, then plot the ljpsi distribution
// if integrate = true, then plot the integrated distribution (cumulative distribution)

void integrateHist(TH1F* hist);

void plot_ljpsi(const char* fdata, const char* fmc) {
   TFile *tfd = TFile::Open(fdata);
   TFile *tfm = TFile::Open(fmc);

   TCanvas *c1 = new TCanvas();
   c1->cd();
   TPad *p1 = new TPad("p1","",0,0,1,0.33);
   p1->SetBottomMargin(2.*gStyle->GetPadBottomMargin());
   p1->SetTopMargin(0);
   p1->SetFrameBorderMode(0);
   p1->SetBorderMode(0);
   p1->SetBorderSize(0);
   p1->SetGridy();
   if (datamc) p1->Draw();
   TPad *p2 = new TPad("p2","",0,0.33,1,1);
   // p2->SetTopMargin(gStyle->GetPadTopMargin()/(1.-0.2));
   p2->SetBottomMargin(0);
   p2->SetFrameBorderMode(0);
   p2->SetBorderMode(0);
   p2->SetBorderSize(0);
   // p2->SetLogy();
   if (datamc) p2->Draw();
   p1->cd();
   if (!datamc) c1->cd();

   TF1 *dg = new TF1("f1","[0]*([1]*TMath::Gaus(x,[2],[3])+(1-[1])*TMath::Gaus(x,[4],[3]*[5]))",-0.2,0.04);
   dg->SetParLimits(0,0,1e15);
   dg->SetParLimits(1,0.5,1);
   dg->SetParLimits(2,-0.05,0.05);
   dg->SetParLimits(3,0,0.1);
   dg->SetParLimits(4,-0.05,0.05);
   dg->SetParLimits(5,0,100);
   dg->SetParameters(1000,0.8,1e-3,0.03,1e-4,2.);

   TIter next(tfd->GetListOfKeys()); TObject *obj;
   while ((obj = next())) {
      obj = ((TKey*) obj)->ReadObj();
      if (TString(obj->ClassName()) == "TH1F") {
         TH1F *hdata = (TH1F*) tfd->Get(obj->GetName());
         TH1F *hmc = (TH1F*) tfm->Get(obj->GetName());
         // hdata->Rebin(6);
         // hmc->Rebin(6);

         if (!integrate) {
            double int_data = hdata->Integral(hdata->FindBin(-0.015), hdata->FindBin(0.015));
            double int_mc = hmc->Integral(hmc->FindBin(-0.015), hmc->FindBin(0.015));
            hmc->Scale(int_data/int_mc);
         } else {
            double int_mc_all = hmc->Integral(0, hmc->GetNbinsX()+1);
            double int_data_all = hdata->Integral(0, hdata->GetNbinsX()+1);
            hmc->Scale(1./int_mc_all);
            hdata->Scale(1./int_data_all);
            double int_data = hdata->Integral(0, hdata->FindBin(0.));
            double int_mc = hmc->Integral(0, hmc->FindBin(0.));
            hdata->Scale(int_mc/int_data); 
            integrateHist(hdata);
            integrateHist(hmc);
         }
         hdata->GetXaxis()->SetTitle("#font[12]{l}_{J/#psi} [mm]");
         hdata->GetYaxis()->SetTitle("Entries");
         if (integrate) hdata->GetYaxis()->SetTitle("c.d.f.");

         // read the bin from the name
         TString thname(hdata->GetName());
         bool isfwd = (thname.Index("fwd") != kNPOS);
         int pos1 = thname.Index("pt")+2;
         int pos2 = thname.Index("-")-pos1;
         TString tsubstr(thname(pos1,pos2));
         float ptmin = atof(tsubstr.Data());
         pos1 = thname.Index("-")+1;
         pos2 = thname.Length();
         tsubstr = TString(thname(pos1,pos2));
         float ptmax = atof(tsubstr.Data());

         if (datamc) p2->cd();
         hdata->Draw();
         hmc->SetLineColor(kRed);
         hmc->SetMarkerColor(kRed);
         hmc->Draw("same");
         gPad->Update();

         // display the position of the ctau cut
         double a=1,b=1;
         if (!isfwd){a=ctaucut_a_mid_pp; b=ctaucut_b_mid_pp;} 
         else {a=ctaucut_a_fwd_pp; b=ctaucut_b_fwd_pp;}
         double cutmin = a + b / ptmin;
         double cutmax = a + b / ptmax;
         double xmin, xmax, ymin, ymax;
         TBox *tb = new TBox(cutmin,gPad->GetUymin(),cutmax,gPad->GetUymax());
         tb->SetFillColor(kBlack);
         tb->SetFillStyle(3003);
         tb->Draw();

         if (integrate) {
            TLine *tl = new TLine(hdata->GetXaxis()->GetXmin(),0.9,hdata->GetXaxis()->GetXmax(),0.9);
            tl->SetLineStyle(3);
            tl->SetLineColor(kBlack);
            tl->SetLineWidth(5);
            tl->Draw();
         }

         if (datamc) {
            p1->cd();
            TH1F *hratio = (TH1F*) hdata->Clone("hratio");
            hratio->Divide(hmc);
            hratio->GetYaxis()->SetTitle("data/MC");
            hratio->GetYaxis()->SetRangeUser(0.5,1.5);
            hratio->GetXaxis()->SetLabelSize(2.*hratio->GetXaxis()->GetLabelSize());
            hratio->GetXaxis()->SetTitleSize(2.*hratio->GetXaxis()->GetTitleSize());
            hratio->GetYaxis()->SetLabelSize(2.*hratio->GetYaxis()->GetLabelSize());
            hratio->GetYaxis()->SetTitleSize(2.*hratio->GetYaxis()->GetTitleSize());
            hratio->Draw();

            p2->cd();
         }

         // dg->SetParameters(1000,0.8,1e-3,0.03,1e-4,2.);
         // dg->SetParameter(0,int_data);
         // dg->SetParLimits(0,0,int_data*100.);
         // hdata->Fit(dg, "LERQ");
         // double sigma_data = dg->GetParameter(3);
         // double dsigma_data = dg->GetParError(3);
         // dg->SetRange(-5.*sigma_data,1.5*sigma_data);
         // hdata->Fit(dg, "LERQ");
         // sigma_data = dg->GetParameter(3);
         // dsigma_data = dg->GetParError(3);
         // dg->SetRange(-5.*sigma_data,1.5*sigma_data);
         // hdata->Fit(dg, "LERQ");
         // sigma_data = dg->GetParameter(3);
         // dsigma_data = dg->GetParError(3);
         // hmc->Fit(dg, "LERQ");
         // double sigma_mc = dg->GetParameter(3);
         // double dsigma_mc = dg->GetParError(3);

         double lx1=0.61, lx2=0.9, ly1=0.6, ly2=0.88;
         if (integrate) {ly1-=0.4; ly2-=0.4;}
         TLegend *tleg = new TLegend(lx1,ly1,lx2,ly2);
         tleg->SetBorderSize(0);
         TString header("#splitline{");
         header += isfwd ? "1.6<|y|<2.4" : "|y|<1.6";
         header += Form("}{%.1f<pt<%.1f GeV/c}",ptmin,ptmax);
         tleg->SetHeader(header);
         // tleg->AddEntry(hdata,Form("pp data (#sigma = %.3f +/- %.3f)",sigma_data,dsigma_data),"lp");
         // tleg->AddEntry(hmc,Form("pp prompt J/#psi mc (#sigma = %.3f +/- %.3f)",sigma_mc,dsigma_mc),"lp");
         tleg->AddEntry(hdata,"data","lp");
         tleg->AddEntry(hmc,"prompt J/#psi mc","lp");
         tleg->Draw();
         c1->SaveAs(Form("%s.pdf",hdata->GetName()));
         c1->SaveAs(Form("%s.png",hdata->GetName()));
      }
   }
}

void integrateHist(TH1F* hist) {
   int nentries = hist->GetNbinsX();
   TH1F *historig = (TH1F*) hist->Clone("historig");
   for (int i=1; i<nentries+2; i++) {
      double val, err;
      val = historig->IntegralAndError(0,i,err);
      hist->SetBinContent(i,val);
      hist->SetBinError(i,err);
   }
   delete historig;
}
