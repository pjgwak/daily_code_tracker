const int nentries = 10000;

void dummydata(const char* input, const char* output, const char* suffix) {
   TFile *f = TFile::Open(input);
   RooWorkspace *ws = (RooWorkspace*) f->Get("workspace");
   // RooAbsPdf *thepdf = ws->pdf(Form("pdfMASS_Tot_%s",suffix));
   RooRealVar *invMass = ws->var("invMass");

   RooAbsPdf *sig1S = ws->pdf(Form("pdfMASS_Jpsi_%s",suffix));
   RooAbsPdf *sig2S = ws->pdf(Form("pdfMASS_Psi2S_%s",suffix));
   RooAbsPdf *pdf_bkgd = ws->pdf(Form("pdfMASS_Bkg_%s",suffix));
   RooRealVar *nsig1f = ws->var(Form("N_Jpsi_%s",suffix));
   RooRealVar *RFrac2Svs1S = ws->var(Form("RFrac2Svs1S_%s",suffix));
   RooRealVar *nbkgd = ws->var(Form("N_Bkg_%s",suffix));
   double nsig1f0 = nsig1f->getVal();
   double nbkgd0 = nbkgd->getVal();
   nsig1f->setVal(0.6*(nsig1f0+nbkgd0));
   nbkgd->setVal(0.4*(nsig1f0+nbkgd0));
   RFrac2Svs1S->setVal(0.2);
   RooFormulaVar *nsig2f = new RooFormulaVar(Form("N_Psi2S_%s",suffix),"@0*@1",RooArgList(*RFrac2Svs1S,*nsig1f));
   RooAddPdf *thepdf = new RooAddPdf (Form("pdfMASS_Tot_%s",suffix),"new total p.d.f.",
         RooArgList(*sig1S,*sig2S,*pdf_bkgd),
         RooArgList(*nsig1f,*nsig2f,*nbkgd));

   RooDataSet *data = thepdf->generate(RooArgSet(*invMass),nentries,RooFit::Extended(),RooFit::Name(Form("dOS_DATA_%s",suffix)),RooFit::Verbose(kTRUE));
   thepdf->fitTo(*data,RooFit::Extended(kTRUE),RooFit::NormRange("MassRange"),RooFit::NumCPU(3));
   RooPlot *frame = invMass->frame(RooFit::Bins(40));
   data->plotOn(frame);
   thepdf->plotOn(frame,RooFit::Range("MassRange"),RooFit::Normalization(nentries, RooAbsReal::NumEvent));
   TCanvas *c1 = new TCanvas();
   frame->Draw();
   c1->SaveAs("dummy.pdf");
   RooWorkspace *wout = new RooWorkspace("workspace");
   wout->import(*data);
   wout->import(*thepdf);
   wout->writeToFile(output);
}
