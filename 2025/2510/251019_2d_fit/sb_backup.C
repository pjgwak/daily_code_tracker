void drawCtauSBPlots(RooWorkspace *ws, RooDataHist *redDataSB, RooDataHist *binDataCtErrSB, RooFitResult *fitCt_Bkg, float lmin, float lmax, double *UnNormChi2_side_t, int *nFitParam_side_t, int *nFullBinsPull_side_t, int *Dof_side_t, double *Chi2_side_t)
{

  char reduceDS[1000];
  string titlestr;
  double unNormChi2;
  int dof;

  RooBinning rb(ws->var("ctau3D")->getBinning().numBins(), ws->var("ctau3D")->getBinning().array());

  TLegend leg11(0.64, 0.65, 0.9, 0.74, NULL, "brNDC");
  leg11.SetFillStyle(0);
  leg11.SetBorderSize(0);
  leg11.SetShadowColor(0);
  leg11.SetMargin(0.2);
  leg11.SetTextSize(0.040);
  leg11.SetTextFont(42);
  // leg11.AddEntry(gfake1, "sideband data", "p");
  // leg11.AddEntry(&hfake11, "background", "l");

  RooPlot *tframe1 = ws->var("ctau3D")->frame(Bins(100));
  double avgBinWidth = rb.averageBinWidth();
  tframe1->GetYaxis()->SetTitle(Form("Counts / (%.0f #mum)", avgBinWidth * 1000));
  tframe1->GetYaxis()->CenterTitle(1);
  redDataSB->plotOn(tframe1, DataError(RooAbsData::SumW2));

  // original
  // ws->pdf("CtBkgTot_PEE")->plotOn(tframe1, ProjWData(RooArgList(*(ws->var("ctau3DErr"))), *binDataCtErrSB, kTRUE), NumCPU(16), Normalization(1, RooAbsReal::NumEvent), LineStyle(7));

  // ----------------------------------------------
  // pjgwak
  ws->pdf("CtBkgPos")->plotOn(tframe1, ProjWData(RooArgSet(*(ws->var("ctau3DErr"))), *binDataCtErrSB), LineColor(kBlue), LineStyle(kDashed), Normalization(binDataCtErrSB->sumEntries() * (1 - ws->var("fracCtBkg3")->getVal()) * ws->var("fracCtBkg2")->getVal() * ws->var("fracCtBkg1")->getVal(), RooAbsReal::NumEvent));
  RooCurve *curvePos = (RooCurve *)tframe1->getObject(tframe1->numItems() - 1);

  // --- CtBkgNeg ---
  ws->pdf("CtBkgNeg")->plotOn(tframe1, ProjWData(RooArgSet(*(ws->var("ctau3DErr"))), *binDataCtErrSB), LineColor(kGreen), LineStyle(kDashed), Normalization(binDataCtErrSB->sumEntries() * (1 - ws->var("fracCtBkg3")->getVal()) * ws->var("fracCtBkg2")->getVal() * (1 - ws->var("fracCtBkg1")->getVal()), RooAbsReal::NumEvent));
  RooCurve *curveNeg = (RooCurve *)tframe1->getObject(tframe1->numItems() - 1);

  // --- CtBkgDbl ---
  ws->pdf("CtBkgDbl")->plotOn(tframe1, ProjWData(RooArgSet(*(ws->var("ctau3DErr"))), *binDataCtErrSB), LineColor(kOrange), LineStyle(kDashed), Normalization(binDataCtErrSB->sumEntries() * (1 - ws->var("fracCtBkg3")->getVal()) * (1 - ws->var("fracCtBkg2")->getVal()), RooAbsReal::NumEvent));
  RooCurve *curveDbl = (RooCurve *)tframe1->getObject(tframe1->numItems() - 1);

  // --- CtPRRes (Resolution) ---
  ws->pdf("CtPRRes")->plotOn(tframe1, ProjWData(RooArgSet(*(ws->var("ctau3DErr"))), *binDataCtErrSB), LineColor(kMagenta), LineStyle(kDashed), Normalization(binDataCtErrSB->sumEntries() * ws->var("fracCtBkg3")->getVal(), RooAbsReal::NumEvent));
  RooCurve *curveRes = (RooCurve *)tframe1->getObject(tframe1->numItems() - 1);

  auto sum12 = new RooCurve("sum12", "", *curvePos, *curveNeg);
  auto sum123 = new RooCurve("sum123", "", *sum12, *curveDbl);
  auto sumAll = new RooCurve("sumAll", "", *sum123, *curveRes);

  sumAll->SetLineColor(kBlack);
  sumAll->SetLineStyle(kSolid);
  sumAll->SetLineWidth(3);
  sumAll->SetMarkerStyle(1);
  sumAll->SetDrawOption("L");

  tframe1->addObject(sumAll);
  // ----------------------------------------------
  

  TCanvas *c3 = new TCanvas("c3", "The Canvas", 200, 10, 600, 750);
  c3->cd();
  TPad *pad1 = new TPad("pad1", "This is pad1", 0.00, 0.3, 1.0, 1.0);
  pad1->SetLeftMargin(0.14);
  pad1->SetRightMargin(0.03);
  pad1->SetTopMargin(0.075);
  pad1->SetBottomMargin(0);
  pad1->Draw();
  TPad *pad2 = new TPad("pad2", "This is pad2", 0.00, 0.00, 1.0, 0.3);
  pad2->SetLeftMargin(0.14);
  pad2->SetRightMargin(0.03);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.30);
  pad2->Draw();

  pad1->cd();
  tframe1->Draw();
  // lumiTextOffset = 0.45;
  // CMS_lumi(c3, opt.isPA, iPos);
  TLatex *t = new TLatex();
  t->SetNDC();
  t->SetTextAlign(32);
  t->SetTextSize(0.040);
  // t->DrawLatex(0.92, 0.89, opt.rapString);
  // t->DrawLatex(0.92, 0.83, opt.ptString);
  // if (opt.EventActivity == 1)
  //   t->DrawLatex(0.92, 0.78, opt.ntrkString);
  // else if (opt.EventActivity == 2)
  //   t->DrawLatex(0.92, 0.78, opt.etString);
  leg11.Draw("same");

  //// *** pull
  TH1 *hdatasb = redDataSB->createHistogram("hdatasb", *ws->var("ctau3D"), Binning(rb));
  RooHist *hpullsb = tframe1->pullHist();
  hpullsb->SetName("hpullSB");
  int nFitPar = fitCt_Bkg->floatParsFinal().getSize();
  double chi2 = 0;
  double *ypullssb = hpullsb->GetY();
  unsigned int nBins = ws->var("ctau3D")->getBinning().numBins();
  unsigned int nFullBins = 0;
  for (unsigned int i = 0; i < nBins; i++)
  {
    if (hdatasb->GetBinContent(i + 1) == 0)
      continue;
    chi2 += ypullssb[i] * ypullssb[i];
    nFullBins++;
  }
  unNormChi2 = chi2;
  *UnNormChi2_side_t = chi2;
  dof = nFullBins - nFitPar;
  chi2 /= (nFullBins - nFitPar);
  int nDOF = ws->var("ctau3D")->getBinning().numBins() - nFitPar;
  *nFitParam_side_t = nFitPar;
  *nFullBinsPull_side_t = nFullBins;
  *Dof_side_t = dof;
  *Chi2_side_t = chi2;

  RooPlot *tframepull = ws->var("ctau3D")->frame(Title("Pull"));
  tframepull->GetYaxis()->SetTitle("Pull");
  tframepull->GetYaxis()->CenterTitle(1);
  tframepull->SetLabelSize(0.04 * 2.5, "XYZ");
  tframepull->SetTitleSize(0.048 * 2.5, "XYZ");
  tframepull->SetTitleOffset(0.47, "Y");
  tframepull->addPlotable(hpullsb, "PX");
  double tframemax = 0;
  if (tframepull->GetMinimum() * -1 > tframepull->GetMaximum())
    tframemax = tframepull->GetMinimum() * -1;
  else
    tframemax = tframepull->GetMaximum();
  tframepull->SetMaximum(tframemax);
  tframepull->SetMinimum(-1 * tframemax);
  tframepull->GetXaxis()->SetTitle("#font[12]{l}_{J/#psi} (mm)");
  tframepull->GetXaxis()->CenterTitle(1);

  pad2->cd();
  tframepull->Draw();
  TLine *line1 = new TLine(-lmin, 0, lmax, 0.);
  line1->SetLineStyle(7);
  line1->Draw();

  TLatex *t2 = new TLatex();
  t2->SetNDC();
  t2->SetTextAlign(22);
  t2->SetTextSize(0.035 * 3);
  sprintf(reduceDS, "#chi^{2}/dof = %.2f/%d", unNormChi2, dof);
  t2->DrawLatex(0.76, 0.90, reduceDS);

  titlestr = "_CtSB_Lin.pdf";
  c3->SaveAs(titlestr.c_str());

  pad1->SetLogy(1);
  double originalmax = tframe1->GetMaximum();
  tframe1->SetMaximum(originalmax * 10);
  tframe1->SetMinimum(0.5);
  titlestr = "_CtSB_Log.pdf";
  c3->SaveAs(titlestr.c_str());

  delete pad1;
  delete pad2;
  delete c3;
}