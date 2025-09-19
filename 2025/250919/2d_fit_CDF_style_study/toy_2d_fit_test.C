#include <TStopwatch.h>
#include <RooRealVar.h>

using namespace RooFit;

RooDataSet *makeFakeDataXY()
{
  TRandom3 trnd{};

  RooRealVar dt("dt", "dt", -10, 10);
  RooRealVar dterr("dterr", "dterr", 0.0001, 5);
  RooArgSet coord(dt, dterr);

  RooDataSet *d = new RooDataSet("d", "d", RooArgSet(dt, dterr));

  for (int i = 0; i < 10000; i++)
  {
    double tmpy = trnd.Gaus(0, 10);
    double tmpx = trnd.Gaus(0.5 * tmpy, 1);
    if (fabs(tmpy) < 10 && fabs(tmpx) < 10)
    {
      dt.setVal(tmpx);
      dterr.setVal(tmpy);
      d->add(coord);
    }
  }

  return d;
}

void toy_2d_fit_test()
{
  cout << "=== start toy_2d_fit_test()";
  TStopwatch t;
  t.Start();

  

  // variables
  RooRealVar dt("dt", "dt", -10, 10);
  RooRealVar dterr("dterr", "dterr", 0.0001, 5);
  RooRealVar tau("tau", "", 2, -100, 10000);

  // build gauss res
  RooRealVar sigma("sigma", "", 1, 0.001, 100);
  RooGaussModel gm("gm", "", dt, RooConst(0), sigma, dterr); // mean = 0

  // construct decay(t, tau) (x) gauss1(t, 0, sigma*dterr)
  RooDecay decay_gm("decay_gm", "", dt, tau, gm, RooDecay::DoubleSided);

  // Obtain fake external experimental dataset with values for x and y - to build toy dterrData 
  RooDataSet *expDataXY = makeFakeDataXY();

  // dataset for dterr -> Should get from experiment result in practice
  std::unique_ptr<RooAbsData> expAbsDataY{expDataXY->reduce(dterr)};
  RooDataSet *dterrData = static_cast<RooDataSet *>(expAbsDataY.get());

  // Toy MC generation of decay(t|dterr)
  RooDataSet *data = decay_gm.generate(dt, ProtoData(*dterrData));

  // Fitting of decay(t|dterr)
  auto fitResult = decay_gm.fitTo(*data, Save(), ConditionalObservables(dterr), NumCPU(32), EvalBackend("legacy"));

  // plotting
  TCanvas c1("c1", "", 800, 800);
  c1.Divide(1, 2);
  TPad *pad1 = (TPad *)c1.cd(1);
  pad1->SetPad(0.0, 0.3, 1.0, 1.0);
  pad1->SetBottomMargin(0.02);

  TPad *pad2 = (TPad *)c1.cd(2);
  pad2->SetPad(0.0, 0.0, 1.0, 0.3);
  pad2->SetTopMargin(0.05);
  pad2->SetBottomMargin(0.25);

  // main plot
  pad1->cd();
  RooPlot *frame = dt.frame();
  data->plotOn(frame);
  decay_gm.plotOn(frame, ProjWData(*data));
  frame->Draw();

  // pull 
  pad2->cd();
  RooHist *hpull = frame->pullHist(); // data - fit / error
  RooPlot *frame_pull = dt.frame();
  frame_pull->addPlotable(hpull, "P");
  frame_pull->SetTitle("");
  frame_pull->GetYaxis()->SetTitle("Pull");
  frame_pull->GetYaxis()->SetNdivisions(505);
  frame_pull->GetYaxis()->SetTitleSize(0.08);
  frame_pull->GetYaxis()->SetLabelSize(0.08);
  frame_pull->GetXaxis()->SetTitleSize(0.1);
  frame_pull->GetXaxis()->SetLabelSize(0.08);
  frame_pull->Draw();

  // --- chi2/ndf ---
  int nFitParam = decay_gm.getParameters(*data)->selectByAttrib("Constant", kFALSE)->getSize();
  double chi2 = frame->chiSquare(nFitParam); // (chi2/ndf)
  std::cout << "Chi2/ndf = " << chi2 << std::endl;
  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.08);
  latex.DrawLatex(0.85, 0.85, Form("#chi^{2}/ndf = %.2f", chi2));

  

  c1.SaveAs("toy_2d_fit_test_withPull.png");

  cout << "=== finish toy_2d_fit_test()";
  t.Stop();
  printf("RealTime=%f seconds, CpuTime=%f seconds\n", t.RealTime(), t.CpuTime());

  fitResult->Print("V");
}