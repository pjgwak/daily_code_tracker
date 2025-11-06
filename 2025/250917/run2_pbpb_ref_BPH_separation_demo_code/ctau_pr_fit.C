#include "TCanvas.h"
#include "TPad.h"
#include "TLatex.h"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooAddPdf.h"
#include "RooGaussModel.h"
#include "RooDecay.h"
#include "RooAddModel.h"
#include "RooExtendPdf.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooExponential.h"
#include "RooHist.h"

using namespace RooFit;

void ctau_pr_fit()
{
  // ROOT::EnableImplicitMT(24);
  float ptLow = 25, ptHigh = 40;
  float yLow = 0, yHigh = 1.6;
  float massLow = 2.6, massHigh = 3.5;
  std::string comp = "CtauFull", region = "SR";
  int cLow = 0, cHigh = 180;

  TStopwatch t;
  t.Start();
  cout << "\n=== Start ctau_pr_fit() ===\n";

  // use rootlogon
  gROOT->Macro("/data/users/pjgwak/input_files/rootlogon.C");

  // make output folder
  gSystem->mkdir("figs", true);
  gSystem->mkdir("roots", true);

  // read input
  TFile *fInput = TFile::Open("/data/hwan/psi2S_RAA_PbPb2018/skimmedFiles/OniaRooDataSet_miniAOD_isMC0_JPsi_cent0_200_Effw0_Accw0_PtW0_TnP0_230721.root");
  if (!fInput || fInput->IsZombie())
  {
    cout << "Error: cannot open input file\n";
    return;
  }

  // read dataset
  RooDataSet *ds = dynamic_cast<RooDataSet *>(fInput->Get("dataset"));
  if (!ds)
  {
    cout << "Error: cannot find RooDataSet\n";
    return;
  }

  // === declare cuts ===
  // --- basic cuts ---
  // acceptance
  TString accCut = "( ((abs(eta1) <= 1.2) && (pt1 >=3.5)) || ((abs(eta2) <= 1.2) && (pt2 >=3.5)) || ((abs(eta1) > 1.2) && (abs(eta1) <= 2.1) && (pt1 >= 5.47-1.89*(abs(eta1)))) || ((abs(eta2) > 1.2)  && (abs(eta2) <= 2.1) && (pt2 >= 5.47-1.89*(abs(eta2)))) || ((abs(eta1) > 2.1) && (abs(eta1) <= 2.4) && (pt1 >= 1.5)) || ((abs(eta2) > 2.1)  && (abs(eta2) <= 2.4) && (pt2 >= 1.5)) )"; // 2018 acceptance cut

  // kinematics cuts
  //  - correct? (<= cBin <) -> maybe (< cBin <=) ??
  TString kineCut = Form( // tmp: no cbin
      "(pt >= %.3f && pt < %.3f && abs(y) >= %.3f && abs(y) < %.3f && "
      " mass >= %.3f && mass < %.3f)",
      ptLow, ptHigh, yLow, yHigh, massLow, massHigh);
  // TString kineCut = Form(
  //     "(pt >= %.3f && pt < %.3f && abs(y) >= %.3f && abs(y) < %.3f && "
  //     " mass >= %.3f && mass < %.3f && cBin >= %d && cBin < %d)",
  //     ptLow, ptHigh, yLow, yHigh, massLow, massHigh, cLow, cHigh);

  // opposite sign -> Must be FLowSkim
  TString osCut = "(recoQQsign == 0)";

  // --- region6 cuts ---
  const TString cutPR = "(abs(ctau3D) < 0.05)";
  const TString cutNP = "(ctau3D >= 0.10 && ctau3D <= 0.80)";
  // const TString cutCtauFull = "(ctau3D >= -0.1 && ctau3D <= 0.5)";
  const TString cutCtauFull = "(ctau3D >= -0.5 && ctau3D <= 2)";

  const TString cutSR = "(mass >= 3.00 && mass <= 3.20)";
  const TString cutLSB = "(mass >= 2.60 && mass <= 2.95)";
  const TString cutRSB = "(mass >= 3.21 && mass <= 3.50)";

  // --- combine cuts ---
  TString compCut;
  if (comp == "PR")
    compCut = cutPR;
  else if (comp == "NP")
    compCut = cutNP;
  else if (comp == "CtauFull")
    compCut = cutCtauFull;
  else
    compCut = "(1)";

  TString regionCut;
  if (region == "SR")
    regionCut = cutSR;
  else if (region == "LSB")
    regionCut = cutLSB;
  else if (region == "RSB")
    regionCut = cutRSB;
  else
    regionCut = "(1)";

  TString fullCut = Form("%s && %s && %s && %s && %s",
                         osCut.Data(),
                         accCut.Data(),
                         kineCut.Data(),
                         compCut.Data(),
                         regionCut.Data());

  // === new dataset with cuts ===
  RooDataSet *data = (RooDataSet *)ds->reduce(Cut(fullCut));
  if (!data || data->numEntries() == 0)
  {
    cout << "[ERROR] reduced dataset is empty. Cut = " << fullCut << "\n";
    return;
  }

  // print out object list
  cout << "\n--- Objects in reduced dataset ---\n";
  data->Print();

  // check weight
  // const bool hasWeight = (data->isWeighted() || data->weightVar() || data->get()->find("weight"));
  const bool hasWeight = false;

  // variables
  // === observables ===
  auto *ctau = dynamic_cast<RooRealVar *>(data->get()->find("ctau3D"));
  auto *ctRes = dynamic_cast<RooRealVar *>(data->get()->find("ctau3DRes")); // per-event res
  // ===== 범위: 데이터는 전체, NP는 ct>0만 기여/정규화 =====
  double ctMin = -0.5, ctMax = 4;
  // ctau->setMin(-0.5);
  // ctau->setMax(4);
  ctau->setRange(ctMin, ctMax);
  ctau->setRange("Full", ctMin, ctMax);
  ctau->setRange("Res",-0.05, 0.05);
  // ctau->setRange("pos", 0, 0.05);


  // === load MC ctau3D Res fit result ===
  // load input
  TFile fResIn(Form("roots/ctau_mc_res_pT%.1f_%.1f_y%.1f_%.1f.root", ptLow, ptHigh, yLow, yHigh), "open");
  RooFitResult *fitResultRes = (RooFitResult *)fResIn.Get("fitResult");
  if (!fitResultRes)
  {
    std::cerr << "fitResult not found in resolution fit!" << std::endl;
    return;
  }

  const RooArgList &params = fitResultRes->floatParsFinal();
  RooRealVar *sigma1Var = (RooRealVar *)params.find("sigma1");
  RooRealVar *r21Var = (RooRealVar *)params.find("r21");
  RooRealVar *r32Var = (RooRealVar *)params.find("r32");
  RooRealVar *fsumVar = (RooRealVar *)params.find("fsum"); // not fixed
  RooRealVar *fg12Var = (RooRealVar *)params.find("fg12"); // not fixed

  if (r21Var)
  {
    r21Var->setConstant(kTRUE);
    std::cout << "Fixed r21 = " << r21Var->getVal() << std::endl;
  }
  if (r32Var)
  {
    r32Var->setConstant(kTRUE);
    std::cout << "Fixed r32 = " << r32Var->getVal() << std::endl;
  }
  // if (fsumVar)
  // {
  //   fsumVar->setConstant(kTRUE);
  //   std::cout << "Fixed r32 = " << fsumVar->getVal() << std::endl;
  // }
  if (fg12Var)
  {
    fg12Var->setConstant(kTRUE);
    std::cout << "Fixed r32 = " << fg12Var->getVal() << std::endl;
  }

  RooRealVar mean("mean", "mean", 0); // fix to zero
  RooRealVar sigma1("sigma1", "sigma1", 0.5, 0.01, 0.9);
  RooFormulaVar sigma2("sigma2", "@0*@1", RooArgList(sigma1, *r21Var));
  RooFormulaVar sigma3("sigma3", "@0*@1", RooArgList(sigma2, *r32Var));

  RooGaussModel g1("g1", "res1", *ctau, mean, sigma1);
  RooGaussModel g2("g2", "res2", *ctau, mean, sigma2);
  RooGaussModel g3("g3", "res3", *ctau, mean, sigma3);

  RooRealVar fsum("fsum", "fractions of g1+g2", 0.4, 0.01, 1);
  RooRealVar fg12("fg12", "g1/(g1+g2)", 0.5, 0.0, 1.0);

  RooFormulaVar fg1("fg1", "@0*@1", RooArgList(*fsumVar, *fg12Var));
  RooFormulaVar fg2("fg2", "@0*(1-@1)", RooArgList(*fsumVar, *fg12Var));
  // RooRealVar fg1("fg1", "fg1", ->getVal(), 0.01, 1);
  // RooRealVar fg2("fg2", "fg2", fg2Var->getVal(), 0.01, 1);
  // fg2.setConstant(kTRUE);

  RooAddModel resModel("resModel", "3-Gaussian resolution",
                       RooArgList(g1, g2, g3),
                       RooArgList(fg1, fg2));

  // --- Prompt = delta ⊗ resolution ---
  RooAbsPdf *prompt = &resModel;

  // --- Nonprompt lifetime ---
  RooRealVar tau("tau", "nonprompt lifetime", 0.45, 0.3, 1);
  RooDecay decay("decay", "nonprompt", *ctau, tau, resModel, RooDecay::SingleSided);

  // === load ctau3D sideband fit result ===
  // load input
  TFile fCtauSideIn(Form("roots/ctau_NP_RSB_pT%.1f_%.1f_y%.1f_%.1f.root", ptLow, ptHigh, yLow, yHigh), "open");
  RooFitResult *fitCtauSide = (RooFitResult *)fCtauSideIn.Get("fitResult");
  if (!fitCtauSide)
  {
    std::cerr << "fitResult not found in ctau NPSB fit!" << std::endl;
    return;
  }

  const RooArgList &paramsCtauSide = fitCtauSide->floatParsFinal();
  RooRealVar *tau_bkg1Var = (RooRealVar *)paramsCtauSide.find("tau1");

  if (tau_bkg1Var)
  {
    tau_bkg1Var->setConstant(kTRUE);
    std::cout << "Fixed tau_bkg1 = " << tau_bkg1Var->getVal() << std::endl;
  }

  // RooRealVar tau_bkg1("tau_bkg1", "short lifetime", 4.2119e-01, 0.1, 1.0);
  RooDecay bkg("bkg", "NP short", *ctau, *tau_bkg1Var, resModel, RooDecay::SingleSided);
  // RooDecay decay_bkg2("decay_bkg2", "NP long", *ctau, tau_bkg2, resModel, RooDecay::SingleSided);
  // RooAddPdf bkg("bkg", "NP sum",
  //                 RooArgList(decay_bkg1, decay_bkg2), RooArgList(fNP_bkg));

  // === Yields ===
  // --- decide Nnp vs Nbkg ratio from NPS mass fit
  // load input
  TFile fMassNPS(Form("roots/mass_NP_MassFull_pT%.1f_%.1f_y%.1f_%.1f.root", ptLow, ptHigh, yLow, yHigh), "open");
  auto *f_mass_bkg = (TParameter<double> *)fMassNPS.Get("f_mass_bkg");
  RooRealVar NbkgRatio("NbkgRatio", "yield nonprompt", f_mass_bkg->GetVal());
  NbkgRatio.setConstant();

  RooRealVar Npr("Npr", "yield prompt", 3.7659e+04, 1, 1000000);
  RooRealVar Nnp("Nnp", "yield nonprompt", 2.8221e+04, 1, 100000);
  RooFormulaVar Nbkg("Nbkg", "@0 * @1", RooArgList(Nnp, NbkgRatio));
  // RooRealVar Nbkg("Nbkg", "yield bkg", 30000, 0, 1e7);

  RooExtendPdf extPR("extPR", "extended prompt", *prompt, Npr);
  RooExtendPdf extNP("extNP", "extended nonprompt", decay, Nnp);
  RooExtendPdf extBkg("extBkg", "extended background", bkg, Nbkg);

  RooAddPdf model("model", "PR+NP+Bkg",
                  RooArgList(extPR, extNP, extBkg));

  // --- Generate toy dataset ---
  // auto data = model.generate(t, 50000);

  // --- Fit ---
  // RooBinning b(-0.5, 4);
  // b.addUniform(5, -0.5, -0.05);
  // b.addUniform(80, -0.05, 1);
  // b.addUniform(5, 0.1, 4);

  auto dh = new RooDataHist("dh", "binned dataset", *ctau, *data);
  // auto dh = new RooDataHist("dh", "binned dataset", *ctau, *data);

  auto fitResult = model.fitTo(*data, Extended(kTRUE), PrintLevel(-1), Save(), Offset(true), NumCPU(32), EvalBackend("legacy"), Strategy(2)); //, BatchMode(true)

  // --- Plot main frame ---
  // --- divided canvas ---
  TCanvas c_ctau("c_ctau", "c_ctau", 800, 800);

  // --- main plot ---
  TPad *pad1 = new TPad("pad1", "pad1", 0.0, 0.25, 1.0, 1.0);
  pad1->SetBottomMargin(0.00001);
  pad1->SetLogy();
  pad1->Draw();
  pad1->cd();

  // --- subrange ---
  ctau->setRange("SR", -0.05, 0.05);
  RooArgSet obs(*ctau);
  auto fPrompt = prompt->createIntegral(obs, NormSet(obs), Range("SR"))->getVal();
  auto fNonprompt = decay.createIntegral(obs, NormSet(obs), Range("SR"))->getVal();
  auto fBkg = bkg.createIntegral(obs, NormSet(obs), Range("SR"))->getVal();

  // nEvt in subrange = (probability in subrange) × (yield in full range)
  double Npr_in = Npr.getVal() * fPrompt;
  double Nnp_in = Nnp.getVal() * fNonprompt;
  double Nbkg_in = Nbkg.getVal() * fBkg;
  double b_fraction = Nnp_in / (Npr_in + Nnp_in + Nbkg_in);

  // double ctMin = -0.1, ctMax = 0.8;
  RooPlot *f_ctau = ctau->frame(Range(ctMin, ctMax), Title("")); // Bins(80)
  data->plotOn(f_ctau, DataError(RooAbsData::SumW2), Name("data"));
  // data->plotOn(f_ctau, DataError(RooAbsData::SumW2), Name("data"), Binning(b));
  model.plotOn(f_ctau, Name("model"), Range("Full"), NormRange("Full"));
  model.plotOn(f_ctau, Components("extPR"), Name("extPR"), LineStyle(kDashed), LineColor(kRed));
  model.plotOn(f_ctau, Components("g1"), Name("g1"),LineStyle(kDotted), LineColor(kSpring));
  model.plotOn(f_ctau, Components("g2"), Name("g2"),LineStyle(kDotted), LineColor(kAzure));
  model.plotOn(f_ctau, Components("g3"), Name("g3"),LineStyle(kDotted), LineColor(kBlack));
  model.plotOn(f_ctau, Components("extNP"), Name("extNP"), LineStyle(kDashed), LineColor(kMagenta));
  model.plotOn(f_ctau, Components("extBkg"), Name("extBkg"), LineStyle(kDashed), LineColor(kOrange));
  model.plotOn(f_ctau, Components("extBkg,extNP"), Name("NPContribution"), LineStyle(kDashed), LineColor(kViolet));

  model.plotOn(f_ctau, Range("SR"), NormRange("SR"), VLines(), FillColor(kBlue - 3), FillStyle(3004), DrawOption("F"), Name("model_sub")); // Components("sig"),
  model.plotOn(f_ctau, Range("SR"), NormRange("SR"), VLines(), Components("extBkg,extNP"), FillColor(kViolet), FillStyle(3004), DrawOption("F"), Name("NP_sub"));
  model.plotOn(f_ctau, Range("SR"), NormRange("SR"), VLines(), Components("extBkg"), FillColor(kOrange), FillStyle(3004), DrawOption("F"), Name("bkg_sub"));

  // y axis: logY style
  double ymin = 1e300, ymax = -1e300;
  RooHist *hdata = (RooHist *)f_ctau->getHist("data"); // use first dataset on f_ctau
  if (hdata)
  {
    for (int i = 0; i < hdata->GetN(); i++)
    {
      double x, y;
      hdata->GetPoint(i, x, y);
      if (y > 0 && y < ymin)
        ymin = y;
      if (y > ymax)
        ymax = y;
    }
  }

  double floor = 1e-3;
  if (ymin <= 0 || ymin == 1e300)
    ymin = floor;

  f_ctau->SetMinimum(ymin * 0.5);
  f_ctau->SetMaximum(ymax * 500.0);

  // title
  f_ctau->GetYaxis()->SetTitle("Events");
  f_ctau->GetXaxis()->SetTitle("");
  f_ctau->Draw("e");

  // --- object legend ---
  TLegend *leg1 = new TLegend(0.44, 0.61, 0.58, 0.92);
  leg1->SetBorderSize(0);
  leg1->SetFillStyle(0);
  leg1->SetTextSize(0.025);

  leg1->AddEntry(f_ctau->findObject("data"), "Data", "lep");
  leg1->AddEntry(f_ctau->findObject("model"), "model", "pe");
  leg1->AddEntry(f_ctau->findObject("extPR"), "Prompt", "pe");
  leg1->AddEntry(f_ctau->findObject("g1"), "Gauss 1", "pe");
  leg1->AddEntry(f_ctau->findObject("g2"), "Gauss 2", "pe");
  leg1->AddEntry(f_ctau->findObject("g3"), "Gauss 3", "pe");
  leg1->Draw("same");

  TLegend *leg2 = new TLegend(0.56, 0.74, 0.65, 0.92);
  leg2->SetBorderSize(0);
  leg2->SetFillStyle(0);
  leg2->SetTextSize(0.025);
  leg2->AddEntry(f_ctau->findObject("NPContribution"), "NP Total", "pe");
  leg2->AddEntry(f_ctau->findObject("extNP"), "NP B #rightarrow J/#psi", "pe");
  leg2->AddEntry(f_ctau->findObject("extBkg"), "NP Continuum Bkg", "pe");
  leg2->Draw("same");

  // --- info latex ---
  TLatex latexInfo;
  latexInfo.SetNDC();
  latexInfo.SetTextSize(0.026);
  latexInfo.SetTextFont(42);

  double x_start = 0.19;
  double y_start = 0.95;
  double y_step = -0.06, y_stepCount = 1;
  latexInfo.DrawLatex(x_start, y_start + y_step * y_stepCount++, "CMS PbPb. #sqrt{s_{NN}} = 5.02 TeV");
  latexInfo.DrawLatex(x_start, y_start + y_step * y_stepCount++, Form("Data, %s %s", comp.c_str(), region.c_str()));
  if (yLow == 0)
    latexInfo.DrawLatex(x_start, y_start + y_step * y_stepCount++, Form("%.1f < p_{T} < %.1f, %.1f < y < %.1f", ptLow, ptHigh, yLow, yHigh));
  else
    latexInfo.DrawLatex(x_start, y_start + y_step * y_stepCount++, Form("%.1f < p_{T} < %.1f, |y| < %.1f", ptLow, ptHigh, yHigh));

  // --- latex parameters ---
  TLatex latexParams;
  latexParams.SetNDC();
  latexParams.SetTextSize(0.025);
  latexParams.SetTextFont(42);

  x_start = 0.76;
  y_step = -0.045;
  y_stepCount = 1;
  latexParams.DrawLatex(x_start, y_start + y_step * y_stepCount++, Form("N_{PR} = %.0f #pm %.0f", Npr.getVal(), Npr.getError()));
  latexParams.DrawLatex(x_start, y_start + y_step * y_stepCount++, Form("N_{NP} = %.0f #pm %.0f", Nnp.getVal(), Nnp.getError()));
  latexParams.DrawLatex(x_start, y_start + y_step * y_stepCount++, Form("N_{Bkg}/N_{NP} = %.2f (fixed)", NbkgRatio.getVal()));
  latexParams.DrawLatex(x_start, y_start + y_step * y_stepCount++, Form("#sigma1 = %.3f #pm %.3f", sigma1.getVal(), sigma1.getError()));
  latexParams.DrawLatex(x_start, y_start + y_step * y_stepCount++, Form("f_{G1+G2} = %.3f #pm %.3f", fsumVar->getVal(), fsumVar->getError()));
  // latexParams.DrawLatex(x_start, y_start + y_step * y_stepCount++, Form("f_{G1}/f_{G1+G2} = %.3f #pm %.3f", fg12Var->getVal(), fg12Var->getError()));
  latexParams.DrawLatex(x_start, y_start + y_step * y_stepCount++, Form("#tau^{NP}_{#psi} = %.3f #pm %.3f", tau.getVal(), tau.getError()));
  // latexParams.DrawLatex(x_start, y_start + y_step * y_stepCount++, Form("#tau^{NP}_{#psi} = %.3f #pm %.3f", tau.getVal(), tau.getError()));
  
  // subrange values
  latexParams.DrawLatex(x_start, y_start + y_step * y_stepCount++, Form("N_{PR, sub} = %.0f", Npr_in));
  latexParams.DrawLatex(x_start, y_start + y_step * y_stepCount++, Form("N_{NP, sub} = %.0f", Nnp_in));
  latexParams.DrawLatex(x_start, y_start + y_step * y_stepCount++, Form("N_{Bkg, sub} = %.0f", Nbkg_in));
  latexParams.DrawLatex(x_start, y_start + y_step * y_stepCount++, Form("f^{PR}_{#psi_{B}} = %.3f", b_fraction));


  // latexParams.DrawLatex(x_start, y_start + y_step * y_stepCount++, Form("f_{CB,R} = %.3f #pm %.3f", f_cbR.getVal(), f_cbR.getError()));
  // latexParams.DrawLatex(x_start, y_start + y_step * y_stepCount++, Form("f_{Gauss} = %.3f #pm %.3f", f_gaus.getVal(), f_gaus.getError()));
  // latexParams.DrawLatex(x_start, y_start + y_step * y_stepCount++, Form("c1 = %.3f #pm %.3f", c1.getVal(), c1.getError()));
  // latexParams.DrawLatex(x_start, y_start + y_step * y_stepCount++, Form("c1 = %.3f #pm %.3f", c2.getVal(), c2.getError()));

  // === pull pad ===
  c_ctau.cd();
  TPad *pad2 = new TPad("pad2", "pad2", 0.0, 0.0, 1.0, 0.25);
  pad2->SetTopMargin(0.00001);
  pad2->SetBottomMargin(0.4);
  pad2->Draw();
  pad2->cd();

  RooHist *hpull = f_ctau->pullHist("data", "model");
  RooPlot *f_pull = ctau->frame(Range(ctMin, ctMax), Title(""));
  f_pull->addPlotable(hpull, "P"); // P: points only

  f_pull->GetYaxis()->SetTitle("Pull");
  f_pull->GetXaxis()->SetTitle("c#tau_{3D} [mm]");
  f_pull->GetXaxis()->CenterTitle();
  f_pull->SetMinimum(-8);
  f_pull->SetMaximum(8);
  f_pull->GetYaxis()->SetNdivisions(505);
  f_pull->GetYaxis()->SetTitleSize(0.12);
  f_pull->GetYaxis()->SetLabelSize(0.10);
  f_pull->GetXaxis()->SetTitleSize(0.15);
  f_pull->GetXaxis()->SetLabelSize(0.10);
  f_pull->Draw();

  // --- draw pull = 0 line ---
  double xmin = ctMin;
  double xmax = ctMax;
  TLine *line = new TLine(xmin, 0.0, xmax, 0.0);
  // line->SetLineColor();
  line->SetLineStyle(2);
  line->Draw("same");

  // --- compute and draw chi square ---
  int nFitParam = fitResult->floatParsFinal().getSize();
  double chi2ndf = f_ctau->chiSquare("model", "data", nFitParam);

  TLatex latex;
  latex.SetNDC(); // use pad coordinates (0~1)
  latex.SetTextSize(0.1);
  latex.DrawLatex(0.82, 0.86, Form("#chi^{2}/ndf = %.2f", chi2ndf));

  c_ctau.SaveAs(Form("figs/ctau_%s_%s_pT%.1f_%.1f_y%.1f_%.1f.png", comp.c_str(), region.c_str(), ptLow, ptHigh, yLow, yHigh));
  c_ctau.SaveAs(Form("figs/ctau_%s_%s_pT%.1f_%.1f_y%.1f_%.1f.pdf", comp.c_str(), region.c_str(), ptLow, ptHigh, yLow, yHigh));

  fitResult->Print("V");

  cout << "Chi2/ndf: " << chi2ndf << "\n";


  // --- print # of events in subrange ---
  std::cout
      << "[-0.05,0.05] yields:\n"
      << "  Prompt     = " << Npr_in << "\n"
      << "  Nonprompt  = " << Nnp_in << "\n"
      << "  Background = " << Nbkg_in << "\n"
      << "  Total      = " << (Npr_in + Nnp_in + Nbkg_in) << "\n";

  // === save results ===
  TFile fout(Form("roots/ctau_%s_%s_pT%.1f_%.1f_y%.1f_%.1f.root", comp.c_str(), region.c_str(), ptLow, ptHigh, yLow, yHigh), "RECREATE");

  TParameter<double> b_frac("b_frac", b_fraction);
  b_frac.Write();
  fitResult->Write("fitResult");
  model.Write();
  fout.Close();

  cout << "=== Done ===\n";
  t.Stop();
  printf("RealTime=%f seconds, CpuTime=%f seconds\n", t.RealTime(), t.CpuTime());
}
