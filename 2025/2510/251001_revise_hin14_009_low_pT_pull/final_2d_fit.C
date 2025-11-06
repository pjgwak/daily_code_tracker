#include <TStopwatch.h>
#include <RooArgList.h>

using namespace RooFit;

// === kinematics === -> Read from config!
float ptLow = 13, ptHigh = 15;
float yLow = 0, yHigh = 2.4;
float massLow = 2.6, massHigh = 3.5;
int cLow = 0, cHigh = 180;
double ctMin = -2, ctMax = 6; // lmin, lmax: 1 for lowpT, 2 for higpT
float errmin = 0.008, errmax = 0.3;
const bool isWeight = true;

void defineMassBkg(RooWorkspace *ws);
void defineMassSig(RooWorkspace *ws);
void defineCtBkg(RooWorkspace *ws);
void defineCtPRRes(RooWorkspace *ws);
void readErrPdf(RooHistPdf *&sigPdf, RooHistPdf *&bkgPdf, TFile *&fInData);
void readInputs(RooDataSet *&data, const bool isWeight);
void setSilent();

void readFitResult(RooFitResult *&fitRes, TFile *&fin,
                   const std::string &fname,
                   const std::string &objname = "fitResult")
{
  string fileNameData = Form("roots/%s_pT%.1f_%.1f_y%.1f_%.1f_isWeight%d.root", fname.c_str(), ptLow, ptHigh, yLow, yHigh, isWeight);
	fin = TFile::Open(fileNameData.c_str(), "READ");
	if (!fin || fin->IsZombie()) {
			std::cerr << "[Error] Cannot open file: " << fileNameData << std::endl;
	}

	fitRes = dynamic_cast<RooFitResult *>(fin->Get(objname.c_str()));
	if (!fitRes) {
			std::cerr << "[Error] RooFitResult object '" << objname
								<< "' not found in file: " << fname << std::endl;
			fin->Close();
			fin = nullptr;
	}
}

void defineCtNP(RooWorkspace *ws)
{
  // --- simple ctau NP pdf ---
  // ws->factory("GaussModel::BResTrue(ctau3Dtrue,mean[0.0],sigmaNPTrue[0.00002,0.0000001,0.1])"); // resolution from B meson (before CtPRRes)
  // ws->factory("TruthModel::BResTrue(ctau3Dtrue)");
  // ws->factory("Decay::CtNPTrue(ctau3Dtrue,coefExpNPTrue[0.1,0.00001,10],BResTrue,RooDecay::SingleSided)");

  // --- two decay model ---
  // ws->factory("TruthModel::BResTrue(ctau3Dtrue)");
  // // ws->factory("GaussModel::BResTrue(ctau3Dtrue,mean[0.0],sigmaNPTrue[0.00002,0.0000001,0.1])"); // resolution from B meson (before CtPRRes)

  ws->factory("coefExpNPTrue1[0.10, 0.001, 1.0]");
  ws->factory("rCoefExpNPTrue[2.0, 1.01, 50.0]");
  ws->factory("expr::coefExpNPTrue2('coefExpNPTrue1/rCoefExpNPTrue',{coefExpNPTrue1,rCoefExpNPTrue})");

  // ws->factory("Decay::CtNPTrue1(ctau3Dtrue, coefExpNPTrue1, BResTrue, RooDecay::SingleSided)");
  // ws->factory("Decay::CtNPTrue2(ctau3Dtrue, coefExpNPTrue2, BResTrue, RooDecay::SingleSided)");

  // ws->factory("SUM::CtNPTrue(fFrac[0.5,0.0,1.0]*CtNPTrue1, CtNPTrue2)");

  // // --- three decay model ---
  // ws->factory("TruthModel::BResTrue(ctau3Dtrue)");
  // ws->factory("coefExpNPTrue1[0.1,0.001,1]");

  // ws->factory("r21[3.0,1,10.0]");   // lambda2 / lambda1
  // ws->factory("r31[6.0,1,20.0]");   // lambda3 / lambda1

  // ws->factory("expr::coefExpNPTrue2('coefExpNPTrue1*r21', coefExpNPTrue1,r21)");
  // ws->factory("expr::coefExpNPTrue3('coefExpNPTrue2*r31', coefExpNPTrue2,r31)");

  // // Decay PDFs
  // ws->factory("Decay::CtNPTrue1(ctau3Dtrue, coefExpNPTrue1, BResTrue, RooDecay::SingleSided)");
  // ws->factory("Decay::CtNPTrue2(ctau3Dtrue, coefExpNPTrue2, BResTrue, RooDecay::SingleSided)");
  // ws->factory("Decay::CtNPTrue3(ctau3Dtrue, coefExpNPTrue3, BResTrue, RooDecay::SingleSided)");

  // ws->factory("SUM::CtNPTrue(fFrac1[0.3,0.0,1.0]*CtNPTrue1, fFrac2[0.3,0.0,1.0]*CtNPTrue2, CtNPTrue3)");

  // --- build NP Pdf convoluted PRRes ---
  // RooFormulaVar sigmaNPResW("sigmaNPResW", "sqrt((@0*@1)**2+(@2)**2)", RooArgList(*(ws->var("sigmaPRResW")), *(ws->var("ctau3DErr")), *(ws->var("sigmaNPTrue"))));
  // ws->import(sigmaNPResW);
  // RooFormulaVar sigmaNPResN("sigmaNPResN", "sqrt((@0*@1)**2+(@2)**2)", RooArgList(*(ws->var("sigmaPRResN")), *(ws->var("ctau3DErr")), *(ws->var("sigmaNPTrue"))));
  // ws->import(sigmaNPResN);
  // ws->factory("GaussModel::GW_NPRes(ctau3D,meanPRResW,sigmaNPResW)");
  // ws->factory("GaussModel::GN_NPRes(ctau3D,meanPRResN,sigmaNPResN)");

	ws->factory("GaussModel::GW_NPRes(ctau3D,meanPRResW,sigmaPRResW,ctau3DErr)");
  ws->factory("GaussModel::GN_NPRes(ctau3D,meanPRResN,sigmaPRResN, ctau3DErr)");
  ws->factory("AddModel::CtNPRes({GW_NPRes,GN_NPRes},{fracRes})");

  // final model
  // float coefExpNPTrueVal = ws->var("coefExpNPTrue")->getVal();
  ws->factory("Decay::CtNPTot1(ctau3D,coefExpNPTrue1,CtNPRes,RooDecay::SingleSided)");
  ws->factory("Decay::CtNPTot2(ctau3D,coefExpNPTrue2,CtNPRes,RooDecay::SingleSided)");
  ws->factory("SUM::CtNPTot(fFracCtNPTot[0.5,0.0,1.0]*CtNPTot1, CtNPTot2)");

  // **** check NPMC Reco
  // 이거 왜 있지????
  // drawCtNPReco(ws, redNPMC, titlestr, opt);
  // return;
}


void fixParams(RooWorkspace *ws, RooFitResult *fitResult, const vector<string> &fixList)
{
  if (!ws || !fitResult)
  {
    cerr << "[Error] fixParams: null input" << "\n";
    return;
  }

	const RooArgList &floatPars = fitResult->floatParsFinal();
	const RooArgList &constPars = fitResult->constPars();

	RooArgList pars;
	pars.add(floatPars);
	pars.add(constPars);

  for (const auto &name : fixList)
  {
    RooRealVar *varWS = ws->var(name.c_str());
    if (!varWS)
    {
      cerr << "[Warning] variable '" << name << "' not found in workspace" << "\n";
      continue;
    }

    RooRealVar *varFit = dynamic_cast<RooRealVar *>(pars.find(name.c_str()));
    if (!varFit)
    {
      cerr << "[Warning] variable '" << name << "' not found in fitResult" << "\n";
      continue;
    }

    // set value and fix
    varWS->setVal(varFit->getVal());
    varWS->setConstant(true);

    std::cout << "[Info] variable '" << name
              << "' set to " << varFit->getVal()
              << " and fixed." << "\n";
  }
}

void drawFinalCtau(RooWorkspace *ws, RooDataSet *redData, RooFitResult *fitRes, std::string outname="finalCtau")
{
  //----------------------------------------------
  RooRealVar *mass = ws->var("mass");
  mass->setRange("massRange", 2.6, 3.5);
  RooRealVar *ctau3DErr = ws->var("ctau3DErr");
  ctau3DErr->setRange("massRange", 0.008, 0.152);
  RooRealVar *ctau3D = ws->var("ctau3D");
  if (!ctau3D) { std::cout << "[Error] no variable ctau3D" << std::endl; return; }


  RooPlot *frame = ctau3D->frame(Title("Background fit"));
  //  frame->updateNormVars(RooArgSet(*(ws->var("ctau3D")), *(ws->var("mass"))));
  // frame->updateNormVars(RooArgSet(*mass, *ctau3D, *ctau3DErr));

  // ws->pdf("CtBkgTot_PEE")->plotOn(frame, Name("fit"), LineColor(kBlue), Normalization(redData->sumEntries(), RooAbsReal::NumEvent), Name("fit"));

  

  // ws->pdf("totPDF_PEE")->plotOn(frame,LineStyle(kDashed), NumCPU(32), LineColor(kBlue), ProjWData(RooArgSet(*(ws->var("ctau3DErr"))), *redData), Name("model"), Normalization(1, RooAbsReal::NumEvent));

  // ws->pdf("totPDF_PEE")->plotOn(frame, Name("model"), NumCPU(32), ProjWData(RooArgList(*ctau3DErr), *redData), LineColor(kBlue), Normalization(1, RooAbsReal::NumEvent), Precision(1e-4)); // redData->sumEntries(),

  // ws->pdf("totPDF_PEE")->plotOn(frame, Name("model"), NumCPU(32), ProjWData(RooArgList(*ctau3DErr), *redData, kTRUE), LineColor(kBlue), Normalization(redData->sumEntries(), RooAbsReal::NumEvent));

  redData->plotOn(frame, Name("data")); //, DataError(RooAbsData::SumW2)

  // redData->sumEntries(), ConditionalObservables(RooArgSet(*ctau3DErr), RooArgList(*)
  // === pull hist ===
  // RooHist *hpull = frame->pullHist("data","model");
  // RooPlot *framePull = ctau3D->frame(Title("Pull distribution"));
  // framePull->addPlotable(hpull,"P");

  TCanvas *c = new TCanvas("c","c",800,800);
  c->Divide(1,2);

  c->cd(1);
  gPad->SetPad(0.0,0.3,1.0,1.0);
  gPad->SetBottomMargin(0.02);
  gPad->SetLogy();
  frame->SetMinimum(0.1);
  frame->Draw();

  // double chi2 = frame->chiSquare("model","data");
  // TLatex lat;
  // lat.SetNDC();
  // lat.SetTextSize(0.05);
  // lat.DrawLatex(0.65,0.85,Form("#chi^{2}/ndf = %.2f",chi2));

  c->cd(2);
  gPad->SetPad(0.0,0.0,1.0,0.3);
  gPad->SetTopMargin(0.05);
  gPad->SetBottomMargin(0.25);
  // framePull->GetYaxis()->SetTitle("Pull");
  // framePull->GetYaxis()->SetTitleSize(0.12);
  // framePull->GetYaxis()->SetLabelSize(0.10);
  // framePull->GetXaxis()->SetTitleSize(0.12);
  // framePull->GetXaxis()->SetLabelSize(0.10);
  // framePull->Draw();

  c->SaveAs(Form("%s.png", outname.c_str()));
  c->SaveAs(Form("%s.pdf", outname.c_str()));
}

void drawFinalMass(RooWorkspace *ws, RooDataSet *redData, RooFitResult *fitRes, std::string outname = "finalCtau")
{
  RooRealVar *mass = ws->var("mass");
  if (!mass)
  {
    std::cout << "[Error] no variable mass" << std::endl;
    return;
  }

  RooPlot *frame = mass->frame(Title("Background fit"));
  //  frame->updateNormVars(RooArgSet(*(ws->var("mass")), *(ws->var("mass"))));

  // ws->pdf("CtBkgTot_PEE")->plotOn(frame, Name("fit"), LineColor(kBlue), Normalization(redData->sumEntries(), RooAbsReal::NumEvent), Name("fit"));
  redData->plotOn(frame, Name("data"));

  // ws->pdf("totPDF_PEE")->plotOn(frame,LineStyle(kDashed), NumCPU(32), LineColor(kBlue), ProjWData(RooArgSet(*(ws->var("massErr"))), *redData), Name("model"), Normalization(1, RooAbsReal::NumEvent));

  //   ws->var("mass")->setBins(150);  // bin 개수 설정
  // auto h_ctau = new RooDataHist("h_ctau", "binned mass data",
  //                    RooArgSet(*(ws->var("mass"))), *redData);

  // 히스토그램을 프레임에 그림
  // h_ctau->plotOn(frame, DataError(RooAbsData::SumW2), Name("data"));
  ws->pdf("totPDF_PEE")->plotOn(frame, Name("model"), NumCPU(32), ProjWData(*redData), LineColor(kBlue), Normalization(1, RooAbsReal::NumEvent));

  // === pull hist ===
  RooHist *hpull = frame->pullHist("data", "model");
  RooPlot *framePull = mass->frame(Title("Pull distribution"));
  framePull->addPlotable(hpull, "P");

  TCanvas *c = new TCanvas("c", "c", 800, 800);
  c->Divide(1, 2);

  c->cd(1);
  gPad->SetPad(0.0, 0.3, 1.0, 1.0);
  gPad->SetBottomMargin(0.02);
  gPad->SetLogy();
  frame->SetMinimum(0.1);
  frame->Draw();

  double chi2 = frame->chiSquare("model", "data");
  TLatex lat;
  lat.SetNDC();
  lat.SetTextSize(0.05);
  lat.DrawLatex(0.65, 0.85, Form("#chi^{2}/ndf = %.2f", chi2));

  c->cd(2);
  gPad->SetPad(0.0, 0.0, 1.0, 0.3);
  gPad->SetTopMargin(0.05);
  gPad->SetBottomMargin(0.25);
  framePull->GetYaxis()->SetTitle("Pull");
  framePull->GetYaxis()->SetTitleSize(0.12);
  framePull->GetYaxis()->SetLabelSize(0.10);
  framePull->GetXaxis()->SetTitleSize(0.12);
  framePull->GetXaxis()->SetLabelSize(0.10);
  framePull->Draw();

  c->SaveAs(Form("%s.png", outname.c_str()));
  c->SaveAs(Form("%s.pdf", outname.c_str()));
}

void final_2d_fit()
{
	cout << "\n=== Start ctau_bkg_fit() ===\n";
  TStopwatch t;
  t.Start();

  // silent
  setSilent();

  // rootlogon -> To global config
  gROOT->Macro("/data/users/pjgwak/input_files/rootlogon.C");

  // === import inputs === 
  // redData, errSig, errBkg, mass
	// NP model, ctauBkg - 이 안에 Res 포함
  cout << "=== Import inputs ===\n";
  RooDataSet *redData = nullptr;
  readInputs(redData, isWeight);

  // get errPdfSig
  TFile *fInErr = nullptr;
  RooHistPdf *sigPdf = nullptr, *bkgPdf = nullptr;
  readErrPdf(sigPdf, bkgPdf, fInErr);

	// ctau true
  TFile *fInCTrue = nullptr;
  RooFitResult *fitCTrue = nullptr;
  readFitResult(fitCTrue, fInCTrue, "ctau_true", "fitCt_NPMC");

	// ctauBkg
	TFile *fInCBkg = nullptr;
  RooFitResult *fitCBkg = nullptr;
  readFitResult(fitCBkg, fInCBkg, "ctau_bkg", "fitCt_Bkg");

	// mass
	TFile *fInMass = nullptr;
  RooFitResult *fitMass = nullptr;
  readFitResult(fitMass, fInMass, "mass", "fitMass");
	

  // === make workspace ===
  auto ws = new RooWorkspace("ws");
  ws->import(*redData);
  ws->import(*bkgPdf);
	ws->import(*sigPdf);

  // === define model ===
	cout << "\n=== build model ===\n";
	// mass
	defineMassSig(ws);
  defineMassBkg(ws);

  char funct[100];
  sprintf(funct, "SUM::MassPDF(NSig[%f,1.0,5000000.0]*%s, NBkg[%f,1.0,50000000.0]*%s)", 1000000., "G1CB1Sig", 1000000., "expBkg"); // factory 문법으로 바꾸기
  ws->factory(funct);

	// ctauBkg
	defineCtPRRes(ws);
  defineCtBkg(ws);

	// ctauTrue
	defineCtNP(ws);

  // // Build 2D PDF (mass x ctau)
	RooProdPdf CtPR_PEE("CtPR_PEE", "CtPDF with PEE", *(ws->pdf("errPdfSig")), Conditional(*(ws->pdf("CtPRRes")), RooArgList(*(ws->var("ctau3D")))));
	RooProdPdf CtBkgTot_PEE("CtBkgTot_PEE", "PDF with PEE", *(ws->pdf("errPdfBkg")), Conditional(*(ws->pdf("CtBkgTot")), RooArgList(*(ws->var("ctau3D")))));
  
	sprintf(funct, "PROD::MassCtPR(%s,CtPRRes)", "G1CB1Sig"); ws->factory(funct);
  sprintf(funct, "PROD::MassCtNP(%s,CtNPTot)", "G1CB1Sig"); ws->factory(funct);
  sprintf(funct, "PROD::MassCtBkg(%s,CtBkgTot)", "expBkg"); ws->factory(funct);

  RooProdPdf MassCtPR_PEE("MassCtPR_PEE", "PDF with PEE", *(ws->pdf("errPdfSig")),
                          Conditional(*(ws->pdf("MassCtPR")), RooArgList(*(ws->var("ctau3D")), *(ws->var("mass")))));
  ws->import(MassCtPR_PEE);
  RooProdPdf MassCtNP_PEE("MassCtNP_PEE", "PDF with PEE", *(ws->pdf("errPdfSig")),
                          Conditional(*(ws->pdf("MassCtNP")), RooArgList(*(ws->var("ctau3D")), *(ws->var("mass")))));
  ws->import(MassCtNP_PEE);
  RooProdPdf MassCtBkg_PEE("MassCtBkg_PEE", "PDF with PEE", *(ws->pdf("errPdfBkg")),
                           Conditional(*(ws->pdf("MassCtBkg")), RooArgList(*(ws->var("ctau3D")), *(ws->var("mass")))));
  ws->import(MassCtBkg_PEE);

  // Total fit PDF = promptPDF + nonpromptPDF + bkgPDF
  // 앞의 둘은 mass 시그널에서 개수를 정하고 마지막은 mass bkg에서 개수를 정한다.
  RooFormulaVar fracBkg("fracBkg", "@0/(@0+@1)", RooArgList(*(ws->var("NBkg")), *(ws->var("NSig"))));
  ws->import(fracBkg);
  ws->factory("RSUM::totPDF_PEE(fracBkg*MassCtBkg_PEE,Bfrac[0.25,0.01, 0.99]*MassCtNP_PEE,MassCtPR_PEE)");
  // ws->factory("SUM::totPDF_PEE(fracBkg*MassCtBkg_PEE, Bfrac[0.25,0.0,1.]*MassCtNP_PEE, MassCtPR_PEE)");


	// === fix parameters ===
	cout << "\n=== fix parameters ===\n";
	// mass
  vector<string> fixList = {"NBkg", "NSig", "coefExp", "fracG1", "meanSig", "sigmaSig1", "alpha", "enne", "sigmaSig2", "alpha1", "enne1", "alpha2", "enne2"};
  fixParams(ws, fitMass, fixList);
	// ctauBkg, Res
  fixList = {"fracCtBkg1", "fracCtBkg2", "fracRes", "lambdap", "lambdasym", "sigmaPRResW", "lambdam"}; // fracCtBkg3, , sigmaPRResN
  fixParams(ws, fitCBkg, fixList);

	// ctauTrue
	fixList = {"rCoefExpNPTrue", "fFracCtNPTot"}; //coefExpNPTrue1
  fixParams(ws, fitCTrue, fixList);

	// === perfrom fit ===
	cout << "\n=== perfrom fit ===\n";
  RooFitResult *fit2D = ws->pdf("totPDF_PEE")->fitTo(*redData, Minos(0), Save(1), SumW2Error(kTRUE), PrintEvalErrors(-1), ConditionalObservables(RooArgSet(*(ws->var("ctau3DErr")))), NumCPU(32));
  // RecoverFromUndefinedRegions(1), 
  // // NumCPU(8), EvalBackend("legacy"),
  // // 8 cpu, No EvalBAckend - 389 s
  // // 16 cpu, No EvalBAckend - 374s
  // // 8 cpu, EvalBAckend(legacy) - 403s
  // // 8 cpu, EvalBAckend(cpu) 378 ss{

	// === draw plot ===
  cout << "\n=== draw plot ===\n";
  // drawFinalMass(ws, redData, fit2D, "fitFinalMass");
  // drawFinalMass(ws, redData, NSigNP_fin, NBkg_fin, fitMass, &UnNormChi2_mass, &nFitParam_mass, &nFullBinsPull_mass, &Dof_mass, &Chi2_mass);
  // drawFinalCtau(ws, redData, binDataCtErr, NSigNP_fin, NBkg_fin, Bfrac_fin, ErrBfrac_fin, fit2D, -ctMin, ctMax, &UnNormChi2_time, &nFitParam_time, &nFullBinsPull_time, &Dof_time, &Chi2_time);
  drawFinalCtau(ws, redData, fit2D, "fitFinalCtau");


  // === save result ===
  // cout << "\n=== save resuts ===\n";
  // gSystem->mkdir("roots", true);
  // TFile fileOut(Form("roots/final_2d_pT%.1f_%.1f_y%.1f_%.1f_isWeight%d.root", ptLow, ptHigh, yLow, yHigh, isWeight), "recreate");
  // fit2D->Write("fit2D");
  // fileOut.Close();

  fit2D->Print("v");

  cout << "\n=== Finish ctau_bkg_fit() ===\n";
  t.Stop();
  printf("RealTime=%f seconds, CpuTime=%f seconds\n", t.RealTime(), t.CpuTime());
}


void setSilent()
{
  RooMsgService::instance().getStream(0).removeTopic(RooFit::Plotting);
  RooMsgService::instance().getStream(1).removeTopic(RooFit::Plotting);
  // RooMsgService::instance().getStream(0).removeTopic(InputArguments);
  // RooMsgService::instance().getStream(1).removeTopic(InputArguments);
  // RooMsgService::instance().getStream(1).removeTopic(Plotting);
  RooMsgService::instance().getStream(1).removeTopic(NumIntegration);
  RooMsgService::instance().getStream(0).removeTopic(RooFit::Integration);
  RooMsgService::instance().getStream(1).removeTopic(RooFit::Integration);
  // RooMsgService::instance().getStream(1).removeTopic(Minimization);
  RooMsgService::instance().getStream(1).removeTopic(Caching);
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING); // only print from WARNING to FATAL
}

void readInputs(RooDataSet *&data, const bool isWeight)
{
  // --- data ---
  string fileNameData = Form("reduced_data/data_pT%.1f_%.1f_y%.1f_%.1f_isWeight%d.root", ptLow, ptHigh, yLow, yHigh, isWeight);
  TFile fInData(fileNameData.c_str());
  cout << "Load: " << fileNameData.c_str() << endl;
  data = (RooDataSet *)fInData.Get("redData");
  if (isWeight)
  {
    cout << "[Info] data is WEIGHTED\n";
  }
  else
    cout << "[Info] data is UN-WEIGHTED" << endl;
  data->SetName("redData");
}

void readErrPdf(RooHistPdf *&sigPdf, RooHistPdf *&bkgPdf, TFile *&fInData)
{
  string fileNameData = Form("roots/ctau_err_pT%.1f_%.1f_y%.1f_%.1f_isWeight%d.root", ptLow, ptHigh, yLow, yHigh, isWeight);
  fInData = TFile::Open(fileNameData.c_str());
  cout << "Load: " << fileNameData.c_str() << endl;
	sigPdf = dynamic_cast<RooHistPdf *>(fInData->Get("errPdfSig"));
  bkgPdf = dynamic_cast<RooHistPdf *>(fInData->Get("errPdfBkg"));
}

void defineCtPRRes(RooWorkspace *ws)
{
  ws->factory("GaussModel::GW_PRRes(ctau3D,meanPRResW[0],sigmaPRResW[2.3, 0.001, 10],one[1],ctau3DErr)");
  ws->factory("GaussModel::GN_PRRes(ctau3D,meanPRResN[0],sigmaPRResN[0.8,0.01, 3],one,ctau3DErr)");
  ws->factory("AddModel::CtPRRes({GW_PRRes,GN_PRRes},{fracRes[0.5,0.01,0.999]})");

  // ws->factory("GaussModel::GW_PRRes(ctau3D,meanPRResW[0],sigmaPRResW[2.3, 0.001, 10],one[1.0],ctau3DErr)");
  // ws->factory("GaussModel::GN_PRRes(ctau3D,meanPRResN[0],sigmaPRResN[0.8,0.001, 10],one,ctau3DErr)");
  // ws->factory("GaussModel::GVW_PRRes(ctau3D,meanPRResVW[0],sigmaPRResVW[0.8,0.001, 10],one,ctau3DErr)");
  // ws->factory("AddModel::CtPRRes({GN_PRRes, GW_PRRes, GVW_PRRes},{fracResN[0.01,0.001,0.999], fracResW[0.01,0.001,0.999]})");
  return;
}

void defineCtBkg(RooWorkspace *ws)
{
  // ws->factory("Decay::CtBkgPos(ctau3D,lambdap[0.42,0.05,1.5],CtPRRes,RooDecay::SingleSided)");
  // ws->factory("Decay::CtBkgNeg(ctau3D,lambdam[0.79,0.02,1.5],CtPRRes,RooDecay::Flipped)");
  // ws->factory("Decay::CtBkgDbl(ctau3D,lambdasym[0.69,0.02,5.0],CtPRRes,RooDecay::DoubleSided)");

  // nominal
  ws->factory("Decay::CtBkgPos(ctau3D,lambdap[0.02,0.001, 2],CtPRRes,RooDecay::SingleSided)");
  ws->factory("Decay::CtBkgNeg(ctau3D,lambdam[0.02,0.001, 2],CtPRRes,RooDecay::Flipped)");
  ws->factory("Decay::CtBkgDbl(ctau3D,lambdasym[0.2,0.001, 2],CtPRRes,RooDecay::DoubleSided)");
  ws->factory("SUM::CtBkgSum1(fracCtBkg1[0.3,0.01,1]*CtBkgPos,CtBkgNeg)");
  ws->factory("SUM::CtBkgSum2(fracCtBkg2[0.1,0.01,1]*CtBkgSum1,CtBkgDbl)");
  ws->factory("SUM::CtBkgTot(fracCtBkg3[0.1,0.01,1]*CtPRRes,CtBkgSum2)");

  // // // --- 2 left decay components ---
  // // 양의 방향 (SingleSided)
  // ws->factory("Decay::CtBkgPos(ctau3D, lambdap[0.02,0.001, 2], CtPRRes, RooDecay::SingleSided)");

  // // 음의 방향 (Flipped)
  // ws->factory("Decay::CtBkgNeg(ctau3D, lambdam[0.02,0.001, 2], CtPRRes, RooDecay::Flipped)");

  // // 새로운 Left 성분 (SingleSided, 독립적인 lambda 사용)
  // ws->factory("Decay::CtBkgLeft(ctau3D, lambdaleft[0.05,0.001, 2], CtPRRes, RooDecay::DoubleSided)");

  // // 대칭(DoubleSided)
  // ws->factory("Decay::CtBkgDbl(ctau3D, lambdasym[0.2,0.001, 2], CtPRRes, RooDecay::DoubleSided)");
  // ws->factory("SUM::CtBkgSum1(fracCtBkg1[0.3,0.01,1]*CtBkgPos, CtBkgNeg)");
  // ws->factory("SUM::CtBkgSum2(fracCtBkgL[0.2,0.01,1]*CtBkgLeft, CtBkgSum1)");

  // ws->factory("SUM::CtBkgSum3(fracCtBkg2[0.1,0.01,1]*CtBkgSum2, CtBkgDbl)");

  // ws->factory("SUM::CtBkgTot(fracCtBkg3[0.1,0.01,1]*CtPRRes, CtBkgSum3)");

  return;
}

void defineMassSig(RooWorkspace *ws)
{
  // --- Gauss + CB ---
  // ws->factory("Gaussian::G1Sig(mass,meanSig[3.0975,3.05,3.15],sigmaSig1[0.03,0.008,0.075])");
  // ws->factory("CBShape::CB1Sig(mass,meanSig,sigmaSig2[0.03,0.0008,0.075],alpha[1.9,1.2,2.8],enne[2.5,1.0,4.0])");
  // ws->factory("SUM::G1CB1Sig(fracG1[0.5,0.01,0.99]*G1Sig,CB1Sig)");

  // --- DCB ---
  ws->factory("CBShape::CB1Sig(mass,meanSig[3.0975,3.05,3.15],sigmaSig1[0.03,0.0008,0.075],alpha1[1.9,1.2,2.8],enne1[2.5,1.0,4.0])");
  ws->factory("CBShape::CB2Sig(mass,meanSig,sigmaSig2[0.03,0.0008,0.075],alpha2[1.9,1.2,2.8],enne2[2.5,1.0,4.0])");

  ws->factory("SUM::G1CB1Sig(fracG1[0.5,0.01,0.99]*CB1Sig,CB2Sig)");
  return;
}

void defineMassBkg(RooWorkspace *ws)
{
  // 1st order polynomial
  ws->factory("Polynomial::polBkg(mass,{coefPol[-0.05,-5.,5.]})");

  // expo
  ws->factory("Exponential::expBkg(mass,coefExp[-1.,-3.,2.])");
  return;
}