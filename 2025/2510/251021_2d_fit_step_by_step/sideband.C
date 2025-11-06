#include <iostream>
#include <TStopwatch.h>
#include <string.h>
#include <RooDataSet.h>
#include <TFile.h>
#include <RooRealVar.h>
#include <RooCBShape.h>
#include <RooAddPdf.h>
#include <RooGaussian.h>
#include <RooFitResult.h>
#include <RooPlot.h>
#include <RooHist.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TAxis.h>
#include <RooChebychev.h>
#include <TH1D.h>
#include <TMath.h>
#include <RooHistPdf.h>
#include <TROOT.h>   // gROOT
#include <TSystem.h> // gSystem
#include <TLine.h>
// #include <>
// #include <>
// #include <>
// #include <>
// #include <>
// #include <>
// #include <>
// #include <>
// #include <>
// #include <>
// #include <>
// #include <>
// #include <>
// #include <>
// #include <>

using namespace RooFit;

using std::cout;
using std::string;

// === define helper functions ===
// --- calculate mass bkg scale value ---
inline double computeScaleF_cheby(const RooArgList &allPars,
                                  int order, // 1~6
                                  double massMin, double massMax,
                                  double sbL_lo, double sbL_hi,
                                  double sig_lo, double sig_hi,
                                  double sbR_lo, double sbR_hi,
                                  const char *coeffPrefix = "sl") // name of chebychev parameters
{
  // --- helpers ---
  auto getVal = [&](const char *name) -> double
  {
    if (auto *obj = allPars.find(name))
    {
      if (auto *rr = dynamic_cast<RooAbsReal *>(obj))
        return rr->getVal();
    }
    return std::numeric_limits<double>::quiet_NaN();
  };

  auto toUnitX = [&](double m) -> double
  {
    return 2.0 * (m - massMin) / (massMax - massMin) - 1.0; // [massMin,massMax] → [-1,1]
  };

  auto chebT = [&](int n, double x) -> double
  {
    if (n == 0)
      return 1.0;
    if (n == 1)
      return x;
    double Tn_2 = 1.0, Tn_1 = x, Tn = 0.0;
    for (int k = 2; k <= n; ++k)
    {
      Tn = 2.0 * x * Tn_1 - Tn_2;
      Tn_2 = Tn_1;
      Tn_1 = Tn;
    }
    return Tn;
  };

  auto chebVal_roofit = [&](double m, const std::vector<double> &c) -> double
  {
    // RooChebychev definition: f = 1*T0 + c1*T1 + ... + cN*TN  (상수항 1 내장)
    const double x = toUnitX(m);
    double v = 1.0;
    for (int i = 1; i <= int(c.size()); ++i)
      v += c[i - 1] * chebT(i, x);
    return v;
  };

  auto integrateSimpson = [&](auto &&f, double a, double b, int nSeg = 2048) -> double
  {
    if (b <= a)
      return 0.0;
    if (nSeg % 2)
      ++nSeg;
    const double h = (b - a) / nSeg;
    double s = f(a) + f(b);
    for (int i = 1; i < nSeg; ++i)
    {
      const double x = a + i * h;
      s += f(x) * ((i % 2) ? 4.0 : 2.0);
    }
    return s * (h / 3.0);
  };

  if (order < 0 || order > 6)
  {
    ::Error("computeScaleF_cheby", "Unsupported Cheby order %d (allowed 0..6).", order);
    return std::numeric_limits<double>::quiet_NaN();
  }

  std::vector<double> c;
  c.reserve(order);
  for (int i = 1; i <= order; ++i)
  {
    char nm[32];
    std::snprintf(nm, sizeof(nm), "%s%d", coeffPrefix, i); // sl1, sl2 ...
    double v = getVal(nm);
    if (!std::isfinite(v))
    {
      ::Error("computeScaleF_cheby", "Coefficient '%s' not found or NaN.", nm);
      return std::numeric_limits<double>::quiet_NaN();
    }
    c.push_back(v);
  }

  // calculation
  auto f = [&](double m)
  { return chebVal_roofit(m, c); };
  const double I_sig = integrateSimpson(f, sig_lo, sig_hi);
  const double I_sb = integrateSimpson(f, sbL_lo, sbL_hi) + integrateSimpson(f, sbR_lo, sbR_hi);

  if (std::abs(I_sb) < 1e-300)
  {
    ::Error("computeScaleF_cheby", "Sideband integral is zero. Check ranges.");
    return std::numeric_limits<double>::quiet_NaN();
  }
  return I_sig / I_sb;
}

// --- subtract histograms ---
// - binScaleBKG : k*SB
// - binSubtrSIG : SIG - k*SB (가중치+오차)
// - dk          : scale factor's error - 0
// - floorW      : min value to avoid negative, 0 values
inline void subtractSidebands(RooDataHist *binSubtrSIG,   // out
                                                    RooDataHist *binScaleBKG,   // out
                                                    const RooDataHist *binSIG,  // in
                                                    const RooDataHist *binSB,   // in
                                                    RooRealVar *obsVar,         // e.g. ctau3DErr
                                                    double k,                   // scaleF
                                                    const std::string &varName, // "ctau3DErr"
                                                    double dk = 0.0,
                                                    double floorW = 1e-6)
{
  if (!binSubtrSIG || !binScaleBKG || !binSIG || !binSB || !obsVar)
  {
    ::Error("subtractSidebands", "Null input.");
    return;
  }
  const int nSIG = binSIG->numEntries();
  const int nSB = binSB->numEntries();
  if (nSIG != nSB)
  {
    ::Error("subtractSidebands", "Different binning: %d vs %d", nSIG, nSB);
    return;
  }

  std::unique_ptr<TH1> hSIG(binSIG->createHistogram("hSIG_tmp", *obsVar));
  std::unique_ptr<TH1> hSB(binSB->createHistogram("hSB_tmp", *obsVar));
  if (!hSIG || !hSB)
  {
    ::Error("subtractSidebands", "createHistogram failed.");
    return;
  }
  const int nb = hSIG->GetNbinsX();
  if (nb != hSB->GetNbinsX() || nb != nSIG)
  {
    ::Error("subtractSidebands", "Bin mismatch: TH1(%d) vs RooDataHist(%d).", nb, nSIG);
    return;
  }

  for (int i = 0; i < nSIG; ++i)
  {
    const RooArgSet *argSIG = binSIG->get(i);
    const RooArgSet *argSB = binSB->get(i);
    if (!argSIG || !argSB)
      continue;
    if (!argSIG->find(varName.c_str()) || !argSB->find(varName.c_str()))
    {
      ::Error("subtractSidebands", "Var '%s' not found in argset.", varName.c_str());
      return;
    }

    const int j = i + 1;                        
    const double wSIG = hSIG->GetBinContent(j); 
    const double eSIG = hSIG->GetBinError(j);

    const double wSB = hSB->GetBinContent(j);
    const double eSB = hSB->GetBinError(j);

    // k*SB
    const double wSBs = k * wSB;
    const double varSBs = (k * k) * (eSB * eSB) + (dk > 0.0 ? (wSB * wSB) * (dk * dk) : 0.0);
    const double eSBs = std::sqrt(std::max(0.0, varSBs));

    // SIG - k*SB
    const double wSUB = wSIG - wSBs;
    const double varSUB = (eSIG * eSIG) + varSBs;
    const double eSUB = std::sqrt(std::max(0.0, varSUB));

    const double wSBs_floor = (wSBs <= floorW ? floorW : wSBs);
    const double wSUB_floor = (wSUB <= floorW ? floorW : wSUB);

    binScaleBKG->set(*argSIG, wSBs_floor, eSBs);
    binSubtrSIG->set(*argSIG, wSUB_floor, eSUB);
  }
}

void sideband()
{
  cout << "=== start sideband() ===\n";
  TStopwatch time;
  time.Start();

  gROOT->Macro("/data/users/pjgwak/input_files/rootlogon.C");

  // set variables
  float ptLow = 6.5, ptHigh = 7.5;
  float yLow = 0, yHigh = 2.4;
  float massLow = 2.6, massHigh = 3.5;
  int cLow = 0, cHigh = 180;

  TString userLabel = "";
  TString figDir = Form("figs%s/pT%.1f_%.1f_y%.1f_%.1f", userLabel.Data(), ptLow, ptHigh, yLow, yHigh);
  gSystem->mkdir(figDir, kTRUE);
  TString rootDir = Form("roots%s/pT%.1f_%.1f_y%.1f_%.1f", userLabel.Data(), ptLow, ptHigh, yLow, yHigh);
  gSystem->mkdir(rootDir, kTRUE);

  // observable ranges
  double ctMin = -1, ctMax = 4; // lmin, lmax: 1 for lowpT, 2 for higpT
  float errmin = 0.0, errmax = 1.0;

  // sideband ranges
  double sbL_lo=2.6, sbL_hi=2.9, sig_lo=2.9, sig_hi=3.3, sbR_lo=3.3, sbR_hi=3.5;

  // read inputs
  string fileNameData = "/data/hwan/psi2S_RAA_PbPb2018/skimmedFiles/OniaRooDataSet_isMC0_JPsi_pp_y0.00_2.40_Effw0_Accw0_PtW1_TnP1_230221.root";
  TFile inputData(fileNameData.c_str());
  RooDataSet *dataPRMC = (RooDataSet *)inputData.Get("dataset");
  dataPRMC->SetName("dataPRMC");

  // apply basic cuts
  RooDataSet *dsRaw = dynamic_cast<RooDataSet *>(inputData.Get("dataset"));
  dsRaw->Print();
  if (!dsRaw)
  {
    cout << "[Error]: Cannot find RooDataSet\n";
    return;
  }

  RooDataSet *dsWeight = new RooDataSet("dsWeight", "dataset with weight", dsRaw, *dsRaw->get(), 0, "weight");
  dsWeight->Print();

  // define cut
  char reduceDS_woCtErr[3000];
  sprintf(reduceDS_woCtErr, // no multiplicty case
          "( ((abs(eta1) <= 1.2) && (pt1 >=3.5)) || ((abs(eta2) <= 1.2) && (pt2 >=3.5)) || ((abs(eta1) > 1.2) && (abs(eta1) <= 2.1) && (pt1 >= 5.47-1.89*(abs(eta1)))) || ((abs(eta2) > 1.2)  && (abs(eta2) <= 2.1) && (pt2 >= 5.47-1.89*(abs(eta2)))) || ((abs(eta1) > 2.1) && (abs(eta1) <= 2.4) && (pt1 >= 1.5)) || ((abs(eta2) > 2.1) && (abs(eta2) <= 2.4) && (pt2 >= 1.5)) )"

          "&& (pt >= %.3f && pt < %.3f && abs(y) >= %.3f && abs(y) < %.3f && mass >= %.3f && mass < %.3f && ctau3D >= %.3f && ctau3D < %.3f)"

          " && (recoQQsign == 0) ",
          ptLow, ptHigh, yLow, yHigh, massLow, massHigh, ctMin, ctMax);

  RooDataSet *dsReduced_tmp = (RooDataSet *)dsWeight->reduce(reduceDS_woCtErr);

  // observable
  auto mass = new RooRealVar("mass", "invariant mass", massLow, massHigh, "GeV/c^{2}");
  auto ctau3D = new RooRealVar("ctau3D", "", ctMin, ctMax, "#sigma_{L_{J/#psi}} [mm]");
  auto ctau3DErr = new RooRealVar("ctau3DErr", "", 0, 1, "#sigma_{L_{J/#psi}} [mm]");
  ctau3DErr->setBins(200);
  auto weight = new RooRealVar("weight", "", 0, 10000, "");
  RooArgSet obs(*mass, *ctau3DErr, *weight);
  auto dsReduced = new RooDataSet("dsReduced", "dataset with local vars", obs, Import(*dsReduced_tmp));
  cout << "\n--- Objects in reduced dataset ---\n";
  dsReduced->Print();

  RooDataSet *redDataSIG_tmp = (RooDataSet *)dsReduced->reduce("mass>2.9 && mass<3.3");
  RooDataSet *redDataSB_tmp = (RooDataSet *)dsReduced->reduce("mass<2.9 || mass>3.3");

  // === start getCtErrRange ===
  // get raw mass fit result -> 그러면 여기서 mass fit 할 필요 X
  RooRealVar sl1("sl1","", 0.0, -1.0, 1.0);
  RooRealVar sl2("sl2","", 0.0, -1.0, 1.0);

  {
    TFile fPrev(Form("%s/mass_fit_pT%.1f_%.1f.root", rootDir.Data(), ptLow, ptHigh), "READ");
    if (fPrev.IsZombie())
    {
      ::Error("cheb-load", "Cannot open previous mass-fit file.");
    }
    else
    {
      RooFitResult *fr = dynamic_cast<RooFitResult *>(fPrev.Get("fitMass"));
      if (!fr)
      {
        ::Error("cheb-load", "RooFitResult 'fitMass' not found in file.");
      }
      else
      {
        auto apply_one = [&](RooRealVar &dest)
        {
          const RooArgList &fl = fr->floatParsFinal();
          const RooArgList &cl = fr->constPars();
          if (auto *a = fl.find(dest.GetName()))
          {
            dest.setVal(static_cast<const RooRealVar *>(a)->getVal());
            dest.setConstant(true);
            return;
          }
          if (auto *a = cl.find(dest.GetName()))
          {
            dest.setVal(static_cast<const RooRealVar *>(a)->getVal());
            dest.setConstant(true);
            return;
          }
          ::Warning("cheb-load", "Parameter '%s' not found in fit result; keep as-is.", dest.GetName());
        };

        apply_one(sl1);
        apply_one(sl2);
        // apply_one(sl3);
      }
    }
  }

  RooArgList allPars_local;
  allPars_local.add(sl1);
  allPars_local.add(sl2);
  double scaleF = computeScaleF_cheby(allPars_local, 2, massLow, massHigh, sbL_lo, sbL_hi, sig_lo, sig_hi, sbR_lo, sbR_hi);
  cout << "scaleF: " << scaleF << "\n";

  // --- subtract sig - bkg ---
  // original RooDataHist
  RooDataHist *binDataCtErrSB = new RooDataHist("binDataCtErrSB", "Data ct error distribution for bkg", *ctau3DErr, *redDataSB_tmp);
  RooDataHist *binDataCtErrSIG = new RooDataHist("binDataCtErrSIG", "Data ct error distribution for sig", *ctau3DErr, *redDataSIG_tmp);

  RooDataHist* binSubtrSIG = new RooDataHist("binSubtrSIG","SIG-SB", RooArgSet(*ctau3DErr));
  RooDataHist* binScaleBKG = new RooDataHist("binScaleBKG","scaled SB", RooArgSet(*ctau3DErr));

  subtractSidebands(binSubtrSIG, binScaleBKG,
                                          binDataCtErrSIG, binDataCtErrSB,
                                          ctau3DErr,
                                          scaleF, "ctau3DErr",
                                          /*dk=*/0.0, /*floorW=*/1e-1);

  // --- find out proper bin ranges ---
  // RooDataHist -> TH1
  TH1 *histDataCtErrSIG = binDataCtErrSIG->createHistogram("histDataCtErrSIG", *ctau3DErr);
  TH1 *histSubtractedSIG = binSubtrSIG->createHistogram("histSubtractedSIG", *ctau3DErr);
  TH1 *histScaledBKG = binScaleBKG->createHistogram("histScaledBKG", *ctau3DErr);

  double minSig = 0.5, maxSig = 0.0, minBkg = 0.5, maxBkg = 0.0;
  double cutValue = 0.2;

  // maximum value bin
  int maxBinSig = histSubtractedSIG->GetMaximumBin();
  int maxBinBkg = histScaledBKG->GetMaximumBin();

  // initial limits
  minSig = histSubtractedSIG->GetBinLowEdge(maxBinSig);
  minBkg = histScaledBKG->GetBinLowEdge(maxBinBkg);

  maxSig = histSubtractedSIG->GetBinLowEdge(maxBinSig + 1);
  maxBkg = histScaledBKG->GetBinLowEdge(maxBinBkg + 1);

  // find out the min
  for (int xbins = maxBinSig; xbins > 0; --xbins)
  {
    if (histSubtractedSIG->GetBinContent(xbins) > cutValue)
    {
      minSig = histSubtractedSIG->GetBinLowEdge(xbins);
    }
    else
      break;
  }
  for (int xbins = maxBinBkg; xbins > 0; --xbins)
  {
    if (histScaledBKG->GetBinContent(xbins) > cutValue)
    {
      minBkg = histScaledBKG->GetBinLowEdge(xbins);
    }
    else
      break;
  }

  // find out max
  for (int xbins = maxBinSig; xbins <= histSubtractedSIG->GetNbinsX(); ++xbins)
  {
    if (histSubtractedSIG->GetBinContent(xbins) > cutValue)
    {
      maxSig = histSubtractedSIG->GetBinLowEdge(xbins + 1);
    }
    else
      break;
  }
  for (int xbins = maxBinBkg; xbins <= histScaledBKG->GetNbinsX(); ++xbins)
  { // ★ fix: maxBinBkg
    if (histScaledBKG->GetBinContent(xbins) > cutValue)
    {
      maxBkg = histScaledBKG->GetBinLowEdge(xbins + 1);
    }
    else
      break;
  }

  double tmpMin = (minSig > minBkg) ? minSig : minBkg;
  double tmpMax = (maxSig < maxBkg) ? maxSig : maxBkg;

  // rounding
  tmpMin = TMath::Floor(tmpMin * 1000.0) / 1000.0;
  tmpMax = TMath::Ceil(tmpMax * 1000.0) / 1000.0;

  // reduce dataset
  char reduceDS[512];
  std::snprintf(reduceDS, sizeof(reduceDS), "ctau3DErr > %.3f && ctau3DErr < %.3f", tmpMin, tmpMax);
  RooDataSet *redDataTmp = (RooDataSet *)dsReduced->reduce(reduceDS);

  if (redDataTmp->sumEntries() < dsReduced->sumEntries() * 0.9)
  {
    delete redDataTmp;
    std::snprintf(reduceDS, sizeof(reduceDS), "ctau3DErr > %.3f && ctau3DErr < %.3f", minSig, maxSig);
    redDataTmp = (RooDataSet *)dsReduced->reduce(reduceDS);
    tmpMin = minSig;
    tmpMax = maxSig;
  }

  // leaset range
  if ((tmpMax - tmpMin) < 0.008)
  {
    std::cout << "getCtErrRange:: Maximum is less than minimum! Possibly few events in this bin.\n";
    tmpMax = tmpMin + 0.008;
  }

  // --- draw raw plots ---
  TCanvas c0("ctau_err", "ctau_err", 800, 800);
  c0.cd();
  c0.SetLogy(1);

  RooPlot *errFrame = ctau3DErr->frame();
  errFrame->SetTitle("");
  errFrame->GetXaxis()->SetTitle("ctau3DErr [mm]");
  errFrame->GetXaxis()->CenterTitle();
  errFrame->GetYaxis()->SetTitle("Counts");
  double maxY = std::max({histDataCtErrSIG->GetMaximum(),
                          histScaledBKG->GetMaximum(),
                          histSubtractedSIG->GetMaximum()});
  errFrame->SetMaximum(maxY * 1.3);
  errFrame->SetMinimum(0.2);
  errFrame->Draw();

  histDataCtErrSIG->SetMarkerColor(kRed);
  histDataCtErrSIG->SetLineColor(kWhite);
  histDataCtErrSIG->SetMarkerStyle(24);
  histDataCtErrSIG->GetXaxis()->CenterTitle(1);
  histDataCtErrSIG->GetYaxis()->CenterTitle(1);
  histDataCtErrSIG->Draw("pe same");

  histScaledBKG->SetMarkerColor(kBlue);
  histScaledBKG->SetLineColor(kWhite);
  histScaledBKG->SetMarkerStyle(24);
  histScaledBKG->Draw("pe same");

  histSubtractedSIG->SetLineColor(kWhite);
  histSubtractedSIG->SetMarkerSize(0.3);
  histSubtractedSIG->Draw("pe same");

  TLatex t;
  t.SetNDC();
  t.SetTextAlign(32);
  t.SetTextSize(0.04);
  t.SetTextColor(kRed);
  char comment[200];
  std::snprintf(comment, sizeof(comment), "Range: %.3f-%.3f (mm)", tmpMin, tmpMax);
  t.DrawLatex(0.92, 0.60, comment);
  t.SetTextColor(kBlack);

  TLegend legsb(0.60, 0.39, 0.90, 0.55, nullptr, "brNDC");
  legsb.SetFillStyle(0);
  legsb.SetBorderSize(0);
  legsb.SetShadowColor(0);
  legsb.SetMargin(0.2);
  legsb.SetTextFont(42);
  legsb.SetTextSize(0.035);
  legsb.AddEntry(histDataCtErrSIG, "sig cands", "p");
  legsb.AddEntry(histScaledBKG, "scaled bkg", "p");
  legsb.AddEntry(histSubtractedSIG, "sig (= cands - bkg)", "p");
  legsb.Draw("same");


  c0.SaveAs(Form("%s/ctau_err_raw_pT%.1f_%.1f.png", figDir.Data(), ptLow, ptHigh));

  errmin = tmpMin;
  errmax = tmpMax;

  // === finish getCtErrRange ===


  // === start building dataset and punzi-term with ctau3DErr cut ===
  // --- make RooDataSet for the next processes ---
  ctau3DErr = new RooRealVar("ctau3DErr", "", errmin, errmax, "#sigma_{L_{J/#psi}} [mm]");
  // ctau3DErr->setBins(50);
  RooArgSet obs2(*mass, *ctau3D, *ctau3DErr, *weight);
  auto dsReduced2 = new RooDataSet("dsReduced2", "dataset with local vars", obs2, Import(*dsReduced_tmp));;

  char reduceDS_wCtErr[3000];
  sprintf(reduceDS_wCtErr, // no multiplicty case
        "ctau3DErr > %.3f && ctau3DErr < %.3f", errmin, errmax);
  auto dsReduced3 = (RooDataSet *)dsReduced2->reduce(reduceDS_wCtErr);
  RooDataSet *redDataSB = (RooDataSet *)dsReduced3->reduce("mass<2.9 || mass>3.3");

  // --- build RooHistFunc---
  RooArgList obs_cterr(*ctau3DErr);
  RooHistPdf sigPdf("sigPdf", "Signal hist func",
                    obs_cterr, *binSubtrSIG, 0);

  RooHistPdf bkgPdf("bkgPdf", "Background hist func",
                    obs_cterr, *binScaleBKG, 0);


  // --- draw RooHistFunc plots ---
  auto findObj = [&](RooPlot *fr, const char *n) -> TObject * { return fr ? fr->findObject(n) : nullptr; };
  // --- Overlay + Pull: SIG ---
  {
    TCanvas c("c_mass", "c_mass", 800, 800);
    TPad pad1("pad1", "pad1", 0.0, 0.25, 1.0, 1.0);
    pad1.SetBottomMargin(0.00001);
    pad1.SetLogy();
    pad1.Draw();
    pad1.cd();

    // --- frame & plot ---
    RooPlot *frSig = ctau3DErr->frame();
    binSubtrSIG->plotOn(frSig, DataError(RooAbsData::SumW2), Name("data"));
    sigPdf.plotOn(frSig, Name("model"), LineColor(kRed));

    // --- dynamic y-range for log scale ---
    double ymin = 1e300, ymax = -1e300;
    if (auto *hdata = dynamic_cast<RooHist *>(frSig->getHist("data")))
    {
      for (int i = 0; i < hdata->GetN(); ++i)
      {
        double x, y;
        hdata->GetPoint(i, x, y);
        if (y > 0 && y < ymin)
          ymin = y;
        if (y > ymax)
          ymax = y;
      }
    }
    if (ymin <= 0 || ymin == 1e300)
      ymin = 1e-3;
    frSig->SetMinimum(ymin * 0.5);
    frSig->SetMaximum(std::max(ymax, ymin) * 1e4);

    frSig->GetYaxis()->SetTitle("Events");
    frSig->GetXaxis()->SetTitle("");
    frSig->Draw("e");

    // --- legend ---
    TLegend leg(0.49, 0.65, 0.70, 0.94);
    {
      leg.SetBorderSize(0);
      leg.SetFillStyle(0);
      leg.SetTextSize(0.03);
      if (auto *o = findObj(frSig, "data"))
        leg.AddEntry(o, "Data", "lep");
      if (auto *o = findObj(frSig, "model"))
        leg.AddEntry(o, "Sig", "pe");

      leg.Draw("same");
    }

    // --- CMS/info latex ---
    {
      TLatex tx;
      tx.SetNDC();
      tx.SetTextSize(0.03);
      tx.SetTextFont(42);
      double x = 0.19, y0 = 0.90, dy = -0.06;
      int k = 0;
      tx.DrawLatex(x, y0 + dy * k++, "CMS pp Ref. #sqrt{s_{NN}} = 5.02 TeV");
      tx.DrawLatex(x, y0 + dy * k++, "Data, J/#psi #rightarrow #mu^{+}#mu^{-}");
      if (yLow == 0)
        tx.DrawLatex(x, y0 + dy * k++, Form("%.1f < p_{T} < %.1f, |y| < %.1f", ptLow, ptHigh, yHigh));
      else
        tx.DrawLatex(x, y0 + dy * k++, Form("%.1f < p_{T} < %.1f, %.1f < y < %.1f", ptLow, ptHigh, yLow, yHigh));
      // tx.DrawLatex(x, y0 + dy * k++, Form("%.3f < c#tau < %.3f", ctMin, ctMax));
    }

    // --- pull pad ---
    c.cd();
    TPad pad2("pad2", "pad2", 0.0, 0.0, 1.0, 0.25);
    pad2.SetTopMargin(0.00001);
    pad2.SetBottomMargin(0.4);
    pad2.Draw();
    pad2.cd();

    RooHist *hpull = frSig->pullHist("data", "model");
    RooPlot *fpull = ctau3DErr->frame(Title(""));
    fpull->addPlotable(hpull, "P");
    fpull->GetYaxis()->SetTitle("Pull");
    fpull->GetXaxis()->SetTitle("ctau3DErr [mm]");
    fpull->GetXaxis()->CenterTitle();
    fpull->SetMinimum(-8);
    fpull->SetMaximum(8);
    fpull->GetYaxis()->SetNdivisions(505);
    fpull->GetYaxis()->SetTitleSize(0.12);
    fpull->GetYaxis()->SetLabelSize(0.10);
    fpull->GetXaxis()->SetTitleSize(0.15);
    fpull->GetXaxis()->SetLabelSize(0.10);
    fpull->Draw();

    TLine line(errmin, 0.0, errmax, 0.0);
    line.SetLineStyle(2);
    line.Draw("same");

    c.SaveAs(Form("%s/ctau_err_sig_pT%.1f_%.1f.png", figDir.Data(), ptLow, ptHigh));
  }

  // --- BKG ---
  {
    TCanvas c("c_mass", "c_mass", 800, 800);
    TPad pad1("pad1", "pad1", 0.0, 0.25, 1.0, 1.0);
    pad1.SetBottomMargin(0.00001);
    pad1.SetLogy();
    pad1.Draw();
    pad1.cd();

    // --- frame & plot ---
    RooPlot *frBkg = ctau3DErr->frame();
    binScaleBKG->plotOn(frBkg, DataError(RooAbsData::SumW2), Name("data"));
    bkgPdf.plotOn(frBkg, Name("model"));

    // --- dynamic y-range for log scale ---
    double ymin = 1e300, ymax = -1e300;
    if (auto *hdata = dynamic_cast<RooHist *>(frBkg->getHist("data")))
    {
      for (int i = 0; i < hdata->GetN(); ++i)
      {
        double x, y;
        hdata->GetPoint(i, x, y);
        if (y > 0 && y < ymin)
          ymin = y;
        if (y > ymax)
          ymax = y;
      }
    }
    if (ymin <= 0 || ymin == 1e300)
      ymin = 1e-3;
    frBkg->SetMinimum(ymin * 0.5);
    frBkg->SetMaximum(std::max(ymax, ymin) * 1e4);

    frBkg->GetYaxis()->SetTitle("Events");
    frBkg->GetXaxis()->SetTitle("");
    frBkg->Draw("e");

    // --- legend ---
    TLegend leg(0.49, 0.65, 0.70, 0.94);
    {
      leg.SetBorderSize(0);
      leg.SetFillStyle(0);
      leg.SetTextSize(0.03);
      if (auto *o = findObj(frBkg, "data"))
        leg.AddEntry(o, "Data", "lep");
      if (auto *o = findObj(frBkg, "model"))
        leg.AddEntry(o, "Bkg", "pe");

      leg.Draw("same");
    }

    // --- CMS/info latex ---
    {
      TLatex tx;
      tx.SetNDC();
      tx.SetTextSize(0.03);
      tx.SetTextFont(42);
      double x = 0.19, y0 = 0.90, dy = -0.06;
      int k = 0;
      tx.DrawLatex(x, y0 + dy * k++, "CMS pp Ref. #sqrt{s_{NN}} = 5.02 TeV");
      tx.DrawLatex(x, y0 + dy * k++, "Data, J/#psi #rightarrow #mu^{+}#mu^{-}");
      if (yLow == 0)
        tx.DrawLatex(x, y0 + dy * k++, Form("%.1f < p_{T} < %.1f, |y| < %.1f", ptLow, ptHigh, yHigh));
      else
        tx.DrawLatex(x, y0 + dy * k++, Form("%.1f < p_{T} < %.1f, %.1f < y < %.1f", ptLow, ptHigh, yLow, yHigh));
      // tx.DrawLatex(x, y0 + dy * k++, Form("%.3f < c#tau < %.3f", ctMin, ctMax));
    }
    
    // --- pull pad ---
    c.cd();
    TPad pad2("pad2", "pad2", 0.0, 0.0, 1.0, 0.25);
    pad2.SetTopMargin(0.00001);
    pad2.SetBottomMargin(0.4);
    pad2.Draw();
    pad2.cd();

    RooHist *hpull = frBkg->pullHist("data", "model");
    RooPlot *fpull = ctau3DErr->frame(Title(""));
    fpull->addPlotable(hpull, "P");
    fpull->GetYaxis()->SetTitle("Pull");
    fpull->GetXaxis()->SetTitle("ctau3DErr [mm]");
    fpull->GetXaxis()->CenterTitle();
    fpull->SetMinimum(-8);
    fpull->SetMaximum(8);
    fpull->GetYaxis()->SetNdivisions(505);
    fpull->GetYaxis()->SetTitleSize(0.12);
    fpull->GetYaxis()->SetLabelSize(0.10);
    fpull->GetXaxis()->SetTitleSize(0.15);
    fpull->GetXaxis()->SetLabelSize(0.10);
    fpull->Draw();

    TLine line(errmin, 0.0, errmax, 0.0);
    line.SetLineStyle(2);
    line.Draw("same");

    c.SaveAs(Form("%s/ctau_err_bkg_pT%.1f_%.1f.png", figDir.Data(), ptLow, ptHigh));
  }

  // === Finish building dataset and punzi-term with ctau3DErr cut ===

  // save result
  // --- save observable ranges as RooValues ---
  RooRealVar rooCtLow("rooCtLow", "", ctMin);
  RooRealVar rooCtHigh("rooCtHigh", "", ctMax);
  RooRealVar rooCtErrLow("rooCtErrLow", "", errmin);
  RooRealVar rooCtErrHigh("rooCtErrHigh", "", errmax);

  TFile output(Form("%s/ctau_err_pT%.1f_%.1f.root", rootDir.Data(), ptLow, ptHigh), "recreate");
  rooCtLow.Write();
  rooCtHigh.Write();
  rooCtErrLow.Write();
  rooCtErrHigh.Write();
  sigPdf.Write();
  bkgPdf.Write();
  redDataSB->Write("redDataSB");
  dsReduced3->Write("dsReduced");
  output.Close();
  // mass range, ctau range, ctauErr range

  cout << "=== finish sideband() ===\n";
  time.Stop();
  printf("RealTime=%f seconds, CpuTime=%f seconds\n", time.RealTime(), time.CpuTime());
}