#include <iostream>
#include <iomanip>
#include <cmath>
#include <RooWorkspace.h>
#include <RooDataSet.h>
#include <RooRealVar.h>
#include <RooArgSet.h>
#include <TFile.h>
#include <TChain.h>
#include <TBranch.h>
#include <TStopwatch.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TMath.h>
#include <RooPlot.h>

static const long MAXTREESIZE = 1000000000000;

// ------------------------------------------------------------------
// acceptance
// ------------------------------------------------------------------
// from Run2 cutsAndBins.h
// bool IsAcceptanceQQ(double pt, double eta)
// {
//   return ((fabs(eta) < 1.2 && pt >= 3.5) ||
//           (1.2 <= fabs(eta) && fabs(eta) < 2.1 && pt >= 5.47 - 1.89 * fabs(eta)) ||
//           (2.1 <= fabs(eta) && fabs(eta) < 2.4 && pt >= 1.5));
// }

bool IsAcceptanceQQ(double pt, double eta)
{
    const double aeta = std::abs(eta);

  if (aeta > 2.4)
    return false;
  if (aeta <= 1.0)
    return pt >= 3.3;
  if (aeta <= 1.3)
    return pt >= ((2.1 - 3.3) / (1.3 - 1.0)) * (aeta - 1.0) + 3.3;
  if (aeta <= 1.7)
    return pt >= ((1.0 - 2.1) / (1.7 - 1.3)) * (aeta - 1.3) + 2.1;

  return pt >= 1.0;

  // 2018 soft muon cut style - for Run3
  // https://twiki.cern.ch/twiki/bin/view/CMS/DileptonMuonSelection
  // return ((fabs(eta) < 0.3 && pt > 3.4) ||
  //         (0.3 < fabs(eta) && fabs(eta) < 1.1 && pt > 3.3) ||
  //         (1.1 < fabs(eta) && fabs(eta) < 1.5 && pt > 9.08 - 5.25*fabs(eta)) ||
  //         (1.5 < fabs(eta) && fabs(eta) < 2.4 && pt > 0.8 && pt>2.4-0.8*fabs(eta)));
}

inline void EnableBranchStatus(TChain *ch, bool isMC)
{
  ch->SetBranchStatus("*", 0);

  auto enableIf = [&](const char *name) {
    if (ch->GetBranch(name)) ch->SetBranchStatus(name, 1);
  };

  // event variables
  enableIf("eventNb");
  // ch->SetBranchStatus("cBin", 1);
  enableIf("nDimu");
  enableIf("zVtx");
  enableIf("weight");

  // dimuon
  enableIf("recoQQsign");
  enableIf("mass");
  enableIf("pt");
  enableIf("y");

  // single muon
  enableIf("pt1");
  enableIf("pt2");
  enableIf("eta1");
  enableIf("eta2");

  // dimuon lifetime
  enableIf("ctau3D");
  enableIf("ctau3DErr");
  enableIf("ctau3DRes");
  // ch->SetBranchStatus("ctau", 1);
  // ch->SetBranchStatus("ctauErr", 1);
  // ch->SetBranchStatus("ctauRes", 1);
  if (isMC) enableIf("ctau3DTrue");

  // polarization
  // ch->SetBranchStatus("cos_theta", 1);
  // ch->SetBranchStatus("cos_theta1", 1);
  enableIf("phi");
  enableIf("phi1");
  // ch->SetBranchStatus("recoQQdca", 1);

  enableIf("cosHX");
  enableIf("phiHX");
  enableIf("cosCS");
  enableIf("phiCS");
  // ch->SetBranchStatus("cos_ep", 1);
  // ch->SetBranchStatus("phi_ep", 1);
}


// ------------------------------------------------------------------
// main macro
// ------------------------------------------------------------------
void skim_to_roodataset(
    std::string dataLabel = "_OO25_HLT_OxyL1SingleMuOpen_v1",
    int nEvt = -1,
    int cLow = 0, int cHigh = 200,
    float massLow = 2.6, float massHigh = 3.5,
    bool oppositeSign = true, // true=OS, false=SS
    bool isMC = true, int MCtype = 2,
    bool fAccW = false, bool fEffW = false,
    bool isTnP = false, bool isPtW = false,
    int hiHFBinEdge = 0)
{
  using namespace RooFit;
  TStopwatch t;
  t.Start();
  std::cout << "\n=== Convert FlowSkim to RooDataSet ===\n";

  // ------------------------------------------------------------------
  // style
  // ------------------------------------------------------------------
  gStyle->SetOptStat(0);
  gROOT->ForceStyle();

  TString DATE = dataLabel.c_str();

  // dimuon label - MUST use the opposite sign
  TString oppositeSignString = oppositeSign ? "OS" : "SS";

  // MC label - Jpsi MC must be PR or NP (CMS doesn't have the Inclusive Jpsi MC)
  TString MCtype_tag = (MCtype == 1) ? "PR" : (MCtype == 2 ? "NP" : "");

  // ------------------------------------------------------------------
  // input
  // ------------------------------------------------------------------
  // FlowSkim input
  TChain *muon_chain = new TChain("myTree");
  if (!isMC)
    muon_chain->Add("skim_roots/skim_OO2025_isMC0_Dimuon_MiniAOD_PromptReco_v1_Jul12_merged.root"); // data
  else {
    if (MCtype==1) muon_chain->Add("skim_roots/skim_OO2025_isMC1_PR_Dimuon_MiniAOD_PromptReco_v1_Jul12_merged.root"); // mc PR
    else if (MCtype==2) muon_chain->Add("skim_roots/skim_OO2025_isMC1_NP_Dimuon_MiniAOD_PromptReco_v1_Jul12_merged.root"); // mc NP
    else muon_chain->Add("skim_roots/skim_OO2025_isMC1_Dimuon_MiniAOD_PromptReco_v1_Jul12_merged.root"); // fallback
  }
  
  
  // ------------------------------------------------------------------
  // Acc×Eff correction inputs
  // ------------------------------------------------------------------
  // TH2D *h_correct[6] = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};
  // TFile *fin_correct = nullptr;

  // should distinguish - Acc
  // if (fAccW)
  // {}
  // Acc doesn't have a angle correction -> pT only

  // Eff
  // if (fEffW)
  // {
  //   fin_correct = TFile::Open("../eff_acc/roots/mc_eff_vs_pt_cent_rap_prompt_pbpb_Jpsi_PtW1_tnp1_250221.root", "READ");
  //   if (fin_correct && !fin_correct->IsZombie())
  //   {
  //     h_correct[0] = (TH2D *)fin_correct->Get("mc_eff_vs_pt_TnP1_PtW1_cent0to20_absy0_1p6");
  //     h_correct[1] = (TH2D *)fin_correct->Get("mc_eff_vs_pt_TnP1_PtW1_cent0to20_absy1p6_2p4");
  //     h_correct[2] = (TH2D *)fin_correct->Get("mc_eff_vs_pt_TnP1_PtW1_cent20_60_absy0_1p6");
  //     h_correct[3] = (TH2D *)fin_correct->Get("mc_eff_vs_pt_TnP1_PtW1_cent20_60_absy1p6_2p4");
  //     h_correct[4] = (TH2D *)fin_correct->Get("mc_eff_vs_pt_TnP1_PtW1_cent60to180_absy0_1p6");
  //     h_correct[5] = (TH2D *)fin_correct->Get("mc_eff_vs_pt_TnP1_PtW1_cent60to180_absy1p6_2p4");
  //   }
  //   else
  //   {
  //     std::cout << "[info] correction file not found. continue with weight=1.\n";
  //   }
  // }

  // ------------------------------------------------------------------
  // bind input branch
  // ------------------------------------------------------------------
  Int_t eventNb, cBin, nDimu;
  const int nMaxDimu = 1000;
  int recoQQsign[nMaxDimu];
  float zVtx, mass[nMaxDimu], pt[nMaxDimu], y[nMaxDimu];
  float phi[nMaxDimu], phi1[nMaxDimu];
  float ctau3D[nMaxDimu], ctau3DErr[nMaxDimu], ctau3DRes[nMaxDimu];
  float ctau[nMaxDimu], ctauErr[nMaxDimu], ctauRes[nMaxDimu];
  float ctau3DTrue[nMaxDimu] = {0};
  double weight;
  float recoQQdca[nMaxDimu];

  // float cos_theta[nMaxDimu], cos_theta1[nMaxDimu], phi[nMaxDimu], phi1[nMaxDimu];
  float cosHX[nMaxDimu], phiHX[nMaxDimu];
  float cosCS[nMaxDimu], phiCS[nMaxDimu];
  // float cos_ep[nMaxDimu], phi_ep[nMaxDimu];
  
  // branch optimization
  EnableBranchStatus(muon_chain, isMC); // turn on the branches you need

  auto hasBranch = [&](const char *name) -> bool {
    return muon_chain->GetBranch(name) != nullptr;
  };

  const bool hasCtau3DTrue = isMC && hasBranch("ctau3DTrue");
  const bool hasPhi = hasBranch("phi");
  const bool hasPhi1 = hasBranch("phi1");
  const bool hasCosHX = hasBranch("cosHX");
  const bool hasPhiHX = hasBranch("phiHX");
  const bool hasCosCS = hasBranch("cosCS");
  const bool hasPhiCS = hasBranch("phiCS");

  // bind branches (guard against missing branches)
  if (hasBranch("eventNb")) muon_chain->SetBranchAddress("eventNb", &eventNb);
  // muon_chain->SetBranchAddress("cBin", &cBin);
  if (hasBranch("nDimu")) muon_chain->SetBranchAddress("nDimu", &nDimu);
  if (hasBranch("recoQQsign")) muon_chain->SetBranchAddress("recoQQsign", recoQQsign);
  if (hasBranch("zVtx")) muon_chain->SetBranchAddress("zVtx", &zVtx);
  if (hasBranch("mass")) muon_chain->SetBranchAddress("mass", mass);
  if (hasBranch("y")) muon_chain->SetBranchAddress("y", y);
  if (hasBranch("pt")) muon_chain->SetBranchAddress("pt", pt);
  // skip pt1/pt2/eta/eta1/eta2 to reduce I/O (not in argSet)
  if (hasBranch("ctau3D")) muon_chain->SetBranchAddress("ctau3D", ctau3D);
  if (hasCtau3DTrue) muon_chain->SetBranchAddress("ctau3DTrue", ctau3DTrue);
  if (hasBranch("ctau3DErr")) muon_chain->SetBranchAddress("ctau3DErr", ctau3DErr);
  if (hasBranch("ctau3DRes")) muon_chain->SetBranchAddress("ctau3DRes", ctau3DRes);
  // muon_chain->SetBranchAddress("ctau", ctau);
  // muon_chain->SetBranchAddress("ctauErr", ctauErr);
  // muon_chain->SetBranchAddress("ctauRes", ctauRes);
  if (hasBranch("weight")) muon_chain->SetBranchAddress("weight", &weight);

  // muon_chain->SetBranchAddress("cos_theta", cos_theta);
  // muon_chain->SetBranchAddress("cos_theta1", cos_theta1);
  if (hasPhi) muon_chain->SetBranchAddress("phi", phi);
  if (hasPhi1) muon_chain->SetBranchAddress("phi1", phi1);
  // muon_chain->SetBranchAddress("recoQQdca", recoQQdca);
  if (hasCosHX) muon_chain->SetBranchAddress("cosHX", cosHX);
  if (hasPhiHX) muon_chain->SetBranchAddress("phiHX", phiHX);
  if (hasCosCS) muon_chain->SetBranchAddress("cosCS", cosCS);
  if (hasPhiCS) muon_chain->SetBranchAddress("phiCS", phiCS);
  // muon_chain->SetBranchAddress("cos_ep", cos_ep);
  // muon_chain->SetBranchAddress("phi_ep", phi_ep);

  // ------------------------------------------------------------------
  // output file
  // ------------------------------------------------------------------
  TFile *out_file = nullptr;
  if (isMC)
  {
    out_file = new TFile(Form("roodataset_roots/RooDataSet_miniAOD_isMC%d_%s_Jpsi_cent%i_%i_Effw%d_Accw%d_PtW%d_TnP%d%s.root",
                              isMC, MCtype_tag.Data(), cLow, cHigh, fEffW, fAccW, isPtW, isTnP, DATE.Data()),
                         "RECREATE");
  }
  else
  {
    out_file = new TFile(Form("roodataset_roots/RooDataSet_miniAOD_isMC%d_Jpsi_cent%i_%i_Effw%d_Accw%d_PtW%d_TnP%d%s.root",
                              isMC, cLow, cHigh, fEffW, fAccW, isPtW, isTnP, DATE.Data()),
                         "RECREATE");
  }

  // ------------------------------------------------------------------
  // Roo variables
  // ------------------------------------------------------------------
  RooRealVar massVar("mass", "mass", 2.6, 3.5, "GeV/c^{2}");
  RooRealVar ptVar("pt", "pt", 0, 100, "GeV/c");
  RooRealVar yVar("y", "rapidity", -20, 20);
  RooRealVar pt1Var("pt1", "pt of muon+", 0, 500, "GeV/c");
  RooRealVar eta1Var("eta1", "eta of muon+", -4, 4);
  RooRealVar pt2Var("pt2", "pt of muon-", 0, 500, "GeV/c");
  RooRealVar eta2Var("eta2", "eta of muon-", -4, 4);
  // RooRealVar cBinVar("cBin", "Centrality", 0, 200);
  RooRealVar evtWeight("weight", "corr weight", 0, 10000);
  RooRealVar recoQQ("recoQQsign", "qq sign", -1, 3);
  RooRealVar ctau3DVar("ctau3D", "c_{#tau}", -10.0, 10.0, "mm");
  RooRealVar ctau3DErrVar("ctau3DErr", "#sigma_{c#tau}", 0, 10.0, "mm");
  RooRealVar ctau3DResVar("ctau3DRes", "c_{#tau} res", -20.0, 20.0);
  // RooRealVar ctauVar("ctau", "c_{#tau}", -10.0, 10.0, "mm");
  // RooRealVar ctauErrVar("ctauErr", "#sigma_{c#tau}", 0, 10.0, "mm");
  // RooRealVar ctauResVar("ctauRes", "c_{#tau} res", -20.0, 20.0);
  RooRealVar ctau3DTrueVar("ctau3DTrue", "c_{#tau} true", -0.5, 10.0, "mm");
  RooRealVar NumDimu("NumDimu", "number of dimuon", 0, 1000);

  // RooRealVar cos_thetaVar("cos_theta", "", -1.0, 1.0);
  // RooRealVar cos_theta1Var("cos_theta1", "", -1.0, 1.0);
  RooRealVar phiVar("phi", "", -5.0, 5.0);
  RooRealVar phi1Var("phi1", "", -5.0, 5.0);
  // RooRealVar recoQQdcaVar("recoQQdca", "", 0, 5.0);
  RooRealVar cosHXVar("cosHX", "", -1.0, 1.0);
  RooRealVar phiHXVar("phiHX", "", -1.0, 1.0);
  RooRealVar cosCSVar("cosCS", "", -5.0, 5.0);
  RooRealVar phiCSVar("phiCS", "", -5.0, 5.0);
  // RooRealVar cos_epVar("cos_ep", "", -1.0, 1.0);
  // RooRealVar phi_epVar("phi_ep", "", -5.0, 5.0);

  RooArgSet argSet(massVar, ptVar, yVar, evtWeight, ctau3DVar, ctau3DErrVar, ctau3DResVar);
  // pt1Var, pt2Var, eta1Var, eta2Var, phiVar, phi1Var, recoQQ, NumDimu
  
  // argSet.add(ctauVar);
  // argSet.add(ctauErrVar);
  // argSet.add(ctauResVar);
  // argSet.add(recoQQdcaVar);
  
  if (hasCosHX) argSet.add(cosHXVar);
  if (hasPhiHX) argSet.add(phiHXVar);
  if (hasCosCS) argSet.add(cosCSVar);
  if (hasPhiCS) argSet.add(phiCSVar);

  // RooArgSet argSet(massVar, ptVar, yVar, pt1Var, pt2Var, eta1Var, eta2Var, evtWeight,
  //                  cBinVar, recoQQ, NumDimu, ctau3DVar, ctau3DErrVar, ctau3DResVar,
  //                  cos_thetaVar, cos_theta1Var, phiVar, phi1Var, cosHXVar, phiHXVar,
  //                  cosCSVar, phiCSVar, cos_epVar, phi_epVar);
  if (hasCtau3DTrue) argSet.add(ctau3DTrueVar);

  // ------------------------------------------------------------------
  // RooDataSet
  // ------------------------------------------------------------------
  // declare RooDataSet
  RooDataSet ds("dataset", "", argSet);

  // ------------------------------------------------------------------
  // counters / bookkeeping
  // ------------------------------------------------------------------
  // error counter
  int nZeroCtauErr = 0;
  int nNonFiniteCtau = 0;
  int nZeroOrNegWeight = 0;
  int nNullHistoAccess = 0;

  // set event number
  if (nEvt == -1)
    nEvt = muon_chain->GetEntries();
  std::cout << "Total events = " << nEvt << "\n";

  // print progress
  const Long64_t totalEvents = nEvt;
  const Long64_t reportEvery =
      (totalEvents >= 10'000'000) ? 200'000 :
      (totalEvents >= 1'000'000) ? 100'000 :
      (totalEvents >= 100'000) ? 20'000 : 10'000;

  // number of dimuon count
  int count_dimu = 0;

  // OS trigger -> Use OS
  const int signWanted = oppositeSign ? 0 : 1; // 0=OS, 1=SS

  // ------------------------------------------------------------------
  // main loop
  // ------------------------------------------------------------------
  for (int i = 0; i < nEvt; ++i)
  {
    muon_chain->GetEntry(i);
    
    // print progress and ETA
    if ((i % reportEvery) == 0 || (i + 1) == nEvt)
    {
      const double elapsed = t.RealTime(); // accumulated running time (s)
      t.Continue();                        // keep accumluation
      const double frac = (totalEvents > 0) ? double(i + 1) / double(totalEvents) : 0.0;
      const double estTotal = (frac > 0) ? (elapsed / frac) : 0.0;
      const double eta = estTotal - elapsed; // ETA (s)

      // formatting
      auto to_hms = [](double sec)
      {
        int h = int(sec / 3600);
        sec -= 3600 * h;
        int m = int(sec / 60);
        sec -= 60 * m;
        int s = int(sec + 0.5);
        std::ostringstream os;
        os << std::setfill('0') << std::setw(2) << h << ":" << std::setw(2) << m << ":" << std::setw(2) << s;
        return os.str();
      };

      std::cout << "["
                << std::fixed << std::setprecision(1)
                << (frac * 100.0) << "%] "
                << (i + 1) << " / " << totalEvents
                << "  | elapsed " << to_hms(elapsed)
                << "  | ETA " << to_hms(std::max(0.0, eta))
                << std::endl;
    }

    // if (TMath::Abs(zVtx) > 15)
    //   continue;
    
    // if (!(cBin >= cLow && cBin < cHigh))
    //   continue;

    // ===== dimuon loop =====
    const int nDimuSafe = std::min(nDimu, nMaxDimu);
    for (int j = 0; j < nDimuSafe; ++j)
    {
      // basic cut: pt<50, sign, mass, |y|<2.4, muon acc
      // applied in the onia_to_skim.C
      // if (!(pt[j] < 50 &&
      //       recoQQsign[j] == signWanted &&
      //       mass[j] > massLow && mass[j] < massHigh &&
      //       TMath::Abs(y[j]) < 2.4 &&
      //       IsAcceptanceQQ(pt1[j], eta1[j]) &&
      //       IsAcceptanceQQ(pt2[j], eta2[j])))
      //   continue;

      // === acc×eff weight (default: 1) ===
      // double weight_acc_eff = 1.0;
      // // 나중에 할 거 - Acc, Eff 따로 켜고 끄기
      // if (fAccW || fEffW)
      // {
      //   const double absy = TMath::Abs(y[j]);
      //   // if (cBin < 20)
      //   // {
      //   //   weight_acc_eff = (absy < 1.6) ? getCorrectionSafe(h_correct[0], pt[j], cos_ep[j])
      //   //                                 : (absy < 2.4 ? getCorrectionSafe(h_correct[1], pt[j], cos_ep[j]) : 1.0);
      //   // }
      //   // else if (cBin >= 20 && cBin < 60)
      //   // {
      //   //   weight_acc_eff = (absy < 1.6) ? getCorrectionSafe(h_correct[2], pt[j], cos_ep[j])
      //   //                                 : (absy < 2.4 ? getCorrectionSafe(h_correct[3], pt[j], cos_ep[j]) : 1.0);
      //   // }
      //   // else if (cBin >= 60 && cBin < 180)
      //   // {
      //   //   weight_acc_eff = (absy < 1.6) ? getCorrectionSafe(h_correct[4], pt[j], cos_ep[j])
      //   //                                 : (absy < 2.4 ? getCorrectionSafe(h_correct[5], pt[j], cos_ep[j]) : 1.0);
      //   // }

      //   weight_acc_eff = 1;
      // }
      // const double weight_final = weight * weight_acc_eff; // 여기도 분리 weight * weight_acc * weight eff

      // // correction error counter 
      // if (!h_correct[0] && (fAccW || fEffW))
      //   nNullHistoAccess++; // can't read correction histogram
      // if (weight_acc_eff <= 0)
      //   nZeroOrNegWeight++; // correction < 0

      // ----- push values into RooVariables -----
      recoQQ.setVal(recoQQsign[j]);
      massVar.setVal(mass[j]);
      ptVar.setVal(pt[j]);
      yVar.setVal(y[j]);
      // pt1/pt2/eta1/eta2 are not stored in argSet here
      // cBinVar.setVal(cBin);
      const double ct = static_cast<double>(ctau3D[j]);
      const double ctErr = static_cast<double>(ctau3DErr[j]);
      if (!std::isfinite(ct) || !std::isfinite(ctErr))
      {
        ++nNonFiniteCtau;
        continue;
      }
      ctau3DVar.setVal(ct);
      ctau3DErrVar.setVal(ctErr);
      
      // prevent 0-division
      double err = ctErr;
      if (err == 0)
        nZeroCtauErr++;
      ctau3DResVar.setVal((err != 0.0) ? (ct / err) : 0.0);
      if (hasCtau3DTrue) ctau3DTrueVar.setVal(ctau3DTrue[j]);


      // ctauVar.setVal(ctau[j]);
      // ctauErrVar.setVal(ctauErr[j]);
      
      // // prevent 0-division
      // err = (double)ctauErr[j];
      // if (err == 0)
      //   nZeroCtauErr++;
      // ctauResVar.setVal((err != 0.0) ? ((double)ctau[j] / err) : 0.0);

      // evtWeight.setVal(weight_final);
      NumDimu.setVal(nDimu);

      // cos_thetaVar.setVal(cos_theta[j]);
      // cos_theta1Var.setVal(cos_theta1[j]);
      if (hasPhi) phiVar.setVal(phi[j]);
      if (hasPhi1) phi1Var.setVal(phi1[j]);
      // recoQQdcaVar.setVal(recoQQdca[j]);
      
      if (hasCosHX) cosHXVar.setVal(cosHX[j]);
      if (hasPhiHX) phiHXVar.setVal(phiHX[j]);
      if (hasCosCS) cosCSVar.setVal(cosCS[j]);
      if (hasPhiCS) phiCSVar.setVal(phiCS[j]);
      // cos_epVar.setVal(cos_ep[j]);
      // phi_epVar.setVal(phi_ep[j]);

      ds.add(argSet);
      ++count_dimu;
    }
  }

  // ------------------------------------------------------------------
  // save meta info
  // ------------------------------------------------------------------
  TString metaInfo;
  metaInfo += Form("data_label=%s;", DATE.Data());
  metaInfo += Form("isMC=%d;", isMC);
  metaInfo += Form("MCtype=%d;", MCtype);
  metaInfo += Form("cLow=%d;", cLow);
  metaInfo += Form("cHigh=%d;", cHigh);
  metaInfo += Form("massLow=%.2f;", massLow);
  metaInfo += Form("massHigh=%.2f;", massHigh);
  metaInfo += Form("oppositeSign=%s;", oppositeSign ? "OS" : "SS");
  metaInfo += Form("fAccW=%d;", fAccW);
  metaInfo += Form("fEffW=%d;", fEffW);
  metaInfo += Form("isTnP=%d;", isTnP);
  metaInfo += Form("isPtW=%d;", isPtW);

  // run time
  TDatime now;
  metaInfo += Form("datetime=%04d-%02d-%02d_%02d:%02d:%02d;",
                   now.GetYear(), now.GetMonth(), now.GetDay(),
                   now.GetHour(), now.GetMinute(), now.GetSecond());

  // save info as ROOT object, TNamed
  TNamed config("skim_config", metaInfo.Data());
  config.Write();

  // record number of events and passed dimuons
  TParameter<int>("nEventsProcessed", nEvt).Write();
  TParameter<int>("nDimuonsSelected", count_dimu).Write();

  
  // ------------------------------------------------------------------
  // save RooDataSet
  // ------------------------------------------------------------------
  ds.Write();
  out_file->Close();

  // === error summary ===
  std::cout << "\n=========== Summary ===========\n";
  std::cout << "Processed events         : " << nEvt << "\n";
  std::cout << "Jpsi candidate           : " << count_dimu << "\n";
  std::cout << "non-finite ctau inputs   : " << nNonFiniteCtau << " times\n";
  std::cout << "ctau3DErr == 0           : " << nZeroCtauErr << " times\n";
  std::cout << "weight <= 0              : " << nZeroOrNegWeight << " times\n";
  std::cout << "null correction hist used: " << nNullHistoAccess << " times\n";
  std::cout << "===============================\n";

  std::cout << "=== Done ===\n";
  t.Stop();
  printf("RealTime=%f seconds, CpuTime=%f seconds\n", t.RealTime(), t.CpuTime());

  // close correction files
  // if (fin_correct)
  // {
  //   fin_correct->Close();
  //   delete fin_correct;
  //   fin_correct = nullptr;
  // }
}
