#include <iostream>
#include <iomanip>
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

// from Run2 cutsAndBins.h
bool IsAcceptanceQQ(double pt, double eta)
{
  return ((fabs(eta) < 1.2 && pt >= 3.5) ||
          (1.2 <= fabs(eta) && fabs(eta) < 2.1 && pt >= 5.47 - 1.89 * fabs(eta)) ||
          (2.1 <= fabs(eta) && fabs(eta) < 2.4 && pt >= 1.5));
}

double getCorrectionSafe(TH2D *h, double pt, double angle)
{
  // return 1 when there is no histogram or contents are empty
  if (!h)
  {
    std::cout << "[Warning] getCorrectionSafe: histogram pointer is null → return 1\n";
    return 1.0;
  }

  // find bins including the target angle and pt
  int bx = h->GetXaxis()->FindBin(angle);
  int by = h->GetYaxis()->FindBin(pt);

  // when angle and pt are out of histogram range
  if (bx < 1 || bx > h->GetNbinsX() || by < 1 || by > h->GetNbinsY())
  {
    std::cout << "[Warning] getCorrectionSafe: value out of histogram range "
              << "(pt=" << pt << ", angle=" << angle << ") -> return 1\n";
    return 1.0;
  }

  double c = h->GetBinContent(bx, by);

  // Ususally, weight should be > 0
  if (c <= 0)
  {
    std::cout << "[Warning] getCorrectionSafe: bin content <= 0 at (bx=" << bx
              << ", by=" << by << ") -> return 1\n";
    return 1.0;
  }

  return 1.0 / c;
}

// use it later - apply the save guard like getCorrectionSafe()
double getAccWeight(TH1D *h = 0, double pt = 0)
{
  if (!h)
    return 1.0;
  int b = h->FindBin(pt);
  double c = h->GetBinContent(b);
  return (c > 0) ? 1.0 / c : 1.0;
}

// use it later
double getEffWeight(TH1D *h = 0, double pt = 0)
{
  if (!h)
    return 1.0;
  int b = h->FindBin(pt);
  double c = h->GetBinContent(b);
  return (c > 0) ? 1.0 / c : 1.0;
}

inline void EnableBranchStatus(TChain *ch, bool isMC)
{
  // 나중에 추가하기 -> bool enablePol = false

  ch->SetBranchStatus("*", 0);

  // event variables
  ch->SetBranchStatus("eventNb", 1);
  // ch->SetBranchStatus("cBin", 1);
  ch->SetBranchStatus("nDimu", 1);
  ch->SetBranchStatus("zVtx", 1);
  ch->SetBranchStatus("weight", 1);

  // dimuon
  ch->SetBranchStatus("recoQQsign", 1);
  ch->SetBranchStatus("mass", 1);
  ch->SetBranchStatus("pt", 1);
  ch->SetBranchStatus("y", 1);

  // single muon
  ch->SetBranchStatus("pt1", 1);
  ch->SetBranchStatus("pt2", 1);
  ch->SetBranchStatus("eta1", 1);
  ch->SetBranchStatus("eta2", 1);

  // dimuon lifetime
  ch->SetBranchStatus("ctau3D", 1);
  ch->SetBranchStatus("ctau3DErr", 1);
  if (isMC) ch->SetBranchStatus("ctau3DTrue", 1);

  // polarization
  // ch->SetBranchStatus("cos_theta", 1);
  // ch->SetBranchStatus("cos_theta1", 1);
  ch->SetBranchStatus("phi", 1);
  ch->SetBranchStatus("phi1", 1);
  // ch->SetBranchStatus("cos_cs", 1);
  // ch->SetBranchStatus("phi_cs", 1);
  // ch->SetBranchStatus("cos_hx", 1);
  // ch->SetBranchStatus("phi_hx", 1);
  // ch->SetBranchStatus("cos_ep", 1);
  // ch->SetBranchStatus("phi_ep", 1);
}


// ===== main macro =====
void skim_to_roodata(
    std::string dataLabel = "test",
    int nEvt = 10000,
    int cLow = 0, int cHigh = 180,
    float massLow = 2.6, float massHigh = 3.5,
    bool dimusign = true, // true=OS, false=SS
    bool isMC = false, int MCtype = 1,
    bool fAccW = false, bool fEffW = false,
    bool isTnP = false, bool isPtW = false,
    int hiHFBinEdge = 0,
    bool doQA = true)
{
  using namespace RooFit;
  TStopwatch t;
  t.Start();
  std::cout << "\n=== Convert FlowSkim to RooDataSet ===\n";

  gStyle->SetOptStat(0);
  gROOT->ForceStyle();

  TString DATE = dataLabel.c_str();

  // dimuon label - MUST use the opposite sign
  TString dimusignString = dimusign ? "OS" : "SS";

  // MC label - Jpsi MC must be PR or NP (CMS doesn't have the Inclusive Jpsi MC)
  TString MCtype_tag = (MCtype == 1) ? "PR" : (MCtype == 2 ? "NP" : "");

  // ===== FlowSkim input =====
  TChain *muon_chain = new TChain("myTree");
  muon_chain->Add("/data/users/pjgwak/work/daily_code_tracker/2025/250826/Pb23_flowSkim/FlowSkim_2023PbPbPromptRecoData_132X_miniAOD_noTrack_250826.root");

  // ===== Acc×Eff correction inputs =====
  TH2D *h_correct[6] = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};
  TFile *fin_correct = nullptr;

  // should distinguish - Acc
  // if (fAccW)
  // {}
  // Acc doesn't have a angle correction -> pT only

  // Eff
  if (fEffW)
  {
    fin_correct = TFile::Open("../eff_acc/roots/mc_eff_vs_pt_cent_rap_prompt_pbpb_Jpsi_PtW1_tnp1_250221.root", "READ");
    if (fin_correct && !fin_correct->IsZombie())
    {
      h_correct[0] = (TH2D *)fin_correct->Get("mc_eff_vs_pt_TnP1_PtW1_cent0to20_absy0_1p6");
      h_correct[1] = (TH2D *)fin_correct->Get("mc_eff_vs_pt_TnP1_PtW1_cent0to20_absy1p6_2p4");
      h_correct[2] = (TH2D *)fin_correct->Get("mc_eff_vs_pt_TnP1_PtW1_cent20_60_absy0_1p6");
      h_correct[3] = (TH2D *)fin_correct->Get("mc_eff_vs_pt_TnP1_PtW1_cent20_60_absy1p6_2p4");
      h_correct[4] = (TH2D *)fin_correct->Get("mc_eff_vs_pt_TnP1_PtW1_cent60to180_absy0_1p6");
      h_correct[5] = (TH2D *)fin_correct->Get("mc_eff_vs_pt_TnP1_PtW1_cent60to180_absy1p6_2p4");
    }
    else
    {
      std::cout << "[info] correction file not found. continue with weight=1.\n";
    }
  }

  // ===== bind input branch =====
  Int_t eventNb, cBin, nDimu;
  const int nMaxDimu = 1000;
  int recoQQsign[nMaxDimu];
  float zVtx, mass[nMaxDimu], pt[nMaxDimu], y[nMaxDimu];
  float pt1[nMaxDimu], pt2[nMaxDimu], eta[nMaxDimu], eta1[nMaxDimu], eta2[nMaxDimu];
  float ctau3D[nMaxDimu], ctau3DErr[nMaxDimu];
  float ctau3DTrue[nMaxDimu] = {0};
  double weight;

  float cos_theta[nMaxDimu], cos_theta1[nMaxDimu], phi[nMaxDimu], phi1[nMaxDimu];
  float cos_cs[nMaxDimu], phi_cs[nMaxDimu];
  float cos_hx[nMaxDimu], phi_hx[nMaxDimu];
  float cos_ep[nMaxDimu], phi_ep[nMaxDimu];
  
  // branch optimization
  EnableBranchStatus(muon_chain, isMC); // turn on the branches you need

  // bind branches
  muon_chain->SetBranchAddress("eventNb", &eventNb);
  // muon_chain->SetBranchAddress("cBin", &cBin);
  muon_chain->SetBranchAddress("nDimu", &nDimu);
  muon_chain->SetBranchAddress("recoQQsign", recoQQsign);
  muon_chain->SetBranchAddress("zVtx", &zVtx);
  muon_chain->SetBranchAddress("mass", mass);
  muon_chain->SetBranchAddress("y", y);
  muon_chain->SetBranchAddress("pt", pt);
  muon_chain->SetBranchAddress("pt1", pt1);
  muon_chain->SetBranchAddress("pt2", pt2);
  muon_chain->SetBranchAddress("eta", eta);
  muon_chain->SetBranchAddress("eta1", eta1);
  muon_chain->SetBranchAddress("eta2", eta2);
  muon_chain->SetBranchAddress("ctau3D", ctau3D);
  if (isMC) muon_chain->SetBranchAddress("ctau3DTrue", ctau3DTrue);
  muon_chain->SetBranchAddress("ctau3DErr", ctau3DErr);
  muon_chain->SetBranchAddress("weight", &weight);

  // muon_chain->SetBranchAddress("cos_theta", cos_theta);
  // muon_chain->SetBranchAddress("cos_theta1", cos_theta1);
  muon_chain->SetBranchAddress("phi", phi);
  muon_chain->SetBranchAddress("phi1", phi1);
  // muon_chain->SetBranchAddress("cos_cs", cos_cs);
  // muon_chain->SetBranchAddress("phi_cs", phi_cs);
  // muon_chain->SetBranchAddress("cos_hx", cos_hx);
  // muon_chain->SetBranchAddress("phi_hx", phi_hx);
  // muon_chain->SetBranchAddress("cos_ep", cos_ep);
  // muon_chain->SetBranchAddress("phi_ep", phi_ep);

  // ===== output file =====
  TFile *out_file = nullptr;
  if (isMC)
  {
    out_file = new TFile(Form("RooDataSet_miniAOD_isMC%d_%s_Jpsi_cent%i_%i_Effw%d_Accw%d_PtW%d_TnP%d_%s.root",
                              isMC, MCtype_tag.Data(), cLow, cHigh, fEffW, fAccW, isPtW, isTnP, DATE.Data()),
                         "RECREATE");
  }
  else
  {
    out_file = new TFile(Form("RooDataSet_miniAOD_isMC%d_Jpsi_cent%i_%i_Effw%d_Accw%d_PtW%d_TnP%d_%s.root",
                              isMC, cLow, cHigh, fEffW, fAccW, isPtW, isTnP, DATE.Data()),
                         "RECREATE");
  }

  // ===== Roo varaibles =====
  RooRealVar massVar("mass", "mass", 2.6, 3.5, "GeV/c^{2}");
  RooRealVar ptVar("pt", "pt", 0, 60, "GeV/c");
  RooRealVar yVar("y", "rapidity", -3, 3);
  RooRealVar pt1Var("pt1", "pt of muon+", 0, 500, "GeV/c");
  RooRealVar eta1Var("eta1", "eta of muon+", -4, 4);
  RooRealVar pt2Var("pt2", "pt of muon-", 0, 500, "GeV/c");
  RooRealVar eta2Var("eta2", "eta of muon-", -4, 4);
  // RooRealVar cBinVar("cBin", "Centrality", 0, 200);
  RooRealVar evtWeight("weight", "corr weight", 0, 10000);
  RooRealVar recoQQ("recoQQsign", "qq sign", -1, 3);
  RooRealVar ctau3DVar("ctau3D", "c_{#tau}", -10.0, 10.0, "mm");
  RooRealVar ctau3DErrVar("ctau3DErr", "#sigma_{c#tau}", 0, 10.0, "mm");
  RooRealVar ctau3DResVar("ctau3DRes", "c_{#tau} res", -100.0, 100.0);
  RooRealVar ctau3DTrueVar("ctau3DTrue", "c_{#tau} true", -10.0, 10.0, "mm");
  RooRealVar NumDimu("NumDimu", "number of dimuon", 0, 1000);

  // RooRealVar cos_thetaVar("cos_theta", "", -1.0, 1.0);
  // RooRealVar cos_theta1Var("cos_theta1", "", -1.0, 1.0);
  RooRealVar phiVar("phi", "", -5.0, 5.0);
  RooRealVar phi1Var("phi1", "", -5.0, 5.0);
  // RooRealVar cos_csVar("cos_cs", "", -1.0, 1.0);
  // RooRealVar phi_csVar("phi_cs", "", -5.0, 5.0);
  // RooRealVar cos_hxVar("cos_hx", "", -1.0, 1.0);
  // RooRealVar phi_hxVar("phi_hx", "", -5.0, 5.0);
  // RooRealVar cos_epVar("cos_ep", "", -1.0, 1.0);
  // RooRealVar phi_epVar("phi_ep", "", -5.0, 5.0);

  RooArgSet argSet(massVar, ptVar, yVar, pt1Var, pt2Var, eta1Var, eta2Var, evtWeight,
                   recoQQ, NumDimu, ctau3DVar, ctau3DErrVar, ctau3DResVar,
                   phiVar, phi1Var);
  // RooArgSet argSet(massVar, ptVar, yVar, pt1Var, pt2Var, eta1Var, eta2Var, evtWeight,
  //                  cBinVar, recoQQ, NumDimu, ctau3DVar, ctau3DErrVar, ctau3DResVar,
  //                  cos_thetaVar, cos_theta1Var, phiVar, phi1Var, cos_csVar, phi_csVar,
  //                  cos_hxVar, phi_hxVar, cos_epVar, phi_epVar);
  if (isMC) argSet.add(ctau3DTrueVar);

  // declare RooDataSet
  RooDataSet ds("dataset", "", argSet);

  // ===== others before main loop =====
  // error counter
  int nZeroCtauErr = 0;
  int nZeroOrNegWeight = 0;
  int nNullHistoAccess = 0;

  // set event number
  if (nEvt == -1)
    nEvt = muon_chain->GetEntries();
  std::cout << "Total events = " << nEvt << "\n";

  // print prgoress
  const Long64_t totalEvents = nEvt;
  const Long64_t reportEvery =
      (totalEvents >= 10'000'000) ? 200'000 :
      (totalEvents >= 1'000'000) ? 100'000 :
      (totalEvents >= 100'000) ? 20'000 : 10'000;

  // number of dimuon count
  int count_dimu = 0;

  // OS trigger -> Use OS
  const int signWanted = dimusign ? 0 : 1; // 0=OS, 1=SS

  // ===== main loop =====
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

    if (TMath::Abs(zVtx) > 15)
      continue;
    
    // 나중에 수정하기 -> 그냥 컷을 Onia to Skim에 다 넣어?
    // if (!(cBin >= cLow && cBin < cHigh))
    //   continue;

    // ===== dimuon loop =====
    for (int j = 0; j < nDimu; ++j)
    {
      // basic cut: pt<50, sign, mass, |y|<2.4, muon acc
      if (!(pt[j] < 50 &&
            recoQQsign[j] == signWanted &&
            mass[j] > massLow && mass[j] < massHigh &&
            TMath::Abs(y[j]) < 2.4 &&
            IsAcceptanceQQ(pt1[j], eta1[j]) &&
            IsAcceptanceQQ(pt2[j], eta2[j])))
        continue;

      // === acc×eff weight (default: 1) ===
      double weight_acc_eff = 1.0;
      // 나중에 할 거 - Acc, Eff 따로 켜고 끄기
      if (fAccW || fEffW)
      {
        const double absy = TMath::Abs(y[j]);
        // if (cBin < 20)
        // {
        //   weight_acc_eff = (absy < 1.6) ? getCorrectionSafe(h_correct[0], pt[j], cos_ep[j])
        //                                 : (absy < 2.4 ? getCorrectionSafe(h_correct[1], pt[j], cos_ep[j]) : 1.0);
        // }
        // else if (cBin >= 20 && cBin < 60)
        // {
        //   weight_acc_eff = (absy < 1.6) ? getCorrectionSafe(h_correct[2], pt[j], cos_ep[j])
        //                                 : (absy < 2.4 ? getCorrectionSafe(h_correct[3], pt[j], cos_ep[j]) : 1.0);
        // }
        // else if (cBin >= 60 && cBin < 180)
        // {
        //   weight_acc_eff = (absy < 1.6) ? getCorrectionSafe(h_correct[4], pt[j], cos_ep[j])
        //                                 : (absy < 2.4 ? getCorrectionSafe(h_correct[5], pt[j], cos_ep[j]) : 1.0);
        // }

        weight_acc_eff = 1;
      }
      const double weight_final = weight * weight_acc_eff; // 여기도 분리 weight * weight_acc * weight eff

      // correction error counter 
      if (!h_correct[0] && (fAccW || fEffW))
        nNullHistoAccess++; // can't read correction histogram
      if (weight_acc_eff <= 0)
        nZeroOrNegWeight++; // correction < 0

      // ----- push values into RooVariables -----
      recoQQ.setVal(recoQQsign[j]);
      massVar.setVal(mass[j]);
      ptVar.setVal(pt[j]);
      yVar.setVal(y[j]);
      pt1Var.setVal(pt1[j]);
      eta1Var.setVal(eta1[j]);
      pt2Var.setVal(pt2[j]);
      eta2Var.setVal(eta2[j]);
      // cBinVar.setVal(cBin);
      ctau3DVar.setVal(ctau3D[j]);
      ctau3DErrVar.setVal(ctau3DErr[j]);
      
      // prevent 0-division
      const double err = (double)ctau3DErr[j];
      if (err == 0)
        nZeroCtauErr++;
      ctau3DResVar.setVal((err != 0.0) ? ((double)ctau3D[j] / err) : 0.0);

      if (isMC) ctau3DTrueVar.setVal(ctau3DTrue[j]);

      evtWeight.setVal(weight_final);
      NumDimu.setVal(nDimu);

      // cos_thetaVar.setVal(cos_theta[j]);
      // cos_theta1Var.setVal(cos_theta1[j]);
      phiVar.setVal(phi[j]);
      phi1Var.setVal(phi1[j]);
      // cos_csVar.setVal(cos_cs[j]);
      // phi_csVar.setVal(phi_cs[j]);
      // cos_hxVar.setVal(cos_hx[j]);
      // phi_hxVar.setVal(phi_hx[j]);
      // cos_epVar.setVal(cos_ep[j]);
      // phi_epVar.setVal(phi_ep[j]);

      ds.add(argSet);
      ++count_dimu;
    }
  }

  // ===== save meta info =====
  TString metaInfo;
  metaInfo += Form("data_label=%s;", DATE.Data());
  metaInfo += Form("isMC=%d;", isMC);
  metaInfo += Form("MCtype=%d;", MCtype);
  metaInfo += Form("cLow=%d;", cLow);
  metaInfo += Form("cHigh=%d;", cHigh);
  metaInfo += Form("massLow=%.2f;", massLow);
  metaInfo += Form("massHigh=%.2f;", massHigh);
  metaInfo += Form("dimusign=%s;", dimusign ? "OS" : "SS");
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

  
  // ===== save RooDataSet =====
  ds.Write();

  // ===== simple qualitiy check =====
  if (doQA)
  {
    // mass
    {
      TCanvas c("c", "QA mass", 800, 600);
      RooPlot *f = massVar.frame(RooFit::Bins(50));
      ds.plotOn(f, RooFit::Cut("mass>2.6 && mass<3.5"));
      f->Draw();
      c.SaveAs(Form("QA_mass_%s.png", DATE.Data()));
    }

    // pt
    {
      TCanvas c("c", "QA pt", 800, 600);
      RooPlot *f = ptVar.frame(RooFit::Bins(50));
      ds.plotOn(f, RooFit::Cut("pt<50"));
      f->Draw();
      c.SaveAs(Form("QA_pt_%s.png", DATE.Data()));
    }

    // y
    {
      TCanvas c("c", "QA y", 800, 600);
      RooPlot *f = yVar.frame(RooFit::Bins(50));
      ds.plotOn(f, RooFit::Cut("abs(y)<2.4"));
      f->Draw();
      c.SaveAs(Form("QA_y_%s.png", DATE.Data()));
    }

    // centrality
    // {
    //   TCanvas c("c", "QA centrality", 800, 600);
    //   RooPlot *f = cBinVar.frame(RooFit::Bins(36)); // 0~180을 5씩
    //   ds.plotOn(f, RooFit::Cut("cBin>=0 && cBin<180"));
    //   f->Draw();
    //   c.SaveAs(Form("QA_centrality_%s.png", DATE.Data()));
    // }
  }

  // close the output file
  out_file->Close();

  // === error summary ===
  std::cout << "\n=========== Summary ===========\n";
  std::cout << "Processed events         : " << nEvt << "\n";
  std::cout << "Jpsi candidate           : " << count_dimu << "\n";
  std::cout << "ctau3DErr == 0           : " << nZeroCtauErr << " times\n";
  std::cout << "weight <= 0              : " << nZeroOrNegWeight << " times\n";
  std::cout << "null correction hist used: " << nNullHistoAccess << " times\n";
  std::cout << "===============================\n";

  std::cout << "=== Done ===\n";
  t.Stop();
  printf("RealTime=%f seconds, CpuTime=%f seconds\n", t.RealTime(), t.CpuTime());

  // close correction files
  if (fin_correct)
  {
    fin_correct->Close();
    delete fin_correct;
    fin_correct = nullptr;
  }
}
