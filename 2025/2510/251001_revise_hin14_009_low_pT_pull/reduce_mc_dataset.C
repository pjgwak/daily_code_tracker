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


void setSilent() {
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

void setVarRange(RooAbsData *data, float lmin, float lmax, float errmin, float errmax)
{
  RooRealVar *mass = (RooRealVar *)data->get()->find("mass");
  RooRealVar *ctau3D = (RooRealVar *)data->get()->find("ctau3D");
  RooRealVar *ctau3DErr = (RooRealVar *)data->get()->find("ctau3DErr");
  RooRealVar *ctau3Dtrue = (RooRealVar *)data->get()->find("ctau3Dtrue");

  if (mass)
    mass->setRange(2.6, 3.5);
  if (ctau3D)
    ctau3D->setRange(lmin,lmax); // PR MC
  if (ctau3Dtrue)
    ctau3Dtrue->setRange(-1, 9);

  return;
}

void readInputs(RooDataSet *&dataPRMC, RooDataSet *&dataNPMC, const bool isWeight)
{
  // calculate average weight (first 100 evets)
  auto printAvgWeight = [](RooDataSet *ds, int nCheck = 100)
  {
    // if (!ds || !ds->weightVar())
    // {
    //   cout << "   [Check] No weight variable found" << endl;
    //   return;
    // }
    int nEntries = ds->numEntries();
    int nLoop = std::min(nCheck, nEntries);

    double sumW = 0.0;
    for (int i = 0; i < nLoop; i++)
    {
      ds->get(i); // read event
      sumW += ds->weight();
    }
    double avgW = (nLoop > 0) ? sumW / nLoop : 0.0;
    cout << "   [Check] Avg weight of first " << nLoop << " events = " << avgW << endl;

    if (fabs(avgW - 1.0) < 1e-6) {
      cout << "\033[31m" // color red
           << "   [WARNING] It seems weighting is 1 - UN-WEIGHTED dataset. Is it right?\n"
           << "\033[0m";
    }
  };

  // --- Prompt MC ---
  string fileNamePrMc = "/data/hwan/psi2S_RAA_PbPb2018/skimmedFiles/OniaRooDataSet_isMC1_Jpsi_pp_y0.00_2.40_Effw0_Accw0_PtW0_TnP0_230717.root";
  TFile fInPRMC(fileNamePrMc.c_str());
  cout << "Load: " << fileNamePrMc.c_str() << endl;
  dataPRMC = (RooDataSet *)fInPRMC.Get("dataset");
  if (isWeight) {
    cout << "[Info] dataPRMC is WEIGHTED\n";

    dataPRMC = new RooDataSet("dataPRMC", "dataset with weight",
                              dataPRMC, *dataPRMC->get(), 0,
                              "weight");
    printAvgWeight(dataPRMC, 100);
  } else
  {
    cout << "[Info] dataPRMC is UN-WEIGHTED" << endl;
    dataPRMC->SetName("dataPRMC");
  }

  // --- NP MC ---
  string fileNameNpMc = "/data/hwan/psi2S_RAA_PbPb2018/skimmedFiles/OniaRooDataSet_JPsi_GENONLY_NonPrompt_y0_2p4_230829.root";
  TFile fInNPMC(fileNameNpMc.c_str());
  cout << "Load: " << fileNameNpMc.c_str() << endl;
  dataNPMC = (RooDataSet *)fInNPMC.Get("dataset");

  // turn off the Weight for MC Gen
  if (isWeight && false) {
    cout << "[Info] dataNPMC is WEIGHTED\n";

    dataNPMC = new RooDataSet("dataNPMC", "dataset with weight",
                              dataNPMC, *dataNPMC->get(), 0,
                              "weight");
    printAvgWeight(dataNPMC, 100);
  } else
  {
    cout << "[Info] dataNPMC is UN-WEIGHTED" << endl;
    cout << "\033[31m" // color red
         << "Hard-coded to turn off the weight for GenOnly MC" << endl
         << "\033[0m";
    dataNPMC->SetName("dataNPMC");
  }
}

void reduce_mc_dataset()
{
  cout << "\n=== Start reduce_mc_dataset() ===\n";
  TStopwatch t;
  t.Start();

  // silent
  setSilent();

  // rootlogon -> To global config
  gROOT->Macro("/data/users/pjgwak/input_files/rootlogon.C");

  // === read inputs ===
  cout << "=== Import inputs ===\n";
  RooDataSet *dataPRMC = nullptr, *dataNPMC = nullptr;
  
  readInputs(dataPRMC, dataNPMC, isWeight);

  // === reduce dataset ===
  cout << "\n=== reduce dataset ===\n";

  // --- make cuts ---
  // muon kinematic cut
  TCut muonCut =
      "((abs(eta1) <= 1.2 && pt1 >= 3.5) || "
      "(abs(eta2) <= 1.2 && pt2 >= 3.5) || "
      "((abs(eta1) > 1.2 && abs(eta1) <= 2.1) && pt1 >= 5.47-1.89*abs(eta1)) || "
      "((abs(eta2) > 1.2 && abs(eta2) <= 2.1) && pt2 >= 5.47-1.89*abs(eta2)) || "
      "((abs(eta1) > 2.1 && abs(eta1) <= 2.4) && pt1 >= 1.5) || "
      "((abs(eta2) > 2.1 && abs(eta2) <= 2.4) && pt2 >= 1.5))";

  // Opposite sign
  TCut pairCut = "recoQQsign == 0";

  // kinematic bins cuts
  TCut kineCut = Form("(pt >= %.3f && pt < %.3f && abs(y) >= %.3f && abs(y) < %.3f && mass >= %.3f && mass < %.3f)",
                      ptLow, ptHigh, yLow, yHigh, massLow, massHigh);

  // ctau3D 
  TCut ctCut = Form("(ctau3D >= %.3f && ctau3D < %.3f)", ctMin, ctMax);

  // MC true ctau
  TCut ctTrueCut = "(ctau3Dtrue >= 0.0001 && ctau3Dtrue < 9)";

  // total cuts
  TCut reduceDS_woCtErr = muonCut && pairCut && kineCut && ctCut; // PR MC
  TCut reduceNpMc = muonCut && kineCut && ctTrueCut;                // NP MC

  // --- reduce ---
  // 필요한 변수만 리턴해서 저장? -> 가장 큰 파일도 650 MB 정도라서 용량에 여유가 있다.
  // -> 코드 복잡도를 줄이기 위해 그냥 통째로 저장
  RooDataSet *redPRMC, *redNPMC;
  redPRMC = (RooDataSet *)dataPRMC->reduce(reduceDS_woCtErr); // no ctauErr cut
  redPRMC->SetName("redPRMC");
  redNPMC = (RooDataSet *)dataNPMC->reduce(reduceNpMc);
  redNPMC->SetName("redNPMC");

  setVarRange(redPRMC, -0.2, 0.2, errmin, errmax);
  setVarRange(redNPMC, ctMin, ctMax, errmin, errmax);

  // === set variable ranges in ws ===  
  // we will use this workspace from now on
  // setWSRange(ws, -ctMin, ctMax, errmin, errmax); -> 이거 global config로 만들고 err pdf 이후로 자동 호출하도록 설정
  // ws->var("mass")->SetTitle("m_{#mu#mu}");
  // ws->var("ctau3D")->SetTitle("#font[12]{l}_{J/#psi}");

  // === save results ===
  cout << "=== save resuts ===\n";
  gSystem->mkdir("reduced_data", true);
  TFile fileOut(Form("reduced_data/mc_pT%.1f_%.1f_y%.1f_%.1f_isWeight%d.root", ptLow, ptHigh, yLow, yHigh, isWeight), "recreate");
  redPRMC->Write();
  redNPMC->Write();
  
  // --- test - to avoid the error ---
  // Error in <TBufferFile::WriteByteCount>: bytecount too large (more than 1073741822)
  // TTree *tree = (TTree *)redPRMC->GetClonedTree();
  // tree->SetName("treePRMC");
  // tree->Write();

  fileOut.Close();

  cout << "=== Finish reduce_mc_dataset() ===\n";
  t.Stop();
  printf("RealTime=%f seconds, CpuTime=%f seconds\n", t.RealTime(), t.CpuTime());
}