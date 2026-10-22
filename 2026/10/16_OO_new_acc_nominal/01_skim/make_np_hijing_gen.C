#include "acceptance.h"
// ============================================================================
// NP HIJING: reco-selected J/psi의 generator transverse pseudoproper ctau 추출
// ============================================================================
// 입력: official OO2025 NP HIJING OniaTree의 hionia/myTree
// 출력: unweighted RooDataSet "dataset" + 생성 출처/선택 조건을 담은 ROOT metadata
// 범위: 1 < reco pT < 30 GeV/c, |reco y| < 2.4 (mid + forward 통합)
//
// 저장소에서 전체 이벤트 생성 (입력 경로: config/input_sources.json):
//   python3 tools/skim.py nonprompt
// 인자: (inputPath, outputPath="", maxEvents=-1)
//   outputPath="" : 이 매크로 옆에 아래 kOutputFileName으로 저장
//   maxEvents=-1   : 전체 이벤트 처리. 양수이면 별도 outputPath가 필요함.
// 기존 출력은 덮어쓰지 않는다. 재생성할 때는 다른 outputPath를 지정하거나
// 기존 파일을 먼저 옮겨야 한다.
//
// 리뷰 시 구분할 물리량:
//   mass/pt/y, ctau/ctauErr : 재구성된 dimuon의 값
//   truthCtau              : 매칭된 Gen J/psi의 Gen_QQ_ctau (이미 mm 단위)
// Gen_QQ_ctau3D나 B-hadron proper lifetime을 사용하지 않는다.
// ctau_np.C의 개별 bin 및 0 < truthCtau <= 6 mm fit cut은 여기서 적용하지 않는다.

#include <RooArgSet.h>
#include <RooDataSet.h>
#include <RooRealVar.h>

#include <TFile.h>
#include <TLorentzVector.h>
#include <TNamed.h>
#include <TParameter.h>
#include <TSystem.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

namespace np_hijing_gen
{
// --------------------------------------------------------------------------
// 공통 선택 설정
// 범위를 수정하면 출력 파일명과 마지막 recoSelection metadata도 함께 수정한다.
// --------------------------------------------------------------------------
constexpr double kMassLow = 2.6;
constexpr double kMassHigh = 3.5;
constexpr double kPtLow = 1.0;
constexpr double kPtHigh = 30.0;
constexpr double kAbsYHigh = 2.4;
constexpr ULong64_t kTriggerMask = 1ULL << 2; // HLT_OxyL1SingleMuOpen_v: bit 2
constexpr Long64_t kReportEvery = 1000000;
constexpr int kNumObservables = 8;

// 기존 RooDataSet 이름 규칙 + reco 선택/조상 매칭/Gen ctau 정의/운동학 범위.
constexpr const char *kOutputFileName =
    "RooDataSet_miniAOD_isMC1_NP_Hijing_globalOff_Jpsi_"
    "RecoSel_BAncestorMatched_GenCtau2D_pt1.00_30.00_y0.00_2.40_"
    "Effw0_Accw0_PtW0_TnP0_OO25_HLT_OxyL1SingleMuOpen_v1.root";

// 원본 forward 데이터와 동일한 b-hadron PDG ID 목록을 사용한다.
// 여기서 조상 계보를 새로 탐색하지 않고, OniaTree의 Gen_QQ_momId를 판정한다.
// abs()를 사용하므로 반입자도 포함된다.
bool isBAncestor(int pdg)
{
  switch (std::abs(pdg))
  {
    case 511:  // B0
    case 521:  // B+
    case 531:  // Bs0
    case 541:  // Bc+
    case 5122: // Lambda_b0
    case 5132: // Xi_b-
    case 5232: // Xi_b0
    case 5332: // Omega_b-
      return true;
    default:
      return false;
  }
}
} // namespace np_hijing_gen

void make_np_hijing_gen(
    const char *inputPath = "",
    const char *outputPath = "",
    Long64_t maxEvents = -1)
{
  using namespace np_hijing_gen;

  // ------------------------------------------------------------------------
  // 1. 실행 인자와 입력 파일 확인
  // ------------------------------------------------------------------------
  if (!inputPath || !*inputPath || !outputPath || maxEvents == 0 || maxEvents < -1)
    throw std::runtime_error("Provide an input path and maxEvents=-1 or a positive event count.");

  // 기본 출력은 현재 실행 디렉터리가 아니라 매크로가 위치한 디렉터리에 둔다.
  const std::string defaultOutput = std::string(gSystem->DirName(__FILE__)) + "/" + kOutputFileName;
  const std::string output = *outputPath ? outputPath : defaultOutput;
  if (maxEvents != -1 && !*outputPath)
    throw std::runtime_error("For a partial run, specify a separate outputPath.");
  if (!gSystem->AccessPathName(output.c_str()))
    throw std::runtime_error("Output already exists: " + output + "; move it aside or choose a new outputPath.");

  std::unique_ptr<TFile> input(TFile::Open(inputPath, "READ"));
  if (!input || input->IsZombie())
    throw std::runtime_error("Cannot open input OniaTree.");

  auto *tree = input->Get<TTree>("hionia/myTree");
  if (!tree)
    throw std::runtime_error("Missing hionia/myTree.");

  // ------------------------------------------------------------------------
  // 2. 필수 브랜치 확인 및 reader 연결
  // ------------------------------------------------------------------------
  // 필요한 브랜치가 없으면 이벤트 처리 전에 실패시킨다.
  // TTreeReaderArray는 입력의 vector/가변 길이 배열을 읽으며 고정 버퍼를 쓰지 않는다.
  const char *branches[] = {
      "HLTriggers",
      "Reco_QQ_4mom_m", "Reco_QQ_4mom_pt", "Reco_QQ_4mom_eta", "Reco_QQ_4mom_phi",
      "Reco_QQ_sign", "Reco_QQ_mupl_idx", "Reco_QQ_mumi_idx", "Reco_QQ_whichGen",
      "Reco_QQ_ctau", "Reco_QQ_ctauErr", "Reco_QQ_ctau3D", "Reco_QQ_ctauErr3D",
      "Reco_mu_4mom_pt", "Reco_mu_4mom_eta", "Reco_mu_whichGen",
      "Reco_mu_isSoftCutBased", "Reco_mu_isTracker", "Reco_mu_trig",
      "Gen_QQ_ctau", "Gen_QQ_momId"};
  for (const auto *name : branches)
  {
    if (!tree->GetBranch(name))
      throw std::runtime_error(std::string("Missing branch: ") + name);
  }

  TTreeReader reader(tree);

  // Event-level trigger.
  TTreeReaderValue<ULong64_t> trigger(reader, "HLTriggers");

  // Reconstructed dimuon: 운동학, 전하, muon/Gen index, 재구성 ctau.
  TTreeReaderArray<Float_t> qqMass(reader, "Reco_QQ_4mom_m");
  TTreeReaderArray<Float_t> qqPt(reader, "Reco_QQ_4mom_pt");
  TTreeReaderArray<Float_t> qqEta(reader, "Reco_QQ_4mom_eta");
  TTreeReaderArray<Float_t> qqPhi(reader, "Reco_QQ_4mom_phi");
  TTreeReaderArray<Short_t> sign(reader, "Reco_QQ_sign");
  TTreeReaderArray<Short_t> plus(reader, "Reco_QQ_mupl_idx");
  TTreeReaderArray<Short_t> minus(reader, "Reco_QQ_mumi_idx");
  TTreeReaderArray<Short_t> qqGen(reader, "Reco_QQ_whichGen");
  TTreeReaderArray<Float_t> qqCtau(reader, "Reco_QQ_ctau");
  TTreeReaderArray<Float_t> qqErr(reader, "Reco_QQ_ctauErr");
  TTreeReaderArray<Float_t> qqCtau3D(reader, "Reco_QQ_ctau3D");
  TTreeReaderArray<Float_t> qqErr3D(reader, "Reco_QQ_ctauErr3D");

  // Reconstructed muon: global-muon 조건은 사용하지 않는다 (globalOff).
  TTreeReaderArray<Float_t> muPt(reader, "Reco_mu_4mom_pt");
  TTreeReaderArray<Float_t> muEta(reader, "Reco_mu_4mom_eta");
  TTreeReaderArray<Short_t> muGen(reader, "Reco_mu_whichGen");
  TTreeReaderArray<Bool_t> soft(reader, "Reco_mu_isSoftCutBased");
  TTreeReaderArray<Bool_t> tracker(reader, "Reco_mu_isTracker");
  TTreeReaderArray<ULong64_t> muTrigger(reader, "Reco_mu_trig");

  // Generator J/psi: qqGen으로 매칭한 후보의 수명과 조상 PDG ID.
  TTreeReaderArray<Float_t> genCtau(reader, "Gen_QQ_ctau");
  TTreeReaderArray<Int_t> genMother(reader, "Gen_QQ_momId");

  // ------------------------------------------------------------------------
  // 3. 출력 RooDataSet 스키마
  // ------------------------------------------------------------------------
  // 원본 데이터와의 호환성을 위해 이름/순서/저장 범위를 유지한다.
  // 아래 RooRealVar 범위는 저장 가능한 범위이며, 물리적 선택 범위는 4절에 있다.
  // mass: GeV/c^2, pt: GeV/c, ctau/ctauErr/truthCtau: mm.
  RooRealVar mass("mass", "mass", 2.5, 3.6);
  RooRealVar pt("pt", "pt", 0., 300.);
  RooRealVar y("y", "y", -3., 3.);
  RooRealVar ctau("ctau", "ctau", -500., 1000.);
  RooRealVar ctauErr("ctauErr", "ctauErr", 0., 50.);
  RooRealVar genIndex("genIndex", "genIndex", -1., 1.e6);
  RooRealVar mother("mother", "mother", -1.e7, 1.e7);
  RooRealVar truthCtau("truthCtau", "truthCtau", -1000., 1000.);
  RooArgSet vars(mass, pt, y, ctau, ctauErr, genIndex, mother, truthCtau);
  RooRealVar *observables[kNumObservables] = {
      &mass, &pt, &y, &ctau, &ctauErr, &genIndex, &mother, &truthCtau};

  // WeightVar를 지정하지 않으며, Gen_weight/Eff/Acc/pT/TnP 가중치를 적용하지 않는다.
  RooDataSet dataset("dataset", "Dimuon matched to a b-hadron ancestor", vars);

  // ------------------------------------------------------------------------
  // 4. 이벤트 및 dimuon 후보 선택
  // ------------------------------------------------------------------------
  // 한 행은 선택된 dimuon 후보 하나이다. 동일 이벤트의 여러 후보도 각각 저장하며,
  // 입력 이벤트/후보 순서를 유지한다. centrality, zVtx, vertex-probability cut은 없다.
  const Long64_t events = maxEvents < 0 ? tree->GetEntries() : std::min(maxEvents, tree->GetEntries());
  Long64_t processed = 0;
  while (processed < events && reader.Next())
  {
    ++processed;
    if (processed % kReportEvery == 0)
      std::cout << "Events " << processed << "/" << events
                << ", candidates " << dataset.numEntries() << std::endl;

    // 4a. 해당 이벤트가 HLT bit 2를 통과해야 한다.
    if ((*trigger & kTriggerMask) == 0)
      continue;

    for (std::size_t i = 0; i < qqMass.GetSize(); ++i)
    {
      // 4b. 후보에 연결된 두 reco muon과 Gen J/psi index가 유효한지 확인한다.
      // 같은 collection의 브랜치 길이는 OniaTree 스키마에서 일치한다고 가정한다.
      const int muPlusIndex = plus[i];
      const int muMinusIndex = minus[i];
      const int genQQIndex = qqGen[i];
      if (muPlusIndex < 0 || muMinusIndex < 0 ||
          muPlusIndex >= int(muPt.GetSize()) || muMinusIndex >= int(muPt.GetSize()) ||
          genQQIndex < 0 || genQQIndex >= int(genCtau.GetSize()) || genQQIndex >= int(genMother.GetSize()))
        continue;

      // 4c. OS dimuon, mass 경계 포함 [2.6, 3.5], pT 경계 제외 (1, 30).
      if (qqMass[i] < kMassLow || qqMass[i] > kMassHigh || sign[i] != 0)
        continue;
      if (!(qqPt[i] > kPtLow && qqPt[i] < kPtHigh))
        continue;

      // 4d. 두 muon 모두 Gen-matched + soft + tracker + acceptance를 요구한다.
      // HLT matching은 두 muon 중 하나 이상이면 된다 (OR 조건).
      if (muGen[muPlusIndex] < 0 || muGen[muMinusIndex] < 0 ||
          !soft[muPlusIndex] || !soft[muMinusIndex] || !tracker[muPlusIndex] || !tracker[muMinusIndex])
        continue;
      if (((muTrigger[muPlusIndex] | muTrigger[muMinusIndex]) & kTriggerMask) == 0)
        continue;
      // The old acceptance is a strict superset; audit_acceptance proves equivalence.
      if (!oo_new_acc::passMuon(muPt[muPlusIndex], muEta[muPlusIndex]) ||
          !oo_new_acc::passMuon(muPt[muMinusIndex], muEta[muMinusIndex]))
        continue;

      // 4e. reco 4-vector로 rapidity를 계산한다. eta로 rapidity cut을 대체하지 않는다.
      // 기존 skim과 정확히 맞추기 위해 Rapidity() 결과를 float 정밀도로 저장한다.
      // |y| < 2.4에는 y=0과 mid/forward 경계도 포함된다.
      TLorentzVector qq;
      qq.SetPtEtaPhiM(qqPt[i], qqEta[i], qqPhi[i], qqMass[i]);
      const float rapidity = qq.Rapidity();
      if (!(std::abs(rapidity) < kAbsYHigh))
        continue;

      // 4f. 기존 skim의 3D ctau 유한성 조건을 유지한다.
      // 이것은 reco 품질 선택이며, 저장하는 ctau 자체는 transverse (2D) 값이다.
      if (!std::isfinite(qqCtau3D[i]) || !std::isfinite(qqErr3D[i]))
        continue;
      if (!isBAncestor(genMother[genQQIndex]))
        continue;

      // 4g. observables와 같은 순서로 값을 채운다. Gen_QQ_ctau는 이미 mm이므로
      // 단위 변환 없이 사용한다. 양의 truthCtau나 특정 fit window는 요구하지 않는다.
      const double values[kNumObservables] = {
          qqMass[i], qqPt[i], rapidity, qqCtau[i], qqErr[i],
          double(genQQIndex), double(genMother[genQQIndex]), genCtau[genQQIndex]};
      bool valid = qqErr[i] > 0;
      for (int j = 0; j < kNumObservables; ++j)
      {
        valid = valid && std::isfinite(values[j]) &&
                values[j] >= observables[j]->getMin() && values[j] <= observables[j]->getMax();
      }
      if (!valid)
        continue;

      for (int j = 0; j < kNumObservables; ++j)
        observables[j]->setVal(values[j]);
      dataset.add(vars);
    }
  }

  // ------------------------------------------------------------------------
  // 5. 처리 완료 확인 및 ROOT 저장
  // ------------------------------------------------------------------------
  // 읽기가 중간에 실패했거나 선택된 후보가 없으면 출력 파일을 만들지 않는다.
  if (processed != events)
    throw std::runtime_error("OniaTree read stopped before all requested events were processed.");
  if (dataset.numEntries() == 0)
    throw std::runtime_error("No selected candidates; no output written.");

  TFile result(output.c_str(), "CREATE");
  if (result.IsZombie())
    throw std::runtime_error("Cannot create output: " + output);
  if (dataset.Write() <= 0)
    throw std::runtime_error("Failed to write dataset.");

  // 별도 JSON 없이 ROOT 안에 출처, 선택, 수명 정의, 처리량을 기록한다.
  TNamed("sourceOniaTree", inputPath).Write();
  TNamed("singleMuonCut", oo_new_acc::description).Write();
  TNamed("acceptanceSourceSha256", oo_new_acc::sourceSha256).Write();
  TNamed("truthSelection",
         "genIndex>=0 && (abs(mother)==511 || abs(mother)==521 || abs(mother)==531 || "
         "abs(mother)==541 || abs(mother)==5122 || abs(mother)==5132 || "
         "abs(mother)==5232 || abs(mother)==5332)").Write();
  TNamed("truthDefinition",
         "truthCtau = Gen_QQ_ctau[Reco_QQ_whichGen], transverse pseudoproper ctau in mm; "
         "mass/pt/y are reconstructed").Write();
  TNamed("recoSelection",
         "OS; mass 2.6..3.5; 1<pt<30; abs(y)<2.4; event HLT and >=1 muon HLT bit 2; "
         "both muons gen-matched, soft, tracker, in acceptance; globalOff; "
         "finite ctau3D/ctauErr3D; ctauErr>0; finite observables in stored ranges; "
         "no centrality, vertex-probability, or zVtx cut; unweighted").Write();
  TParameter<Long64_t>("sourceEvents", tree->GetEntries()).Write();
  TParameter<Long64_t>("processedEvents", processed).Write();
  TParameter<int>("entries", dataset.numEntries()).Write();
  result.Close();

  std::cout << "Wrote " << dataset.numEntries() << " candidates to " << output << std::endl;
}
