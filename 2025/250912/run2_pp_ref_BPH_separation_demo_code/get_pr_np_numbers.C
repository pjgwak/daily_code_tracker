#include <TStopwatch.h>

using namespace RooFit;

void get_pr_np_numbers()
{
  float ptLow = 25, ptHigh = 40;
  float yLow = 0, yHigh = 1.6;
  float massLow = 2.6, massHigh = 3.5;
  std::string comp = "", region = "";
  int cLow = 0, cHigh = 180;

  // load inputs
  // mass PR Full - # of mass PRS, f_PR_Bkg (tmp)
  TFile *fMassPrFull = TFile::Open(Form("roots/mass_PR_MassFull_pT%.1f_%.1f_y%.1f_%.1f.root", ptLow, ptHigh, yLow, yHigh));
  auto *nSig_pr = (TParameter<double> *)fMassPrFull->Get("n_sig_sub");
  auto *f_pr_bkg = (TParameter<double> *)fMassPrFull->Get("f_mass_bkg");
  

  // mass NP Full - # of mass NPS, f_NP_Bkg (tmp)
  TFile *fMassNpFull = TFile::Open(Form("roots/mass_NP_MassFull_pT%.1f_%.1f_y%.1f_%.1f.root", ptLow, ptHigh, yLow, yHigh));
  auto *nSig_np = (TParameter<double> *)fMassNpFull->Get("n_sig_sub");
  auto *f_np_bkg = (TParameter<double> *)fMassNpFull->Get("f_mass_bkg");

  // mass PRSB - # of mass PRSB (LSB onlt for test)
  TFile *fMassPrLsb = TFile::Open(Form("roots/mass_PR_LSB_pT%.1f_%.1f_y%.1f_%.1f.root", ptLow, ptHigh, yLow, yHigh));
  RooFitResult *resultMassPrLsb = (RooFitResult *)fMassPrLsb->Get("fitResult");
  const RooArgList &paramsMassPrLsb = resultMassPrLsb->floatParsFinal();
  RooRealVar *NsigMassPrLsb = (RooRealVar *)paramsMassPrLsb.find("N");

  // mass NPSB - # of mass NPSB (?) - (LSB onlt for test)
  TFile *fMassNpLsb = TFile::Open(Form("roots/mass_NP_LSB_pT%.1f_%.1f_y%.1f_%.1f.root", ptLow, ptHigh, yLow, yHigh));
  RooFitResult *resultMassNpLsb = (RooFitResult *)fMassNpLsb->Get("fitResult");
  const RooArgList &paramsMassNpLsb = resultMassNpLsb->floatParsFinal();
  RooRealVar *NsigMassNpLsb = (RooRealVar *)paramsMassNpLsb.find("N");

  // mass PR Center - f_PR_Bkg -> skip
  // mass NP Center - f_NP_Bkg-> skip

  // ctauFit - f_PR_JpsiB
  TFile *fCtauPr = TFile::Open(Form("roots/ctau_CtauFull_SR_pT%.1f_%.1f_y%.1f_%.1f.root", ptLow, ptHigh, yLow, yHigh));
  auto *b_frac = (TParameter<double> *)fCtauPr->Get("b_frac");

  cout << "\n--- Componetns ---\n";
  cout << "nMassPrs:     " << nSig_pr->GetVal() << "\n";
  cout << "f_Pr_bkg:     " << f_pr_bkg->GetVal() << "\n";
  cout << "nMassNps:     " << nSig_np->GetVal() << "\n";
  cout << "f_Np_bkg:     " << f_np_bkg->GetVal() << "\n";
  cout << "nMass PR LSB: " << NsigMassPrLsb->getVal() << "\n";
  cout << "nMass NP LSB: " << NsigMassNpLsb->getVal() << "\n";
  cout << "f_PR_JpsiB:   " << b_frac->GetVal() << "\n";

  cout << "\n--- # of Jpsi ---\n";
  cout << "PR: " << nSig_pr->GetVal() - (b_frac->GetVal() * (nSig_np->GetVal() - (f_np_bkg->GetVal() * NsigMassNpLsb->getVal()))) - (f_pr_bkg->GetVal() * NsigMassPrLsb->getVal()) << "\n";
  cout << "NP: " << nSig_np->GetVal() - (f_np_bkg->GetVal() * NsigMassNpLsb->getVal()) << "\n";
}