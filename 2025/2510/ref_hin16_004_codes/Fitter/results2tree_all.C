#ifndef results2tree_all_C
#define results2tree_all_C

#include <iostream>

#include "TString.h"

#include "results2tree.C"

using namespace std;

void results2tree_all(bool is16004 = false)
{
  TString app("");
  if (is16004) app = "16004";
  else app = "16025";
  
  TString workdir("");
  // Nominal
  workdir = Form("DataFits_%s_2D_2CB_polBkg_nominal",app.Data());
//  results2tree(workdir.Data());
  
  // Background variations
  //workdir = Form("DataFits_%s_2D_2CB_polBkg_minvRange",app.Data());
  //results2tree(workdir.Data());
  workdir = Form("DataFits_%s_2D_2CB_polBkg_LLR25",app.Data());
 // results2tree(workdir.Data());
  workdir = Form("DataFits_%s_2D_2CB_polBkg_LLR10",app.Data());
  //results2tree(workdir.Data());
  workdir = Form("DataFits_%s_2D_2CB_expBkg",app.Data());
  //results2tree(workdir.Data());
  
  // Signal variations
  workdir = Form("DataFits_%s_2D_2CB_polBkg_constrained",app.Data());
  results2tree(workdir.Data());
return;  
workdir = Form("DataFits_%s_2D_CBG_polBkg",app.Data());
  results2tree(workdir.Data());
  
  // CtauErr variations
  workdir = Form("DataFits_%s_2D_2CB_polBkg_ctauErrTot",app.Data());
  results2tree(workdir.Data());
  
  // CtauTrue variations
  workdir = Form("DataFits_%s_2D_2CB_polBkg_ctauRecoTemp",app.Data());
  results2tree(workdir.Data());
  
  // CtauRes variations
  workdir = Form("DataFits_%s_2D_2CB_polBkg_ctauResPromptMC",app.Data());
  results2tree(workdir.Data());
  
  // CtauBkg variations
  workdir = Form("DataFits_%s_2D_2CB_polBkg_ctauBkgTemp",app.Data());
  results2tree(workdir.Data());
  
  return;
}

#endif
