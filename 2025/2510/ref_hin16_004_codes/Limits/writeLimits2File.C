#ifndef writeLimits2File_C
#define writeLimits2File_C

#include "TString.h"

#include <vector>
#include <sstream>
#include <fstream>

#include "../Fitter/Macros/Utilities/resultUtils.h"

#include "RooStats/HypoTestInverterResult.h"

using namespace std;

void writeLimits2File(
                   const char* ACTag // ACTag: where to look for limits files in Limits/Output/

)
{
  // list of files
  vector<TString> theFiles = limitsFileList(ACTag,"");
  if ( theFiles.empty() )
  {
    cout << "# [Error]: No limits files found" << endl;
    return;
  }
  
  // File to save results
  string limitsFileName = "";
  ofstream file;
  
  // Loop over files
  int cnt=1;
  for (vector<TString>::const_iterator it=theFiles.begin(); it!=theFiles.end(); it++)
  {
    // Get the calculator type, CL and sys from file name
    TString sfileName(*it);
    sfileName.Remove(0,sfileName.Last('/')+1);
    TString sSys("woSys");
    if ( sfileName.Index("wSys") != kNPOS ) sSys = "wSys";
    
    TString scalculatorType(sfileName);
    scalculatorType.Remove(scalculatorType.First('_'),scalculatorType.Sizeof());
    
    TString sCL(sfileName);
    sCL.Remove(0,sCL.Last('_')+1);
    sCL.Remove(sCL.First('.'),sCL.Sizeof());
    
    if (cnt == 1) // Create the file to write the limits
    {
      limitsFileName = string("csv/") + "cLimits_" + string(sCL) + "_" + string(ACTag) + "_" + string(scalculatorType) + "_" + string(sSys) + ".csv";
      file.open(limitsFileName.c_str());
    }
    
    sCL.Prepend("0.");
    if (cnt == 1) file << sCL.Data() << endl;
    
    cout << ">>>>>>> Writting " << Form("%.2f",sCL.Atof()*100) << "% " << "limits on " << limitsFileName.c_str() << " for analysis bin " << cnt << endl;
    
    // Get the bin info from file name
    TString sbinInfo(sfileName);
    sbinInfo.Remove(0,sbinInfo.Index("_pt")+1);
    sbinInfo.Remove(sbinInfo.Last('_'),sbinInfo.Sizeof());
    
          // pt
    TString sPt(sbinInfo);
    sPt.Remove(sPt.First('_'),sPt.Sizeof());
    sPt.ReplaceAll("pt","");
    int sSize = sPt.Sizeof();
    int nPos = 0;
    if ( sSize == 7 ) nPos = 3;
    else if ( sSize == 6 || sSize == 5 ) nPos = 2;
    
    TString sptMin(sPt);
    sptMin.Remove(nPos,sptMin.Sizeof());
    TString sptMax(sPt);
    sptMax.Remove(0,nPos);
    
    sptMin.Insert(sptMin.Sizeof()-2,'.');
    sptMax.Insert(sptMax.Sizeof()-2,'.');
    
    
          // rap
    TString sRap(sbinInfo);
    sRap.Remove(0,sRap.First('_')+1);
    sRap.Remove(sRap.First('_'),sRap.Sizeof());
    sRap.ReplaceAll("rap","");
    
    sSize = sRap.Sizeof();
    if (sRap.Index("0") == 0) nPos = 1;
    else if (sSize == 5) nPos = 2;
    TString sRapMin(sRap);
    sRapMin.Remove(nPos,sRapMin.Sizeof());
    TString sRapMax(sRap);
    sRapMax.Remove(0,nPos);
    
    sRapMin.Insert(sRapMin.Sizeof()-2, nPos == 1 ? "0." : ".");
    sRapMax.Insert(sRapMax.Sizeof()-2,'.');
    
    
          // cent
    TString sCent(sbinInfo);
    sCent.Remove(0,sCent.Last('_')+1);
    sCent.ReplaceAll("cent","");
    
    TString sCentMin(sCent);
    TString sCentMax(sCent);
    sSize = sCent.Sizeof();
    if (sCent.Index("0") == 0) nPos = 1;
    else if (sSize == 6 || sSize == 5) nPos = 2;
    else if (sSize == 7) nPos = 3;
    sCentMin.Remove(nPos,sCentMin.Sizeof());
    sCentMax.Remove(0,nPos);
    
    
    // Read the HypotestInverterResult to get the limits
    TFile* f = TFile::Open(*it,"READ");
    RooStats::HypoTestInverterResult* r = static_cast<RooStats::HypoTestInverterResult*>(f->FindObjectAny("result_RFrac2Svs1S_PbPbvsPP"));
    if (!r) r = static_cast<RooStats::HypoTestInverterResult*>(f->FindObjectAny("result_RFrac2Svs1S_PbPbvsPP_P"));
    if (!r) {
       cout << "[ERROR] Could not find HypoTestInverterResult in file "<< *it << endl;
       continue;
    }
    
    // Write results to file
    file << sRapMin.Data() << ", " << sRapMax.Data() << ", "
       << sptMin.Data() << ", " << sptMax.Data() << ", "
       << sCentMin.Data() << ", " << sCentMax.Data() << ", "
       << r->LowerLimit() << ", " << r->UpperLimit() << endl;
    
    f->Close(); delete f;
    
    cnt++;
  } // loop on the files
  file.close();
  cout << "Closed " << limitsFileName << endl << endl;
}

#endif // #ifndef writeLimits2File_C
