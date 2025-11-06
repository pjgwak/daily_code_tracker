#ifndef createPPPbPbWorkspaces_C
#define createPPPbPbWorkspaces_C

#include "TString.h"
#include "TSystem.h"

#include <vector>
#include <map>
#include <cstdlib>
#include <iostream>
#include <sstream>

#include "../Fitter/Macros/Utilities/resultUtils.h"
#include "../Fitter/Systematics/syst.h"
#include "../Fitter/Systematics/syst.C"
#include "combinedWorkspace_4WS.C"

using namespace std;

// use the LLR batch system?
const bool usebatch=true;

void createPPPbPbWorkspaces_4WS(
                   const char* workDirName_pass, // workDirName: usual tag where to look for files in Output
                   const char* workDirName_fail, // workDirName: usual tag where to look for files in Output
                   bool doSyst=false,
                   const char* ACTag="Nominal", // Analysis Case tag (e.g. Nominal fits = "Nominal")
                   int nCPU=2
)
{
  // list of files
  vector<TString> theFiles_PbPb_pass = fileList(workDirName_pass,"PbPb","DATA","../Fitter");
  vector<TString> theFiles_PP_pass = fileList(workDirName_pass,"PP","DATA","../Fitter");
  vector<TString> theFiles_PbPb_fail = fileList(workDirName_fail,"PbPb","DATA","../Fitter");
  vector<TString> theFiles_PP_fail = fileList(workDirName_fail,"PP","DATA","../Fitter");
  if ( (theFiles_PbPb_pass.size() < 1) || (theFiles_PP_pass.size() < 1) )
  {
    cout << "#[Error]: No files found in " << workDirName_pass << endl;
    return;
  }
  if ( (theFiles_PbPb_fail.size() < 1) || (theFiles_PP_fail.size() < 1) )
  {
    cout << "#[Error]: No files found in " << workDirName_fail << endl;
    return;
  }
  
  // read systematics
  map<anabin, syst> syst_2R;
  if ( doSyst )
  {
    syst_2R = readSyst_all_prompt("","../Fitter",workDirName_pass,workDirName_fail);
    if ( syst_2R.empty() )
    {
      cout << "#[Error]: No systematics files found" << endl;
      return;
    }
  }
  
  int cnt=1;
  for (vector<TString>::const_iterator it_PbPb_pass=theFiles_PbPb_pass.begin(); it_PbPb_pass!=theFiles_PbPb_pass.end(); it_PbPb_pass++)
  {
    cout << endl;
    cout << "############ Merging workspaces for analysis bin " << cnt << " ..." << endl;
    cout << "PbPb pass workspace " << cnt << " / " << theFiles_PbPb_pass.size() << ": " << *it_PbPb_pass << endl;
    cout << endl;
    
    double systval(0.);
    anabin thebinPbPb = binFromFile(*it_PbPb_pass);
    if ( doSyst ) systval = syst_2R[thebinPbPb].value_dR;
    
    bool foundPPws = false;
    for (vector<TString>::const_iterator it_PP_pass=theFiles_PP_pass.begin(); it_PP_pass!=theFiles_PP_pass.end(); it_PP_pass++)
    {
      if ( !foundPPws && isSameBinPPPbPb(*it_PbPb_pass, *it_PP_pass) )
      {
        foundPPws = true;
        
        cout << "PP pass workspace " << cnt << " / " << theFiles_PP_pass.size() << ": " << *it_PP_pass << endl;
        cout << endl;
        
        bool foundPbPbfailws = false;
        for (vector<TString>::const_iterator it_PbPb_fail=theFiles_PbPb_fail.begin(); it_PbPb_fail!=theFiles_PbPb_fail.end(); it_PbPb_fail++)
        {
          anabin thebinPbPb_fail = binFromFile(*it_PbPb_fail);
          if (!foundPbPbfailws && (thebinPbPb_fail == thebinPbPb) )
          {
            foundPbPbfailws = true;
            
            cout << "PbPb fail workspace " << cnt << " / " << theFiles_PbPb_fail.size() << ": " << *it_PbPb_fail << endl;
            cout << endl;
            
            bool foundPPfailws = false;
            for (vector<TString>::const_iterator it_PP_fail=theFiles_PP_fail.begin(); it_PP_fail!=theFiles_PP_fail.end(); it_PP_fail++)
            {
              if ( !foundPPfailws && isSameBinPPPbPb(*it_PbPb_fail, *it_PP_fail) )
              {
                foundPPfailws = true;
                
                cout << "PP fail workspace " << cnt << " / " << theFiles_PP_fail.size() << ": " << *it_PP_fail << endl;
                cout << endl;
                
                TString binName(*it_PbPb_pass);
                binName.Remove(0,binName.Last('/')+1);
                binName.Remove(0,binName.Index("_pt")+1);
                binName.Remove(binName.First('.'),binName.Length());
                if ( doSyst ) binName.Prepend("wSyst_");
                else  binName.Prepend("woSyst_");
                
                if (!usebatch) {
                  combinedWorkspace_4WS(*it_PbPb_pass, *it_PbPb_fail, *it_PP_pass, *it_PP_fail, Form("combined_PbPbPP_PnP_workspace_%s.root",binName.Data()), systval, ACTag, nCPU);
                } else {
                  TString exports;
                  exports += Form("export it_PbPb_pass=%s; ",it_PbPb_pass->Data());
                  exports += Form("export it_PbPb_fail=%s; ",it_PbPb_fail->Data());
                  exports += Form("export it_PP_pass=%s; ",it_PP_pass->Data());
                  exports += Form("export it_PP_fail=%s; ",it_PP_fail->Data());
                  exports += Form("export binName=%s; ", binName.Data());
                  exports += Form("export systval=%f; ", systval);
                  exports += Form("export ACTag=%s; ", ACTag);
                  exports += Form("export nCPU=%i; ", nCPU);
                  exports += Form("export pwd_=%s; ", gSystem->pwd());
                  TString command("qsub -k oe -q cms@llrt3 ");
                  command += Form("-l nodes=1:ppn=%i ", nCPU);
                  command += "-l walltime=1:0:0 "; // 1h is enough and should make jobs run faster?
                  command += "-N ${binName} ";
                  command += "-V ";
                  command += Form("-o %s ", gSystem->pwd());
                  command += Form("-v it_PbPb_pass,it_PbPb_fail,it_PP_pass,it_PP_fail,binName,systVal,ACTag,nCPU,pwd_ ");
                  command += "runbatch_workspace_4WS.sh";
                  TString command_full = exports + command;
                  cout << command_full.Data() << endl;
                  system(command_full.Data());
                }
                
                cout << ">>>>>>>> Combined workspace created for bin " << binName.Data() << endl;
              }
              else continue;
            }
            if ( !foundPPfailws ) cout << "# [Error]: No PP fail workspace found for " << *it_PbPb_fail << endl;
          }
          else continue;
          
          if ( !foundPbPbfailws ) cout << "# [Error]: No PbPb fail workspace found for " << *it_PbPb_pass << endl;
        }
      }
      else continue;
    }
    if ( !foundPPws ) cout << "# [Error]: No PP pass workspace found for " << *it_PbPb_pass << endl;
    
    cnt++;
  } // loop on the files
}

#endif // #ifndef createPPPbPbWorkspaces_C
