#ifndef syst_h
#define syst_h

#include "../../Fitter/Macros/Utilities/bin.h"
#include "../../Fitter/Macros/Utilities/resultUtils.h"
#include "TString.h"
#include "TSystemFile.h"
#include "TSystemDirectory.h"
#include "TCollection.h"
#include "TList.h"
#include <map>
#include <string>
#include <fstream>
#include <iostream>

#include "RooWorkspace.h"
#include "RooRealVar.h"

struct syst {
   string name;
   double value;
};

using namespace std;

vector<TString> fileList_syst(const char* token, const char* prependPath="");
map<anabin, syst> readSyst(const char* systfile, bool setSystZero=false);
map<anabin, syst> combineSyst(vector< map<anabin, syst> > theSysts, string name="Total");
map<anabin, syst> readSyst_all(const char* collSystem="", const char* poiname="NJpsi", const char* nameTag="16025", bool readEffSyst=false, bool readSigSyst=true,const char* prependPath="", bool doPrintTex=false, const char* texName="Systematics/systs.tex");
void printTex(vector< map<anabin, syst> > theSysts, const char* texName="Systematics/systs.tex", anabin mask=anabin(0,0,0,0,0,0));
map<anabin, vector<syst> > vm2mv(vector< map<anabin,syst> > v);
RooWorkspace* getWorkspaceFromBin(anabin thebin, const char* workDirName, const char* token="", const char* path2Fitter="");
void expandsyst(const char* input, const char* output, bool ptdep=false, bool rapdep=false, bool centdep=false);

#endif // ifndef syst_h
