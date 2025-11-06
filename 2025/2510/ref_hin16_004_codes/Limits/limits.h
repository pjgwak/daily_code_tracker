#ifndef limits_h
#define limits_h

#include <algorithm>
#include "../Fitter/Macros/Utilities/bin.h"
#include "../Fitter/Systematics/syst.h"
#include "../Fitter/Macros/Utilities/resultUtils.h"
#include "../Fitter/Macros/Utilities/monster.h"
#include "TString.h"
#include "TSystemFile.h"
#include "TSystemDirectory.h"
#include <map>
#include <string>
#include <fstream>
#include <iostream>

struct limits {
   double cl;
   pair<double,double> val;
};

using namespace std;

vector<TString> fileList_lim(const char* prependPath="");
map<anabin, limits> readLimits(const char* limitsfile);
void printTex(vector< map<anabin, limits> > theLims, const char* texName="tex/limits.tex");
void readLimits_all(const char* prependPath="", const char* texName="tex/limits.tex", const char* token="_68");
void limitsFromFits(const char* workDir, const char* output="csv/Limits_68_fromfits.csv", bool withsyst=false, int mode=0, const char* workDirFail="");
// mode = 0 (pass), 1 (prompt) or 2 (non-prompt)
map<anabin, vector<limits> > vm2mv(vector< map<anabin,limits> > v);


vector<TString> fileList_lim(const char* prependPath) {
   vector<TString> ans;

   TString basedir("csv/");
   if ( strcmp(prependPath,"") ) basedir.Prepend(Form("%s/",prependPath));
   TSystemDirectory dir("files",basedir);

   TList *files = dir.GetListOfFiles();

   if (files) {
      TIter next(files);
      TSystemFile *file;
      TString fname;

      while ((file=(TSystemFile*)next())) {
         fname = file->GetName();
         if (fname.EndsWith(".csv")) {
            ans.push_back(basedir+fname);
         }
      }
   }

   return ans;
}

map<anabin, limits> readLimits(const char* limitsfile) {
   map<anabin, limits> ans;

   ifstream file(limitsfile);
   if (!(file.good())) return ans;

   string limitscl; getline(file,limitscl);
   double cl = atof(limitscl.c_str());

   string line;
   double rapmin=0, rapmax=0, ptmin=0, ptmax=0, centmin=0, centmax=0, value=0;
   limits lim;
   lim.cl = cl;

   while (file.good()) {
      getline(file,line);
      if (line.size()==0) break;
      TString tline(line.c_str());
      TString t; Int_t from = 0, cnt=0;
      while (tline.Tokenize(t, from , ",")) {
         t.Strip(TString::kBoth,' ');
         value = atof(t.Data());
         if (cnt==0) rapmin = atof(t.Data());
         else if (cnt==1) rapmax = value;
         else if (cnt==2) ptmin = value;
         else if (cnt==3) ptmax = value;
         else if (cnt==4) centmin = value;
         else if (cnt==5) centmax = value;
         else if (cnt==6) lim.val.first = value;
         else if (cnt==7) lim.val.second = value;
         else if (cnt>7) {
            cout << "Warning, too many fields, I'll ignore the last ones." << endl;
            continue;
         }
         cnt++;
      }
      anabin thebin(rapmin, rapmax, ptmin, ptmax, centmin, centmax);
      ans[thebin] = lim;
   }

   return ans;
}

void readLimits_all(const char* prependPath, const char* texName, const char* token) {
   vector< map<anabin, limits> > limmap_all;
   vector<TString> filelist = fileList_lim(prependPath);

   for (vector<TString>::const_iterator it=filelist.begin(); it!=filelist.end(); it++) {
      if (it->Index(token) == kNPOS) continue; 
      cout << "Reading file " << *it << endl;
      map<anabin,limits> limmap = readLimits(it->Data());
      limmap_all.push_back(limmap);
   }

   printTex(limmap_all, texName);
}

void printTex(vector< map<anabin, limits> > theLims, const char* texName) {
   unsigned int nlimits = theLims.size();

   ofstream file(texName);
   file << "\\begin{tabular}{|ccc|"; 
   for (unsigned int i=0; i<nlimits; i++) {
      file << "c|";
   }
   file << "}" << endl;
   file << "\\hline" << endl;
   file << "$|y|$ & \\pt & Centrality";
   for (unsigned int i=0; i<nlimits; i++) file << " & " << Form("%.0f%s C.L.",100.*theLims[i].begin()->second.cl,"\\%");
   file<< "\\\\" << endl;
   file << "\\hline" << endl;

   anabin oldbin(0,0,0,0,0,0);
   map<anabin, vector<limits> > themap = vm2mv(theLims);
   map<anabin, vector<limits> >::const_iterator itm;
   for (itm=themap.begin(); itm!=themap.end(); itm++) {
      vector<limits> v = itm->second;
      if (v.size() != nlimits) {
         cout << "Error, not all limitsematics have the same size. Aborting." << endl;
         file.close();
         return;
      }
      anabin thebin = itm->first;
      if (thebin.rapbin() == oldbin.rapbin()) {
         file << " - ";
      } else {
         if (itm != themap.begin()) file << "\\hline" << endl;
         file.precision(2); file.setf(ios::fixed);
         file << thebin.rapbin().low() << " $< |y| < $ " << thebin.rapbin().high();
      }
      file << " & ";
      if (thebin.ptbin() == oldbin.ptbin()) {
         file << " - ";
      } else {
         file.precision(2); file.setf(ios::fixed);
         file << thebin.ptbin().low() << " $< \\pt < $ " << thebin.ptbin().high();
      }
      file << " & ";
      file.unsetf(ios::fixed);
      file << thebin.centbin().low()/2 << "-" << thebin.centbin().high()/2 << "\\% ";
      file.precision(2);
      file.setf(ios::fixed);

      for (unsigned int i=0; i<nlimits; i++) {
         file << " & $[" << v[i].val.first << ", " << v[i].val.second << "]$";
      }
      file << " \\\\" << endl;

      oldbin = thebin;
   }

   file << "\\hline" << endl;
   file << "\\end{tabular}" << endl;
   file.close();
   cout << "Closed " << texName << endl;
}

void limitsFromFits(const char* workDir, const char* output, bool withsyst, int mode, const char* workDirFail) {
   set<anabin> thebins = allbins();
   const char* ppp = "../Fitter";

   map<anabin, syst> systs_all;
   if (withsyst) {
      if (mode==0) systs_all = readSyst_all_pass("",ppp,workDir);
      if (mode==1) systs_all = readSyst_all_prompt("",ppp,workDir,workDirFail);
      if (mode==2) systs_all = readSyst_all_nonprompt("",ppp,workDir,workDirFail);
   }

   ofstream file(output);
   file << "0.683" << endl;

   set<anabin>::const_iterator it;
   for (it=thebins.begin(); it!=thebins.end(); it++) {
      anabin thebin = *it;
      double val, err;
      if (mode==0) {
         val = doubleratio_pass_nominal(workDir,thebin,ppp);
         err = doubleratio_pass_stat(workDir,thebin,ppp);
      }
      if (mode==1) {
         val = doubleratio_prompt_nominal(workDir,workDirFail,thebin,ppp);
         err = doubleratio_prompt_stat(workDir,workDirFail,thebin,ppp);
      }
      if (mode==2) {
         val = doubleratio_nonprompt_nominal(workDir,workDirFail,thebin,ppp);
         err = doubleratio_nonprompt_stat(workDir,workDirFail,thebin,ppp);
      }

      file << thebin.rapbin().low() << ", " << thebin.rapbin().high() << ", "
         << thebin.ptbin().low() << ", " << thebin.ptbin().high() << ", "
         << thebin.centbin().low() << ", " << thebin.centbin().high() << ", ";
      if (!withsyst) {
         file << val - err << ", " << val + err << endl;
      } else {
         file << val - sqrt(pow(err,2) + pow(systs_all[thebin].value_dR,2)) << ", " 
            << val + sqrt(pow(err,2) + pow(systs_all[thebin].value_dR,2)) << endl;
      }
   }

   file.close();
}

map<anabin, vector<limits> > vm2mv(vector< map<anabin,limits> > v) {
   map<anabin, vector<limits> > ans;

   vector<map<anabin,limits> >::const_iterator itv;
   map<anabin,limits>::const_iterator itm;
   for (itv=v.begin(); itv!=v.end(); itv++) {
      for (itm=itv->begin(); itm!=itv->end(); itm++) {
         anabin thebin = itm->first;
         if (ans.find(thebin)==ans.end()) ans[thebin] = vector<limits>();
         ans[thebin].push_back(itm->second);
      }
   }

   return ans;
}

#endif // ifndef limits_h
