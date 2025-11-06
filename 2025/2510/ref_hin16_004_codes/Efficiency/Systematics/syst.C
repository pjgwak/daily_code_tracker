#ifndef syst_C
#define syst_C

#include "syst.h"
#include "../../Fitter/Macros/Utilities/monster.h"

vector<TString> fileList_syst(const char* token, const char* prependPath) {
   vector<TString> ans;

   TString basedir("Systematics/csv/");
   if ( strcmp(prependPath,"") ) basedir.Prepend(Form("%s/",prependPath));
   TSystemDirectory dir(token,basedir);

   TList *files = dir.GetListOfFiles();

   if (files) {
      TIter next(files);
      TSystemFile *file;
      TString fname;

      while ((file=(TSystemFile*)next())) {
         fname = file->GetName();
         if (fname.EndsWith(".csv") && (TString(token) == "" || fname.Index(token) != kNPOS)) {
            ans.push_back(basedir+fname);
         }
      }
   }

   return ans;
};


map<anabin, syst> readSyst(const char* systfile, bool setSystZero) {
   map<anabin, syst> ans;

   ifstream file(systfile);
   if (!(file.good())) return ans;

   string systname; getline(file,systname);

   string line;
   double rapmin=0, rapmax=0, ptmin=0, ptmax=0, centmin=0, centmax=0, value=0;

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
         else if (cnt>6) {
            cout << "Warning, too many fields, I'll take the last one." << endl;
            continue;
         }
         cnt++;
      }
      anabin thebin(rapmin, rapmax, ptmin, ptmax, centmin, centmax);
      syst thesyst; thesyst.value = (setSystZero ? 0. : value); thesyst.name = systname;
      ans[thebin] = thesyst;
   }

   file.close();
  
   return ans;
};

map<anabin, syst> combineSyst(vector< map<anabin, syst> > theSysts, string name) {
   map<anabin, syst> ans;

   vector< map<anabin, syst> >::const_iterator it;
   for (it=theSysts.begin(); it!=theSysts.end(); it++) {
      map<anabin, syst>::const_iterator it2;
      for (it2=it->begin(); it2!=it->end(); it2++) {
         anabin thebin = it2->first;
         syst thesyst = it2->second;
         thesyst.name = name;

         // if we already have a syst for this bin, sum quadractically the existing syst and the new syst
         if (ans.find(thebin) != ans.end()) {
           thesyst.value = sqrt(pow(thesyst.value,2) + pow(ans[thebin].value,2));
         }
         ans[thebin] = thesyst;
      }
   }

   return ans;
};

map<anabin, syst> readSyst_all(const char* collSystem, const char* poiname, const char* nameTag, bool readEffSyst, bool readSigSyst, const char* prependPath, bool doPrintTex, const char* texName) {
   // poiname = NJpsi or BJpsi
   TString spoiname(poiname);
   TString sprependPath(prependPath);
   TString snameTag(nameTag);
  
   vector<TString> filelist;
  
   if (readEffSyst)
   {
     filelist.push_back(Form("%sSystematics/csv/stat_%s_%s_PbPb_acc.csv",sprependPath.Data(),snameTag.Data(),spoiname.Data()));
     filelist.push_back(Form("%sSystematics/csv/stat_%s_%s_PbPb_eff.csv",sprependPath.Data(),snameTag.Data(),spoiname.Data()));
     filelist.push_back(Form("%sSystematics/csv/stat_%s_%s_PbPb_muid_stat.csv",sprependPath.Data(),snameTag.Data(),spoiname.Data()));
     filelist.push_back(Form("%sSystematics/csv/stat_%s_%s_PbPb_sta_stat.csv",sprependPath.Data(),snameTag.Data(),spoiname.Data()));
     filelist.push_back(Form("%sSystematics/csv/stat_%s_%s_PbPb_trg_stat.csv",sprependPath.Data(),snameTag.Data(),spoiname.Data()));
     filelist.push_back(Form("%sSystematics/csv/stat_%s_%s_PbPb_ptWeighting.csv",sprependPath.Data(),snameTag.Data(),spoiname.Data()));
     filelist.push_back(Form("%sSystematics/csv/syst_%s_%s_PbPb_muid_syst.csv",sprependPath.Data(),snameTag.Data(),spoiname.Data()));
     filelist.push_back(Form("%sSystematics/csv/syst_%s_%s_PbPb_sta_syst.csv",sprependPath.Data(),snameTag.Data(),spoiname.Data()));
     filelist.push_back(Form("%sSystematics/csv/syst_%s_%s_PbPb_trg_syst.csv",sprependPath.Data(),snameTag.Data(),spoiname.Data()));
     filelist.push_back(Form("%sSystematics/csv/syst_%s_%s_PbPb_trg_binned.csv",sprependPath.Data(),snameTag.Data(),spoiname.Data()));
     
     filelist.push_back(Form("%sSystematics/csv/stat_%s_%s_PP_acc.csv",sprependPath.Data(),snameTag.Data(),spoiname.Data()));
     filelist.push_back(Form("%sSystematics/csv/stat_%s_%s_PP_eff.csv",sprependPath.Data(),snameTag.Data(),spoiname.Data()));
     filelist.push_back(Form("%sSystematics/csv/stat_%s_%s_PP_muid_stat.csv",sprependPath.Data(),snameTag.Data(),spoiname.Data()));
     filelist.push_back(Form("%sSystematics/csv/stat_%s_%s_PP_sta_stat.csv",sprependPath.Data(),snameTag.Data(),spoiname.Data()));
     filelist.push_back(Form("%sSystematics/csv/stat_%s_%s_PP_trg_stat.csv",sprependPath.Data(),snameTag.Data(),spoiname.Data()));
     filelist.push_back(Form("%sSystematics/csv/stat_%s_%s_PP_ptWeighting.csv",sprependPath.Data(),snameTag.Data(),spoiname.Data()));
     filelist.push_back(Form("%sSystematics/csv/syst_%s_%s_PP_muid_syst.csv",sprependPath.Data(),snameTag.Data(),spoiname.Data()));
     filelist.push_back(Form("%sSystematics/csv/syst_%s_%s_PP_sta_syst.csv",sprependPath.Data(),snameTag.Data(),spoiname.Data()));
     filelist.push_back(Form("%sSystematics/csv/syst_%s_%s_PP_trg_syst.csv",sprependPath.Data(),snameTag.Data(),spoiname.Data()));
     filelist.push_back(Form("%sSystematics/csv/syst_%s_%s_PP_trg_binned.csv",sprependPath.Data(),snameTag.Data(),spoiname.Data()));
   }
  
   vector< map<anabin, syst> > systmap_all;

   for (unsigned int i=0; i<filelist.size(); i++) {
      if (!filelist[i].Contains(collSystem)) continue;
      cout << filelist[i].Data() << endl;
      map<anabin,syst> systmap = readSyst(filelist[i].Data());
      systmap_all.push_back(systmap);
   }

   map<anabin,syst> ans = combineSyst(systmap_all, "Total");
   systmap_all.push_back(ans);

   if (doPrintTex) {
      printTex(systmap_all, TString(texName).ReplaceAll(".tex","_cent.tex").Data(), anabin(0,2.4,6.5,50,0,0));
      printTex(systmap_all, TString(texName).ReplaceAll(".tex","_pt.tex").Data(), anabin(0,2.4,0,0,0,200));
      printTex(systmap_all, TString(texName).ReplaceAll(".tex","_rap.tex").Data(), anabin(0,0,6.5,50,0,200));
      printTex(systmap_all, TString(texName).ReplaceAll(".tex","_ptrap.tex").Data(), anabin(0,-2.4,-6.5,-50,0,200));
      printTex(systmap_all, TString(texName).ReplaceAll(".tex","_centrap.tex").Data(), anabin(0,-2.4,6.5,50,0,-200));
      printTex(systmap_all, TString(texName).ReplaceAll(".tex","_ptcent.tex").Data(), anabin(0,2.4,-6.5,-50,0,-200));
      printTex(systmap_all, TString(texName).ReplaceAll(".tex","_fwd.tex").Data(), anabin(1.8,2.4,0,0,0,-200));
   }

   return ans;
};

void printTex(vector< map<anabin, syst> > theSysts, const char* texName, anabin mask) {
   // for the different components BIN (in rapidity bin, pt bin, centrality bin):
   // - if BIN = (0,0), then do not filter
   // - if BIN = (xx,yy), print ONLY bins that match these conditions
   // - if BIN = (-xx,-yy), print all bins EXCEPT those that match these conditions
   unsigned int nsyst = theSysts.size();

   ofstream file(texName);
   file << "\\begin{tabular}{|ccc|"; 
   for (unsigned int i=0; i<nsyst; i++) {
      if (i==nsyst-1) file << "|";
      file << "c|";
   }
   file << "}" << endl;
   file << "\\hline" << endl;
   file << "$|y|$ & \\pt & Centrality";
   for (unsigned int i=0; i<nsyst; i++) file << " & " << theSysts[i].begin()->second.name;
   file<< "\\\\" << endl;
   file << "\\hline" << endl;

   anabin oldbin(0,0,0,0,0,0);
   map<anabin, vector<syst> > themap = vm2mv(theSysts);
   map<anabin, vector<syst> >::const_iterator itm;
   for (itm=themap.begin(); itm!=themap.end(); itm++) {
      vector<syst> v = itm->second;
     if (v.size() != nsyst) {cout << v.size() << " ; " << nsyst << endl;
       for (auto e : v) cout << "v Name: " << e.name << " value: " << e.value << endl;
       for (auto e : theSysts) cout << "theSysts mapsize: " << e.size() << endl;
         cout << "Error, not all systematics have the same size. Aborting." << endl;
         file.close();
         return;
      }
      anabin thebin = itm->first;

      // filter
      if ((mask.rapbin().low()>0 || mask.rapbin().high()>0) && (thebin.rapbin() != mask.rapbin())) continue;
      if ((mask.rapbin().low()<0 || mask.rapbin().high()<0) && (thebin.rapbin() == binF(-mask.rapbin().low(),-mask.rapbin().high()))) continue;
      if ((mask.ptbin().low()>0 || mask.ptbin().high()>0) && (thebin.ptbin() != mask.ptbin())) continue;
      if ((mask.ptbin().low()<0 || mask.ptbin().high()<0) && (thebin.ptbin() == binF(-mask.ptbin().low(),-mask.ptbin().high()))) continue;
      if ((mask.centbin().low()>0 || mask.centbin().high()>0) && (thebin.centbin() != mask.centbin())) continue;
      if ((mask.centbin().low()<0 || mask.centbin().high()<0) && (thebin.centbin() == binI(-mask.centbin().low(),-mask.centbin().high()))) continue;

      if (thebin.rapbin() == oldbin.rapbin()) {
         file << " - ";
      } else {
         if (itm != themap.begin()) file << "\\hline" << endl;
         file.precision(1); file.setf(ios::fixed);
         file << thebin.rapbin().low() << " $< |y| < $ " << thebin.rapbin().high();
      }
      file << " & ";
      if (thebin.ptbin() == oldbin.ptbin()) {
         file << " - ";
      } else {
         file.precision(1); file.setf(ios::fixed);
         file << thebin.ptbin().low() << " $< \\pt < $ " << thebin.ptbin().high();
      }
      file << " & ";
      file.unsetf(ios::fixed);
      file << thebin.centbin().low()/2 << "-" << thebin.centbin().high()/2 << "\\% ";
      file.precision(1);
      file.setf(ios::fixed);

      for (unsigned int i=0; i<nsyst; i++) {
         file << " & " << 100.*v[i].value;
      }
      file << " \\\\" << endl;

      oldbin = thebin;
   }

   file << "\\hline" << endl;
   file << "\\end{tabular}" << endl;
   file.close();
   cout << "Closed " << texName << endl;
}

map<anabin, vector<syst> > vm2mv(vector< map<anabin,syst> > v) {
   map<anabin, vector<syst> > ans;

   vector<map<anabin,syst> >::const_iterator itv;
   map<anabin,syst>::const_iterator itm;
   for (itv=v.begin(); itv!=v.end(); itv++) {
      for (itm=itv->begin(); itm!=itv->end(); itm++) {
         anabin thebin = itm->first;
         if (ans.find(thebin)==ans.end()) ans[thebin] = vector<syst>();
         ans[thebin].push_back(itm->second);
      }
   }

   return ans;
}


RooWorkspace* getWorkspaceFromBin(anabin thebin, const char* workDirName, const char* token, const char* path2Fitter){
  
  // List of files in workDirName
  vector<TString> theFiles = fileList(workDirName,token,"DATA",path2Fitter);
  if ( (theFiles.size() < 1) )
  {
    cout << "#[Error]: No files found in " << workDirName << endl;
    return 0x0;
  }
  
  if (!strcmp(token,"PP"))
  {
    bin<int> bincent(0,200);
    thebin.setcentbin(bincent);
  }
  
  // Find the name of the file containing the workspace
  int cnt=1;
  bool wsFound(false);
  vector<TString>::const_iterator it;
  for (it=theFiles.begin(); it!=theFiles.end(); it++)
  {
    anabin thebinTest = binFromFile(*it);
    if (thebinTest==thebin)
    {
      wsFound = true;
      break;
    }
  }
  if (!wsFound)
  {
    cout << "No workspace root file found in " << workDirName << endl;
    return 0x0;
  }
  
  // Read file
  TFile *f = TFile::Open(*it,"READ");
  if (!f)
  {
    cout << "#[Error]: Unable to read file " << *it << endl;
    return 0x0;
  }
  
  // Retrieve workspace from file
  RooWorkspace* ws = (RooWorkspace*) f->Get("workspace");
  if (!ws)
  {
    cout << "#[Error]: Unable to retrieve workspace" << endl;
    return 0x0;
  }
  
  RooWorkspace* wsClone = (RooWorkspace*)(ws->Clone());
  f->Close(); delete f;
  return wsClone;
}

void expandsyst(const char* input, const char* output, bool ptdep, bool rapdep, bool centdep) {
   set<anabin> sb = allbins();
   map<anabin, syst> ms = readSyst(input);

   ofstream of(output);
   if (ms.size()>0) of << ms.begin()->second.name << endl;
   else {
      cout << "Error, empty input " << input << endl; 
      of.close(); 
      return;
   }

   for (set<anabin>::const_iterator it=sb.begin(); it!=sb.end(); it++) {
      double value=0;
      for (map<anabin,syst>::const_iterator it2=ms.begin(); it2!=ms.end(); it2++) {
         bool okpt = (!ptdep) || (it->ptbin() == it2->first.ptbin());
         bool okrap = (!rapdep) || (it->rapbin() == it2->first.rapbin());
         bool okcent = (!centdep) || (it->centbin() == it2->first.centbin());
         if (okpt && okrap && okcent) {
            value = it2->second.value;
            break;
         }
      }
      of << it->rapbin().low() << "," << it->rapbin().high() << ","
         << it->ptbin().low() << "," << it->ptbin().high() << ","
         << it->centbin().low() << "," << it->centbin().high() << ","
         << value << endl;
   }

   of.close();
}

#endif // ifndef syst_C
