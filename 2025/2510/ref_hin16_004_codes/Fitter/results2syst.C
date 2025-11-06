#include "Macros/Utilities/bin.h"
#include "Macros/Utilities/resultUtils.h"
#include "Macros/Utilities/texUtils.h"
#include "results2tree.C"

#include "TString.h"
#include "RooRealVar.h"
#include "TMath.h"
#include "TSystem.h"

#include <map>
#include <vector>
#include <iostream>
#include <fstream>

using namespace std;

//#ifndef poiname_check
//#define poiname_check
//const char* poiname = "N_Jpsi";
//// const char* poiname = "eff";
//#endif

//////////////////
// DECLARATIONS //
//////////////////

// methods for computing the systematic
double rms(vector<double> v, bool isrelative);
double maxdiff(vector<double> v, bool isrelative);
void printTex(map<anabin, vector<double> > mapvals,
              map<anabin, vector<double> > mapchi2,
              map<anabin, vector<int> > mapndof,
              map<anabin, double> mapsyst,
              map<anabin, double> maperr,
              vector<string> vnames,
              string texName,
              anabin mask);

///////////////////
// MAIN FUNCTION //
///////////////////

void results2syst(const char* workDirNames, const char* systFileName, const char* systHeader, int method, const char* collTag="", bool relativeSyst=true, const char* poiname = "N_Jpsi", bool testChi2 = true, bool is16004 = false) {
  // workDirNames: of the form "dir1,dir2,dir3,..."
  // systFileName: "syst_blabla.csv" (do NOT specify the full path, it will be assigned automatically to Systematics/csv/)
  // systHeader:   this will be the header of the systematics file. A few words describing what this systematic is.
  // method:       0 -> RMS, 1 -> max difference to the first work dir (= nominal)
  // collTag:      can be "PP", "PbPb" or "" (for the ratio PbPb/PP)
  // relativeSyst: does the uncertainty scale with the central value? (false -> no, true -> yes)
  // poiname: it can be N_Jpsi, b_Jpsi, eff, effnp, acc, accnp
  
  map<anabin, vector<double> > mapvals, mapchi2;
  map<anabin, vector<int> > mapndof;
  map<anabin, double> maperr;
  TString workDirNamesStr(workDirNames);
  vector<string> vnames;
  TString workDirName; Int_t from = 0;
  int cnt=0;
  set<anabin> thebins;
  if (is16004) thebins = allbins16004();
  else thebins = allbins();
  
  TString spoiname(poiname);
  bool isPrompt = spoiname.Contains("_prompt");
  bool isNonPrompt = spoiname.Contains("_nonprompt");
  
  while (workDirNamesStr.Tokenize(workDirName, from , ",")) {
    vnames.push_back(workDirName.Data());
    for (set<anabin>::const_iterator it=thebins.begin(); it!=thebins.end(); it++) {
      
      anabin trbin = *it;
      if (string(collTag)=="PP" && trbin.centbin()!=binI(0,200)) continue;
      
      // if this is the first time we see this bin: create it
      if (mapvals.find(trbin)==mapvals.end()) {
        mapvals[trbin] = vector<double>();
        mapchi2[trbin] = vector<double>();
        mapndof[trbin] = vector<int>();
      }
      
      double val;
      double valB(-99.); // To store b-fraction in case we want prompt/nonprompt Jpsi systematics
      double chi2 = poiFromBin(workDirName.Data(),collTag,"chi2_invMass",trbin);
      double ndof = poiFromBin(workDirName.Data(),collTag,"ndof_invMass",trbin);
      
      // in the case of a really bad chi2, print a warning -- except if we are looking at efficiencies
      if (!(TString(poiname).Contains("eff")) && (ndof==-999 || TMath::Prob(chi2,ndof)<1e-10) && testChi2) {
        double ymin = trbin.rapbin().low();
        double ymax = trbin.rapbin().high();
        double ptmin = trbin.ptbin().low();
        double ptmax = trbin.ptbin().high();
        int centmin = trbin.centbin().low();
        int centmax = trbin.centbin().high();
        cout << "Warning, bad chi2: ";
        cout << " for " << workDirName.Data() << ", " << collTag << ", " << ymin << "<|y|<" << ymax << ", " << ptmin << "<pt<" << ptmax << ", " << centmin/2 << "-" << centmax/2 << "%. ";
        cout << "p(" << chi2 << "," << ndof << ")=" << TMath::Prob(chi2,ndof) << endl;
        mapvals[trbin].push_back(-999);
      } else {
        if (isPrompt || isNonPrompt)
        {
          valB = poiFromBin(workDirName.Data(),collTag,"b_Jpsi",trbin);
          val = poiFromBin(workDirName.Data(),collTag,"N_Jpsi",trbin);
          
          if (isPrompt) val = val*(1.-valB);
          if (isNonPrompt) val = val*valB;
        }
        else {val= poiFromBin(workDirName.Data(),collTag,poiname,trbin);}
        
        mapvals[trbin].push_back(val);
      }
      mapchi2[trbin].push_back(chi2);
      mapndof[trbin].push_back(ndof);
      if (cnt==0)
      {
        if (isPrompt || isNonPrompt)
        {
          double bNjpsiCorr = poiFromBin(workDirName.Data(),collTag,"correl_N_Jpsi_vs_b_Jpsi",trbin);
          double bErr = poiErrFromBin(workDirName.Data(),collTag,"b_Jpsi",trbin);
          double jpsi = poiFromBin(workDirName.Data(),collTag,"N_Jpsi",trbin);
          double jpsiErr = poiErrFromBin(workDirName.Data(),collTag,"N_Jpsi",trbin);
          
          double error(0.);
          if (isPrompt) {error = val*sqrt(pow(jpsiErr/jpsi,2)
                                          - 2.*bNjpsiCorr*jpsiErr*bErr/(jpsi*valB)
                                          + pow(bErr/valB,2));}
          if (isNonPrompt){error = val*sqrt(pow(jpsiErr/jpsi,2)
                                            + 2.*bNjpsiCorr*jpsiErr*bErr/(jpsi*valB)
                                            + pow(bErr/valB,2));}
          
          maperr[trbin] = error;
        }
        else maperr[trbin] = poiErrFromBin(workDirName.Data(),collTag,poiname,trbin);
      }
    }
    cnt++;
  }
  
  cout << "Done parsing files. Now computing the systematic." << endl;
  
  map<anabin,double> mapsyst;
  
  string fullSystFileName = string("Systematics/csv/") + systFileName;
  ofstream file(fullSystFileName.c_str());
  file << systHeader << endl;
  map<anabin, vector<double> >::const_iterator it;
  for (it=mapvals.begin(); it!=mapvals.end(); it++) {
    anabin thebin = it->first;
    vector<double> v = it->second;
    double value;
    if (method==0) value=rms(v, relativeSyst);
    else if (method==1) value=maxdiff(v, relativeSyst);
    else value=0;
    file << thebin.rapbin().low() << ", " << thebin.rapbin().high() << ", "
    << thebin.ptbin().low() << ", " << thebin.ptbin().high() << ", "
    << thebin.centbin().low() << ", " << thebin.centbin().high() << ", "
    << value << endl;
    mapsyst[thebin] = value;
  }
  file.close();
  cout << "Closed " << fullSystFileName << endl << endl;
  
  cout << "At last, let's summarize the differences in a TeX table." << endl;
  
  string texName(fullSystFileName);
  myReplace(texName,"csv","tex");
  gSystem->mkdir("Systematics/tex", kTRUE);
  
  printTex(mapvals, mapchi2, mapndof, mapsyst, maperr, vnames, TString(texName).ReplaceAll(".tex","_cent.tex").Data(), anabin(0,2.4,6.5,50,0,0));
  printTex(mapvals, mapchi2, mapndof, mapsyst, maperr, vnames, TString(texName).ReplaceAll(".tex","_pt.tex").Data(), anabin(0,2.4,0,0,0,200));
  printTex(mapvals, mapchi2, mapndof, mapsyst, maperr, vnames, TString(texName).ReplaceAll(".tex","_rap.tex").Data(), anabin(0,0,6.5,50,0,200));
  printTex(mapvals, mapchi2, mapndof, mapsyst, maperr, vnames, TString(texName).ReplaceAll(".tex","_ptrap.tex").Data(), anabin(0,-2.4,-6.5,-50,0,200));
  printTex(mapvals, mapchi2, mapndof, mapsyst, maperr, vnames, TString(texName).ReplaceAll(".tex","_centrap.tex").Data(), anabin(0,-2.4,6.5,50,0,-200));
  printTex(mapvals, mapchi2, mapndof, mapsyst, maperr, vnames, TString(texName).ReplaceAll(".tex","_ptcent.tex").Data(), anabin(0,2.4,-6.5,-50,0,-200));
  printTex(mapvals, mapchi2, mapndof, mapsyst, maperr, vnames, TString(texName).ReplaceAll(".tex","_fwd.tex").Data(), anabin(1.8,2.4,0,0,0,-200));
}

double rms(vector<double> v, bool isrelative) {
  if (v.size()==0 || v[0]==0) return 0;
  double s=0,s2=0;
  for (unsigned int i=0; i<v.size(); i++) {
    if (v[i]==-999) continue;
    s+=v[i];
    s2+=v[i]*v[i];
  }
  double ans = sqrt(s2-(s*s));
  if (isrelative) ans = ans/v[0];
  return ans;
}

double maxdiff(vector<double> v, bool isrelative) {
  if (v.size()==0 || v[0]==0) return 0;
  double maxdiff=0;
  for (unsigned int i=1; i<v.size(); i++) {
    if (v[i]==-999) continue;
    maxdiff=max(maxdiff,fabs(v[i]-v[0]));
  }
  double ans = maxdiff;
  if (isrelative) ans = ans/v[0];
  return ans;
}

void printTex(map<anabin, vector<double> > mapvals,
              map<anabin, vector<double> > mapchi2,
              map<anabin, vector<int> > mapndof,
              map<anabin, double> mapsyst,
              map<anabin, double> maperr,
              vector<string> vnames,
              string texName,
              anabin mask) {
  ofstream texfile(texName.c_str());
  texfile << "\\begin{tabular}{|ccc|";
  for (unsigned int i=0; i<vnames.size()+1; i++) {
    if (i==vnames.size()) texfile << "|";
    texfile << "c|";
  }
  texfile << "}" << endl;
  texfile << "\\hline" << endl;
  texfile << "$|y|$ & \\pt & Centrality";
  for (unsigned int i=0; i<vnames.size(); i++) texfile << " & " << latexSafe(vnames[i]);
  texfile << " & Systematic";
  texfile<< "\\\\" << endl;
  texfile << "\\hline" << endl;
  
  anabin oldbin(0,0,0,0,0,0);
  map<anabin, vector<double> >::const_iterator itm;
  for (itm=mapvals.begin(); itm!=mapvals.end(); itm++) {
    anabin thebin = itm->first;
    
    // filter
    if ((mask.rapbin().low()>0 || mask.rapbin().high()>0) && (thebin.rapbin() != mask.rapbin())) continue;
    if ((mask.rapbin().low()<0 || mask.rapbin().high()<0) && (thebin.rapbin() == binF(-mask.rapbin().low(),-mask.rapbin().high()))) continue;
    if ((mask.ptbin().low()>0 || mask.ptbin().high()>0) && (thebin.ptbin() != mask.ptbin())) continue;
    if ((mask.ptbin().low()<0 || mask.ptbin().high()<0) && (thebin.ptbin() == binF(-mask.ptbin().low(),-mask.ptbin().high()))) continue;
    if ((mask.centbin().low()>0 || mask.centbin().high()>0) && (thebin.centbin() != mask.centbin())) continue;
    if ((mask.centbin().low()<0 || mask.centbin().high()<0) && (thebin.centbin() == binI(-mask.centbin().low(),-mask.centbin().high()))) continue;
    
    vector<double> vval, vchi2; vector<int> vndof;
    vval = itm->second;
    vchi2 = mapchi2[thebin];
    vndof = mapndof[thebin];
    
    // skip problematic bins
    if (vval.size()==0 || vval[0]<=0 || vchi2[0]<=0 || vndof[0]<=0) continue;
    
    if (thebin.rapbin() == oldbin.rapbin()) {
      texfile << " - ";
    } else {
      if (itm != mapvals.begin()) texfile << "\\hline" << endl;
      texfile.precision(1); texfile.setf(ios::fixed);
      texfile << thebin.rapbin().low() << " $< |y| < $ " << thebin.rapbin().high();
    }
    texfile << " & ";
    if (thebin.ptbin() == oldbin.ptbin()) {
      texfile << " - ";
    } else {
      texfile.precision(1); texfile.setf(ios::fixed);
      texfile << thebin.ptbin().low() << " $< \\pt < $ " << thebin.ptbin().high();
    }
    texfile << " & ";
    texfile.unsetf(ios::fixed);
    texfile << thebin.centbin().low()/2 << "-" << thebin.centbin().high()/2 << "\\% ";
    texfile.precision(1);
    texfile.setf(ios::fixed);
    
    // first print the nominal value
    texfile << " & " << vval[0] << " (" << 100.*maperr[thebin]/vval[0] << "\\%, " << vchi2[0] << "/" << vndof[0] << ")";
    
    // then the alternative values
    for (unsigned int i=1; i<vval.size(); i++) {
      if (vval[i] != -999) texfile << " & " << 100.*(vval[i]-vval[0])/vval[0] << "\\% (" << vchi2[i] << "/" << vndof[i] << ")";
      else texfile << " & - ";
    }
    
    // and at last the systematic
    texfile << " & " << 100.*mapsyst[thebin] << "\\%";
    texfile << " \\\\" << endl;
    
    oldbin = thebin;
  }
  
  texfile << "\\hline" << endl;
  texfile << "\\end{tabular}" << endl;
  texfile.close();
  cout << "Closed " << texName << endl;
}
