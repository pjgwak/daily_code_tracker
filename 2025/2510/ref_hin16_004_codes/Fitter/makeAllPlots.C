#include "results2syst.C"
#include "plotVars.C"
#include "plotResults.C"
#include "Systematics/syst.h"
#include "Macros/drawMassPlot.C"

#include "TFile.h"
#include "TString.h"
#include "RooWorkspace.h"
#include "RooRealVar.h"

// flags
const bool doSysts = false;        // compute the systematics
const bool printSysts = false;     // print the systematics summary table
const bool plotMassPlots = false; 
const bool plotAllVars = true;     // plot the dependance of all vars with pt, centrality, y
const bool plotAllResults = false; // and plot the results!
const bool plotSystResults = false;// plot results taking syst fits as nominal too

void makeAllMassPlots(const char* workDirName, const char* DSTag="DATA", bool paperStyle=true, bool incSS=false);

void makeAllPlots(const char* nominalDir, const char* systDirsSig, const char* systDirsBkg) {
   string allDirs = string(nominalDir) + "," + string(systDirsSig) + "," + string(systDirsBkg);
   string allDirsSig = string(nominalDir) + "," + string(systDirsSig);
   string allDirsBkg = string(nominalDir) + "," + string(systDirsBkg);
   string systDirs = string(systDirsSig) + "," + string(systDirsBkg);

   // first, systematics
   if (doSysts) {
      results2syst(allDirsSig.c_str(), "syst_PP_fit_sig.csv", "Sig. shape (PP)", 1, "PP");
      results2syst(allDirsSig.c_str(), "syst_PbPb_fit_sig.csv", "Sig. shape (PbPb)", 1, "PbPb");
      results2syst(allDirsBkg.c_str(), "syst_PP_fit_bkg.csv", "Bkg. shape (PP)", 1, "PP");
      results2syst(allDirsBkg.c_str(), "syst_PbPb_fit_bkg.csv", "Bkg. shape (PbPb)", 1, "PbPb");
   }

   // print the systematics table
   // if (printSysts) {
      // readSyst_all("PP", "", true, "Systematics/tex/syst_PP_all.tex");
      // readSyst_all("PbPb", "", true, "Systematics/tex/syst_PbPb_all.tex");
   // }

   // draw mass plots
   if (plotMassPlots) {
      // cout << "Not yet implemented!" << endl;
      makeAllMassPlots(nominalDir);
   }

   // draw plots for all variables
   if (plotAllVars) {
      plotPt(nominalDir,"RFrac2Svs1S",""); 
      plotCent(nominalDir,"RFrac2Svs1S",""); 
      plotRap(nominalDir,"RFrac2Svs1S",""); 
      plotPt(nominalDir,"sigma1_Jpsi",""); 
      plotCent(nominalDir,"sigma1_Jpsi",""); 
      plotRap(nominalDir,"sigma1_Jpsi",""); 
      plotPt(nominalDir,"f_Jpsi",""); 
      plotCent(nominalDir,"f_Jpsi",""); 
      plotRap(nominalDir,"f_Jpsi",""); 
      plotPt(nominalDir,"rSigma21_Jpsi",""); 
      plotCent(nominalDir,"rSigma21_Jpsi",""); 
      plotRap(nominalDir,"rSigma21_Jpsi",""); 
      plotPt(nominalDir,"nparbkg","",false); 
      plotCent(nominalDir,"nparbkg","",false); 
      plotRap(nominalDir,"nparbkg","",false); 

      if (printSysts) {
         plotFiles(allDirs.c_str(), "RFrac2Svs1S", "pt", 0, 1.6, 6.5, 30, 0, 200, "PP");
         plotFiles(allDirs.c_str(), "RFrac2Svs1S", "pt", 1.6, 2.4, 3, 30, 0, 200, "PP");
         plotFiles(allDirs.c_str(), "RFrac2Svs1S", "pt", 0, 1.6, 6.5, 30, 0, 200, "PbPb");
         plotFiles(allDirs.c_str(), "RFrac2Svs1S", "pt", 1.6, 2.4, 3, 30, 0, 200, "PbPb");
         plotFiles(allDirs.c_str(), "RFrac2Svs1S", "cent", 0, 1.6, 6.5, 30, 0, 200, "PbPb");
         plotFiles(allDirs.c_str(), "RFrac2Svs1S", "cent", 1.6, 2.4, 3, 30, 0, 200, "PbPb");
         plotFiles(allDirs.c_str(), "normchi2", "pt", 0, 1.6, 6.5, 30, 0, 200, "PP",false);
         plotFiles(allDirs.c_str(), "normchi2", "pt", 1.6, 2.4, 3, 30, 0, 200, "PP",false);
         plotFiles(allDirs.c_str(), "normchi2", "pt", 0, 1.6, 6.5, 30, 0, 200, "PbPb",false);
         plotFiles(allDirs.c_str(), "normchi2", "pt", 1.6, 2.4, 3, 30, 0, 200, "PbPb",false);
         plotFiles(allDirs.c_str(), "normchi2", "cent", 0, 1.6, 6.5, 30, 0, 200, "PbPb",false);
         plotFiles(allDirs.c_str(), "normchi2", "cent", 1.6, 2.4, 3, 30, 0, 200, "PbPb",false);
         plotFiles(allDirs.c_str(), "chi2prob", "pt", 0, 1.6, 6.5, 30, 0, 200, "PP",false);
         plotFiles(allDirs.c_str(), "chi2prob", "pt", 1.6, 2.4, 3, 30, 0, 200, "PP",false);
         plotFiles(allDirs.c_str(), "chi2prob", "pt", 0, 1.6, 6.5, 30, 0, 200, "PbPb",false);
         plotFiles(allDirs.c_str(), "chi2prob", "pt", 1.6, 2.4, 3, 30, 0, 200, "PbPb",false);
         plotFiles(allDirs.c_str(), "chi2prob", "cent", 0, 1.6, 6.5, 30, 0, 200, "PbPb",false);
         plotFiles(allDirs.c_str(), "chi2prob", "cent", 1.6, 2.4, 3, 30, 0, 200, "PbPb",false);
      }
   }

   // plot results
   if (plotAllResults) {
      plotPt(nominalDir);
      plotCent(nominalDir);
   }
   if (plotSystResults) {
      TString workDirNamesStr(systDirs);
      TString workDirName; Int_t from = 0;
      while (workDirNamesStr.Tokenize(workDirName, from , ",")) {
         plotPt(workDirName.Data());
         plotCent(workDirName.Data());
      }
   }

   // we're done! and now, a script so that everything is in the right place
   cout << "We're done! Let's copy everything to the note now." << endl << endl << endl;
   cout << "### Copy the lines below to a script ###" << endl << endl;
   cout << "#!/bin/bash" << endl << endl;
   cout << "# settings" << endl;
   cout << "basedir=" << gSystem->pwd() << "/" << endl;
   cout << "notedir=/path/to/the/notes/notes/AN-16-067/trunk/" << endl << endl;

   // copy stuff...
   cout << "cp $basedir/Systematics/tex/*tex $notedir/tables" << endl;
   cout << "cp $basedir/Output/" << nominalDir << "/plot/RESULT/pdf/*pdf $notedir/figures/Results" << endl;
   cout << "cp $basedir/Output/" << nominalDir << "/tex/*tex $notedir/tables/Results" << endl;
   TString workDirNamesStr(systDirs);
   TString workDirName; Int_t from = 0;
   while (workDirNamesStr.Tokenize(workDirName, from , ",")) {
      cout << "mkdir -p $notedir/figures/Results/" << workDirName.Data() << " $notedir/tables/Results/" << workDirName.Data() << endl;
      cout << "cp $basedir/Output/" << workDirName.Data() << "/plot/RESULT/pdf/*pdf $notedir/figures/Results/" << workDirName.Data() << endl;
      cout << "cp $basedir/Output/" << workDirName.Data() << "/tex/*tex $notedir/tables/Results/" << workDirName.Data() << endl;
   }

   cout << endl;
   cout << "echo \"Done copying files. Now do:\"" << endl;
   cout << "echo \"cd  $notedir\"" << endl;
   cout << "echo \"cd ../../..\"" << endl;
   cout << "echo \"eval \\`notes/tdr runtime -sh\\`\"" << endl;
   cout << "echo \"cd -\"" << endl;
   cout << "echo \"tdr --style=an b AN-16-067\"" << endl;
   cout << "echo \"cp ../../tmp/AN-16-067_temp.pdf AN-16-067.pdf\"" << endl;
   cout << "echo \"svn commit\"" << endl;
   cout << "echo \"# freeze note; get approved; publish\"" << endl << endl;
   cout << "### Copy the lines above to a script ###" << endl;
}

void makeAllMassPlots(const char* workDirName, const char* DSTag, bool paperStyle, bool incSS) {
   // list of files
   vector<TString> theFiles = fileList(workDirName,"",DSTag);

   for (vector<TString>::const_iterator it=theFiles.begin(); it!=theFiles.end(); it++) {
      TFile *f = TFile::Open(*it);
      TString t; Int_t from = 0;
      bool catchjpsi=false, catchpsip=false, catchbkg=false;
      Char_t jpsiName[128], psipName[128], bkgName[128];
      while (it->Tokenize(t, from, "_")) {
         if (catchjpsi) {strcpy(jpsiName, t.Data()); catchjpsi=false;}
         if (catchpsip) {strcpy(psipName, t.Data()); catchpsip=false;}
         if (catchbkg) {strcpy(bkgName, t.Data()); catchbkg=false;}
         if (t=="Jpsi") catchjpsi=true;
         if (t=="Psi2S") catchpsip=true;
         if (t=="Bkg") catchbkg=true;
      }
      RooWorkspace *myws = (RooWorkspace*) f->Get("workspace");
      string outputDir = string("Output/") + string(workDirName) + string("/");
      struct InputOpt opt;
      opt.pp.RunNb.Start   = 262157; opt.PbPb.RunNb.Start = 262620;
      opt.pp.RunNb.End     = 262328; opt.PbPb.RunNb.End   = 263757;
      opt.pp.TriggerBit    = (int) PP::HLT_HIL1DoubleMu0_v1; 
      opt.PbPb.TriggerBit  = (int) HI::HLT_HIL1DoubleMu0_v1; 
      struct KinCuts cut;
      RooRealVar *invMass = myws->var("invMass");
      cut.dMuon.M.Min = invMass->getMin();
      cut.dMuon.M.Max = invMass->getMax();
      
      // if this is not a MB bin, skip it
      anabin thebin = binFromFile(*it);
      if (!(thebin==anabin(0,1.6,6.5,30,0,200) || thebin==anabin(1.6,2.4,3,30,0,200))) continue;

      cut.Centrality.Start = thebin.centbin().low();
      cut.Centrality.End = thebin.centbin().high();
      cut.dMuon.Pt.Min = thebin.ptbin().low();
      cut.dMuon.Pt.Max = thebin.ptbin().high();
      cut.dMuon.AbsRap.Min = thebin.rapbin().low();
      cut.dMuon.AbsRap.Max = thebin.rapbin().high();
      bool isPbPb = (it->Index("PbPb")>0);
      string plotLabel = "";
      bool incJpsi = true;//!paperStyle;
      bool incPsi2S = true;//!paperStyle;
      if (incJpsi)  { plotLabel = plotLabel + Form("_Jpsi_%s", jpsiName);  } 
      if (incPsi2S) { plotLabel = plotLabel + Form("_Psi2S_%s", psipName); }
      bool incBkg = true;
      bool cutCtau = true;
      bool doSimulFit = false;
      cout << *it << " " << it->Index("PbPb") << endl;
      bool setLogScale = !paperStyle;
      bool zoomPsi = paperStyle;
      int nBins = 46;
      bool getMeanPT = false;

      drawMassPlot(*myws, outputDir, opt, cut, plotLabel, DSTag, isPbPb, incJpsi, incPsi2S, incBkg, cutCtau, doSimulFit, false, setLogScale, incSS, zoomPsi, nBins, getMeanPT, paperStyle, false);

      delete myws;
      delete f;
   }
}
