// a macro that simply extracts the ratio of Psi2S/Jpsi from a list of files containing fit results

#include "TSystem.h"
#include "TSystemFile.h"
#include "TSystemDirectory.h"
#include "TROOT.h"
#include "TFile.h"
#include "TMath.h"
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "RooAbsData.h"
#include <iostream>
#include <fstream>
#include <dirent.h>
#include <vector>
#include <map>
#include <set>

using namespace std;

// let's define a few useful small class
struct model_t {
   string fileName;
   string binName;
   string modelName;
   int npar;
   double nll;
   int cnt;

   bool operator<(model_t other) const { // needed for std::set
      if (npar<other.npar) return true;
      if (npar>other.npar) return false;
      return fileName<other.fileName;
   }
};
typedef vector<model_t> vecModels_t;
typedef set<model_t> setModels_t;

vector<string> printNLL(map< string, setModels_t > content, string outputDir, string type, double pvalcut, map< string, vector<string> >& winnerModelNames) ;
void setLines(vector<string>& strLin, vector<string> lin); 
void printLines(vector<string> strLin, ofstream& fout); 
bool findFiles(string dirPath, vector<string>& fileNames); 
bool existDir(string dir);
void splitString(string stringOriginal, const string Key, string& stringWithoutKey, string& stringWithKey); 
bool readFiles(string dirPath, vector<string> fileNames, map<string, setModels_t>& content, string type);
bool extractNLL(string fileName, model_t& value);

bool keepLines(string InputFile, vector< int >& lineIndexToKeep, vector< vector< string > > winnerLabels);
bool reduceInputFile(string InputFile, string OutputFile, vector< int > lineIndexToKeep);
bool readInputFile(string FileName, vector< string >& content, int nRow);
void findSubDir(vector<string>& dirlist, string dirname, string ext="");
bool extractWinnerLabel(string fileName, string& winnerLabel);


void printLLRStudy(
                   string outputDirPath="Test",
                   string inputDirPath="",
                   string TAG = "DATA",
                   bool   cutSideBand = false,
                   string type = "Bkg",           // Type of the LLR test, available options are: "Bkg" , "Jpsi" and "Psi2S"
                   double pvalcut = 5.            // cut pvalue, in %  
                   ) 
{

  map< string, vector<string> > DIR;
  DIR["main"].push_back(gSystem->ExpandPathName(gSystem->pwd()));
  
  if (inputDirPath.find("/")==std::string::npos)  { inputDirPath = DIR["main"][0]+"/Input/"+outputDirPath+"/"; }
  if (outputDirPath.find("/")==std::string::npos) { outputDirPath = DIR["main"][0]+"/Output/"+outputDirPath+"/"; }

  DIR["output"].push_back(outputDirPath);
  DIR["outputLLR"].push_back(outputDirPath);
  DIR["inputLLR"].push_back(outputDirPath);
  if (existDir(DIR["output"][0])==false){ 
    cout << "[ERROR] Output directory: " << DIR["output"][0] << " does not exist!" << endl; 
    return;
  } else {
    findSubDir(DIR["output"], DIR["output"][0]);
    findSubDir(DIR["outputLLR"], DIR["outputLLR"][0], "LLR/");
    findSubDir(DIR["inputLLR"], DIR["inputLLR"][0], "LLR/Input/");
  }
  DIR["input"].push_back(inputDirPath);
  for(uint j = 1; j < DIR["output"].size(); j++) {
    string subdir = DIR["output"][j];
    subdir.replace(subdir.find(DIR["output"][0]), std::string(DIR["output"][0]).length(), DIR["input"][0]);
    if (existDir(subdir.c_str())==false){
      cout << "[ERROR] Input directory: " << subdir.c_str() << " does not exist!" << endl;
      return;
    }
    DIR["input"].push_back(subdir);
  }

  for(uint j = 0; j < DIR["output"].size(); j++) {
    if (DIR["output"].size()>1 && j==0) continue; // First entry is always the main input directory
    vector<string> fileNames;
    string dirPath = Form("%smass%s/%s/result/", DIR["output"][j].c_str(), (cutSideBand?"SB":""), TAG.c_str());
    if (!findFiles(dirPath, fileNames)) { return; } 
    cout << "[INFO] Creating " << ((type=="Bkg")?"Background":"Signal") << " Study summary!" << endl;
  
    // Group the files based on their background model
    map<string, setModels_t> content;
    if (!readFiles(dirPath, fileNames, content, type)) { return; }
   
    string plotDir = Form("%smass%s/%s/plot/", DIR["output"][j].c_str(), (cutSideBand?"SB":""), TAG.c_str());
    string outputDir = Form("%smass%s/%s/", DIR["outputLLR"][j].c_str(), (cutSideBand?"SB":""), TAG.c_str());
  
    if (existDir(outputDir)==false){ 
      cout << "[INFO] Output directory: " << outputDir << " does not exist, will create it!" << endl;
      if (existDir(outputDir)==false){ gSystem->mkdir(outputDir.c_str(), kTRUE); }
      if (existDir(DIR["inputLLR"][j].c_str())==false) { gSystem->mkdir(DIR["inputLLR"][j].c_str(), kTRUE);  }
      if (existDir(outputDir+"plot/pdf/")==false) { gSystem->mkdir((outputDir+"plot/pdf/").c_str(), kTRUE);  }
      if (existDir(outputDir+"plot/png/")==false) { gSystem->mkdir((outputDir+"plot/png/").c_str(), kTRUE);  }
      if (existDir(outputDir+"plot/root/")==false){ gSystem->mkdir((outputDir+"plot/root/").c_str(), kTRUE); }
      if (existDir(outputDir+"result/")==false){ gSystem->mkdir((outputDir+"result/").c_str(), kTRUE); }
      if (existDir(outputDir+"tex/")==false) { gSystem->mkdir((outputDir+"tex/").c_str(), kTRUE);  }
    }

    // Loop over each kinematic bin and compute the LLR/AIC tests
    map< string, vector< vector< string > > > winnerLabels; 
    map< string, vector< string > > winnerModelNames;
    vector<string> tmpVec; tmpVec.clear();
    winnerModelNames["PbPb"] = tmpVec; winnerModelNames["PP"] = tmpVec;
    map< string, int> i; i["PbPb"]=0; i["PP"]=0;
    vector<string> bestModelFiles = printNLL(content, outputDir, type, pvalcut, winnerModelNames); 
    cout << "[INFO] " << ((type=="Bkg")?"Background":"Signal") << " Study summary file done!" << endl; 
    
    cout << "The files for the best models are: " << endl;
    for (vector<string>::iterator it = bestModelFiles.begin(); it != bestModelFiles.end(); it++) {
      cout << *it << endl;
 
      string bestModelFile = *it;
      bestModelFile.erase(bestModelFile.find(".root"), bestModelFile.size());    

      gSystem->CopyFile((dirPath+bestModelFile+".root").c_str(), (outputDir+"result/"+bestModelFile+".root").c_str());

      bestModelFile.erase(0, bestModelFile.find("FIT_")+ string("FIT_").length());

      string outputFileName = bestModelFile;
      string tmp = bestModelFile; tmp.erase(0, tmp.find("Bkg_")+string("Bkg_").length()); tmp.erase(0, tmp.find("_"));
      outputFileName = outputFileName.erase(outputFileName.find("Bkg_")-1, outputFileName.size())+tmp;
     
      gSystem->CopyFile((plotDir+"pdf/PLOT_"+bestModelFile+".pdf").c_str(), (outputDir+"plot/pdf/PLOT_"+outputFileName+".pdf").c_str());
      gSystem->CopyFile((plotDir+"png/PLOT_"+bestModelFile+".png").c_str(), (outputDir+"plot/png/PLOT_"+outputFileName+".png").c_str());
      gSystem->CopyFile((plotDir+"root/PLOT_"+bestModelFile+".root").c_str(), (outputDir+"plot/root/PLOT_"+outputFileName+".root").c_str());

      string winnerLabel; bool isPbPb = false;
      if ( bestModelFile.find("PbPb")!=std::string::npos ) { isPbPb = true; }
      if (!extractWinnerLabel(dirPath+"FIT_"+bestModelFile+".root", winnerLabel)) { return; }
      vector < string > tmpLabels;
      tmpLabels.push_back(winnerLabel);
      tmpLabels.push_back( winnerModelNames[(isPbPb?"PbPb":"PP")].at(i[(isPbPb?"PbPb":"PP")]) );
      winnerLabels[(isPbPb?"PbPb":"PP")].push_back(tmpLabels);
      i[(isPbPb?"PbPb":"PP")]++;                        
    }

    // PRODUCING NEW INPUT FILES

    string InputFile;

    map<string, map<string, bool>> VARMAP = {
      {"MASS", 
       {
         {"BKG",   true}, 
         {"JPSI",  true}, 
         {"PSI2S", false}
       }
      },
      {"CTAU", 
       {
         {"BKG",   true}, 
         {"JPSI",  true}, 
         {"PSI2S", false},
         {"RES",   true},
         {"TRUE",  true},
       }
      }
    };
    map<string, bool> COLMAP = {{"PbPb", true}, {"PP", true}};
 
    map< string, vector< int > > lineIndexToKeep;
    InputFile = Form("InitialParam_MASS_BKG_%s.csv", "PbPb");
    keepLines(DIR["input"][j]+InputFile, lineIndexToKeep["PbPb"], winnerLabels["PbPb"]);
    InputFile = Form("InitialParam_MASS_BKG_%s.csv", "PP");
    keepLines(DIR["input"][j]+InputFile, lineIndexToKeep["PP"], winnerLabels["PP"]);

    typedef map<string, map<string, bool>>::iterator var_type;
    typedef map<string, bool>::iterator it_type;
    for(var_type VAR = VARMAP.begin(); VAR != VARMAP.end(); VAR++) {
      map<string, bool> PARMAP = VAR->second;
      string name1 = "InitialParam_" + VAR->first + "_";
      for(it_type PAR = PARMAP.begin(); PAR != PARMAP.end(); PAR++) {
        if (PAR->second) {
          string name2 = name1 + PAR->first + "_";
          for(it_type COL = COLMAP.begin(); COL != COLMAP.end(); COL++) {
            if(COL->second) {
              string name3 = name2 + COL->first + ".csv";
              string InputFile = (DIR["input"][j] + name3);
              string OutputFile = (DIR["inputLLR"][j] + name3);
              if (InputFile.find("PbPb")!=std::string::npos) {
                if (lineIndexToKeep.size()>0) reduceInputFile(InputFile, OutputFile, lineIndexToKeep["PbPb"]);
              } 
              else { if (lineIndexToKeep.size()>0) reduceInputFile(InputFile, OutputFile, lineIndexToKeep["PP"]); }
            }
          }
        }
      }
    }
  }

};


vector<string> printNLL(map< string, setModels_t > content, string outputDir, string type, double pvalcut, map< string, vector<string> >& winnerModelNames) 
{ 
  vector<string> ans;

  string outputFile = Form("%s/LLRTest_%s.txt", outputDir.c_str(), type.c_str());
  ofstream fout( outputFile );
  string outputFileTexTable = Form("%s/tex/LLRTest_%s_TexTables.txt", outputDir.c_str(), type.c_str());
  ofstream foutTexTable( outputFileTexTable );
  map< string, setModels_t>::iterator contIt;

  for ( contIt = content.begin(); contIt != content.end(); contIt++) {
      
    string      binName = contIt->first;
    setModels_t binCont = contIt->second;
    
    cout << "Analzing Kinematic Bin: " << binName << endl;
    cout << " " << endl;
    fout << "Analzing Kinematic Bin: " << binName << endl;
    fout << " " << endl;

    vecModels_t modelNLLB;
    setModels_t::iterator modelRow;

    for (modelRow = binCont.begin(); modelRow != binCont.end(); modelRow++) {  

      vector<string> strLin;
      unsigned int iB=0;

      setModels_t::iterator modelCol;
      for (modelCol = binCont.begin(); modelCol != binCont.end(); modelCol++) {

        string  modelNameA = modelCol->modelName;
        int     nParA      = modelCol->npar;
        double  NLLA       = modelCol->nll;
        double  AICA       = 2*(nParA + NLLA);

        if (modelRow==binCont.begin()) {
          modelNLLB.push_back(*modelCol);
        }

        vector<string> lin;

        if (nParA>=modelRow->npar) {
          string  modelNameB = modelNLLB[iB].modelName;
          int     nParB      = modelNLLB[iB].npar;
          double  NLLB       = modelNLLB[iB].nll;
          double  AICB       = 2*(nParB + NLLB);
          if (modelNameA==modelNameB) {
            lin.push_back("|------------------------");
            lin.push_back("| "+modelNameA);
            lin.push_back(Form("|    NLL: %.2f  ", NLLA));
            lin.push_back(Form("|    AIC: %.2f  ", AICA));
            lin.push_back("|------------------------");
          } else if (nParA >= nParB) {
            double  diffNLL    = -2.0*(NLLA - NLLB);
            double  diffNPar   =  2.0*(nParA-nParB);
            double  probChi2   = 100.*TMath::Prob(diffNLL, diffNPar);
            if (diffNLL<0) probChi2 = 100.;
            if (probChi2>pvalcut && (nParA-nParB)<=2) modelNLLB[iB].cnt++;
            lin.push_back("| "+modelNameA);
            lin.push_back(Form("|    NLL: %.2f  ", NLLA));
            lin.push_back(Form("|    Diff: %.2f  ", diffNLL));
            lin.push_back(Form("|    Prob: %.1f%s   ", probChi2, "%"));
            lin.push_back(Form("|    AIC: %.2f  ", -(AICA-AICB)));
            lin.push_back("|------------------------");
          } 
          setLines(strLin, lin); iB++;
        } 
      }
      for (unsigned int j=0; j<strLin.size(); j++) { strLin[j] = strLin[j] + "|"; }
      printLines(strLin, fout);
    }      

    // which is the best model for this bin?
    string bestModelFile="NOTFOUND", bestModelName=""; int minok=999; int maxpar=0;
    for (vecModels_t::iterator it=modelNLLB.begin(); it!=modelNLLB.end();it++) {
       if (it->npar > maxpar) maxpar = it->npar;
       if (it->cnt>=2 && it->npar<minok) {
          bestModelFile=it->fileName;
          bestModelName=it->modelName;
          minok = it->npar;
       }
    }
    if (minok==999) { // sometimes the best model is one of the two highest orders...
       for (vecModels_t::iterator it=modelNLLB.begin(); it!=modelNLLB.end();it++) {
          int npar = it->npar;
          if (it->cnt>=maxpar-npar && npar<minok) {
             bestModelFile=it->fileName+" WARNING, HIGH ORDER";
             bestModelName=it->modelName;
             minok = it->npar;
          }
       }
    }

    cout << endl << " And the winner is... " << bestModelFile << endl << endl << endl;
    fout << endl << " And the winner is... " << bestModelFile << endl << endl << endl;
    ans.push_back(bestModelFile);
    bool isPbPb = false;
    if (bestModelFile.find("PbPb")!=std::string::npos) isPbPb=true;
    winnerModelNames[(isPbPb?"PbPb":"PP")].push_back(bestModelName);

    vector<string> TexTable;
    TexTable.push_back("\\begin{table}");
    TexTable.push_back("\\centering");
    string iniLinLatex = "\\begin{tabular}{ c  c";
    string header = "N & NLL ";
    for (unsigned int iM=0; iM<binCont.size(); iM++) { 
      if ((binCont.size()-1-iM)>=2) {
        iniLinLatex = iniLinLatex + " c ";
        header = header + "& " + Form("p(H0: N = %d)",iM) + " ";
      }
    }
    iniLinLatex = iniLinLatex + "}";
    header = header + "\\\\";
    TexTable.push_back(iniLinLatex);
    TexTable.push_back(header);
    TexTable.push_back("\\hline");
    
    vector<string> strLinLatex;
    unsigned int iR=0;
    int nParMax = 0;
    setModels_t::iterator model;
    for (model = binCont.begin(); model != binCont.end(); model++) { 
      if (model->npar > nParMax) { nParMax = model->npar; }
    }
    for (modelRow = binCont.begin(); modelRow != binCont.end(); modelRow++) {  

      unsigned int iA=-1;

      setModels_t::iterator modelCol;
      for (modelCol = binCont.begin(); modelCol != binCont.end(); modelCol++) {

        string  modelNameA = modelCol->modelName;
        int     nParA      = modelCol->npar;
        double  NLLA       = modelCol->nll;
        double  AICA       = 2*(nParA + NLLA);

        if (modelRow==binCont.begin()) {
          modelNLLB.push_back(*modelCol);
          strLinLatex.push_back( Form("%d & %.2f", ++iA, NLLA) );
        }

        if (nParA<modelRow->npar) {
          int     nParRow = modelRow->npar;
          double  NLLRow  = modelRow->nll;
          if ((nParMax - nParA)>=2) {
            double  diffNLL    = -2.0*(NLLRow - NLLA);
            double  diffNPar   =  2.0*(nParRow-nParA);
            double  probChi2   = 100.*TMath::Prob(diffNLL, diffNPar);
            if (diffNLL<0) probChi2 = 100.;
            if (bestModelFile==modelCol->fileName && (nParRow-nParA)<=2) {
              strLinLatex[iR] = strLinLatex[iR] + " & " + Form("\\textbf{%.1f%s}", probChi2, "$\\%$");
            } else {
              strLinLatex[iR] = strLinLatex[iR] + " & " + Form("%.1f%s", probChi2, "$\\%$");
            } 
          }
        } else {
          if ((nParMax - nParA)>=2) {
            strLinLatex[iR] = strLinLatex[iR] + " & ";
          }
        } 
      }
      iR++;
    }      
    for (unsigned int j=0; j<strLinLatex.size(); j++) {  
      strLinLatex[j] = strLinLatex[j] + "\\\\";
      TexTable.push_back(strLinLatex[j]);
    }
    TexTable.push_back("\\end{tabular}");
    TexTable.push_back(Form("\\label{tab:LLRTEST_%s}", binName.c_str()));
    string rapStr, centStr, ptStr, modelStr, colStr;
    if (bestModelFile.find("ExpChebychev")!=std::string::npos){ modelStr = "exponential chebychev polynomials"; }
    else if (bestModelFile.find("Chebychev")!=std::string::npos){ modelStr = "chebychev polynomials"; }
    else { modelStr = "chebychev polynomials"; }
    
    if (binName.find("rap06")!=std::string::npos){ rapStr = "$\\abs{y} <$ 0.6"; }
    else if (binName.find("rap024")!=std::string::npos){ rapStr = "$\\abs{y} <$ 2.4"; }
    else if (binName.find("rap0612")!=std::string::npos){ rapStr = "0.6 $\\leq \\abs{y} <$ 1.2"; }
    else if (binName.find("rap1218")!=std::string::npos){ rapStr = "1.2 $\\leq \\abs{y} <$ 1.8"; }
    else if (binName.find("rap1824")!=std::string::npos){ rapStr = "1.8 $\\leq \\abs{y} <$ 2.4"; }
    else if (binName.find("rap0004")!=std::string::npos){ rapStr = "$\\abs{y} <$ 0.4"; }
    else if (binName.find("rap0408")!=std::string::npos){ rapStr = "0.4 $\\leq \\abs{y} <$ 0.8"; }
    else if (binName.find("rap0812")!=std::string::npos){ rapStr = "0.8 $\\leq \\abs{y} <$ 1.2"; }
    else if (binName.find("rap1216")!=std::string::npos){ rapStr = "1.2 $\\leq \\abs{y} <$ 1.6"; }
    else if (binName.find("rap1620")!=std::string::npos){ rapStr = "1.6 $\\leq \\abs{y} <$ 2.0"; }
    else if (binName.find("rap2024")!=std::string::npos){ rapStr = "2.0 $\\leq \\abs{y} <$ 2.4"; }
    
    if (binName.find("pt30500")!=std::string::npos){ ptStr = "3.0 $\\leq \\PT <$ 50.0 $\\GeVc$"; }
    else if (binName.find("pt65500")!=std::string::npos){ ptStr = "6.5 $\\leq \\PT <$ 50.0 $\\GeVc$"; }
    else if (binName.find("pt30300")!=std::string::npos){ ptStr = "3.0 $\\leq \\PT <$ 30.0 $\\GeVc$"; }
    else if (binName.find("pt3045")!=std::string::npos){ ptStr = "3.0 $\\leq \\PT <$ 4.5 $\\GeVc$"; }
    else if (binName.find("pt3065")!=std::string::npos){ ptStr = "3.0 $\\leq \\PT <$ 6.5 $\\GeVc$"; }
    else if (binName.find("pt4555")!=std::string::npos){ ptStr = "4.5 $\\leq \\PT <$ 5.5 $\\GeVc$"; }
    else if (binName.find("pt5565")!=std::string::npos){ ptStr = "5.5 $\\leq \\PT <$ 6.5 $\\GeVc$"; }
    else if (binName.find("pt6575")!=std::string::npos){ ptStr = "6.5 $\\leq \\PT <$ 7.5 $\\GeVc$"; }
    else if (binName.find("pt6585")!=std::string::npos){ ptStr = "6.5 $\\leq \\PT <$ 8.5 $\\GeVc$"; }
    else if (binName.find("pt7585")!=std::string::npos){ ptStr = "7.5 $\\leq \\PT <$ 8.5 $\\GeVc$"; }
    else if (binName.find("pt8595")!=std::string::npos){ ptStr = "8.5 $\\leq \\PT <$ 9.5 $\\GeVc$"; }
    else if (binName.find("pt95110")!=std::string::npos){ ptStr = "9.5 $\\leq \\PT <$ 11.0 $\\GeVc$"; }
    else if (binName.find("pt110130")!=std::string::npos){ ptStr = "11.0 $\\leq \\PT <$ 13.0 $\\GeVc$"; }
    else if (binName.find("pt110150")!=std::string::npos){ ptStr = "11.0 $\\leq \\PT <$ 15.0 $\\GeVc$"; }
    else if (binName.find("pt130150")!=std::string::npos){ ptStr = "13.0 $\\leq \\PT <$ 15.0 $\\GeVc$"; }
    else if (binName.find("pt150200")!=std::string::npos){ ptStr = "15.0 $\\leq \\PT <$ 20.0 $\\GeVc$"; }
    else if (binName.find("pt150500")!=std::string::npos){ ptStr = "15.0 $\\leq \\PT <$ 50.0 $\\GeVc$"; }
    else if (binName.find("pt200300")!=std::string::npos){ ptStr = "20.0 $\\leq \\PT <$ 30.0 $\\GeVc$"; }
    else if (binName.find("pt200500")!=std::string::npos){ ptStr = "20.0 $\\leq \\PT <$ 50.0 $\\GeVc$"; }
    else if (binName.find("pt300500")!=std::string::npos){ ptStr = "30.0 $\\leq \\PT <$ 50.0 $\\GeVc$"; }
    
    if (binName.find("cent0200")!=std::string::npos){ centStr = "centratility bin 0-100$\\%$"; }
    else if (binName.find("cent020")!=std::string::npos){ centStr = "centratility bin 0-10$\\%$"; }
    else if (binName.find("cent010")!=std::string::npos){ centStr = "centratility bin 0-5$\\%$"; }
    else if (binName.find("cent1020")!=std::string::npos){ centStr = "centratility bin 5-10$\\%$"; }
    else if (binName.find("cent2030")!=std::string::npos){ centStr = "centratility bin 10-15$\\%$"; }
    else if (binName.find("cent3040")!=std::string::npos){ centStr = "centratility bin 15-20$\\%$"; }
    else if (binName.find("cent4050")!=std::string::npos){ centStr = "centratility bin 20-25$\\%$"; }
    else if (binName.find("cent5060")!=std::string::npos){ centStr = "centratility bin 25-30$\\%$"; }
    else if (binName.find("cent6070")!=std::string::npos){ centStr = "centratility bin 30-35$\\%$"; }
    else if (binName.find("cent7080")!=std::string::npos){ centStr = "centratility bin 35-40$\\%$"; }
    else if (binName.find("cent8090")!=std::string::npos){ centStr = "centratility bin 40-45$\\%$"; }
    else if (binName.find("cent90100")!=std::string::npos){ centStr = "centratility bin 45-50$\\%$"; }
    else if (binName.find("cent100120")!=std::string::npos){ centStr = "centratility bin 50-60$\\%$"; }
    else if (binName.find("cent120140")!=std::string::npos){ centStr = "centratility bin 60-70$\\%$"; }
    else if (binName.find("cent140200")!=std::string::npos){ centStr = "centratility bin 70-100$\\%$"; }
    else if (binName.find("cent2040")!=std::string::npos){ centStr = "centratility bin 10-20$\\%$"; }
    else if (binName.find("cent4060")!=std::string::npos){ centStr = "centratility bin 20-30$\\%$"; }
    else if (binName.find("cent6080")!=std::string::npos){ centStr = "centratility bin 30-40$\\%$"; }
    else if (binName.find("cent80100")!=std::string::npos){ centStr = "centratility bin 40-50$\\%$"; }
    else if (binName.find("cent100200")!=std::string::npos){ centStr = "centratility bin 50-1000$\\%$"; }
    else if (binName.find("cent2060")!=std::string::npos){ centStr = "centratility bin 10-30$\\%$"; }
    else if (binName.find("cent60200")!=std::string::npos){ centStr = "centratility bin 30-100$\\%$"; }
    else if (binName.find("cent040")!=std::string::npos){ centStr = "centratility bin 0-20$\\%$"; }
    else if (binName.find("cent4080")!=std::string::npos){ centStr = "centratility bin 20-40$\\%$"; }
    else if (binName.find("cent80200")!=std::string::npos){ centStr = "centratility bin 40-100$\\%$"; }

    if (binName.find("PP")!=std::string::npos){ colStr = "pp"; }
    else if (binName.find("PbPb")!=std::string::npos){ colStr = "PbPb"; }
    TexTable.push_back(Form("\\caption{Negative loglikelihoods for fits with %s of orders 0-6 of %s data in %s %s. In addition the p-values of the LLR-test for the null-hypothesis are listed. Tests of which the null-hypothesis cannot be rejected for two consecutive orders are highlighted in bold, together with the corresponding order.}", modelStr.c_str(), colStr.c_str(), rapStr.c_str(), (colStr=="pp" ? Form("and %s", ptStr.c_str()) : Form(", %s and %s", ptStr.c_str(), centStr.c_str()))));
    TexTable.push_back("\\end{table}");
    printLines(TexTable, foutTexTable);
    foutTexTable << endl; foutTexTable << endl;

  } // bin loop

  return ans;
};


void setLines(vector<string>& strLin, vector<string> lin) 
{
  const  string  empty  = "                          ";
  if (strLin.size() < lin.size()) {
    strLin = lin;  
  } else {
    for (unsigned int i = 0; i < lin.size(); i++) {
      strLin.at(i) = strLin.at(i) + lin.at(i);
    }
  }
  for (unsigned int i = 0; i < lin.size(); i++) {
    strLin.at(i).append(empty, 0, (25-lin.at(i).length()));
  }
    
  return;
};


void printLines(vector<string> strLin, ofstream& fout) 
{
  for (vector<string>::iterator line = strLin.begin(); line < strLin.end(); line++) {
      cout << *line << endl;
      fout << *line << endl;
  }
};


bool extractNLL(string fileName, model_t& value) 
{
  string binName = fileName;
  binName.erase(0, binName.find("FIT_")+ string("FIT_").length());

  value.npar = -1;
  value.nll = 1e99;
  TFile *f = new TFile( fileName.c_str() );
  if (!f) {
    cout << "[Error] " << fileName << " not found" << endl; return false;
  }
  RooWorkspace *ws = (RooWorkspace*) f->Get("workspace");
  if (!ws) {
    f->Close(); delete f;
    cout << "[ERROR] Workspace not found in " << fileName << endl; return false;
  }    
  RooAbsReal *nll = NULL;
  double NLL = 0;
  int npar = 0;
  string dsName = ws->allData().front()->GetName();

  if (binName.find("PP")!=std::string::npos) {
    nll = ws->pdf("pdfMASS_Tot_PP")->createNLL(*ws->data(dsName.c_str()));
    npar = ws->pdf("pdfMASS_Bkg_PP")->getParameters(*ws->data(dsName.c_str()))->getSize();
  }
  if (binName.find("PbPb")!=std::string::npos) {
    nll = ws->pdf("pdfMASS_Tot_PbPb")->createNLL(*ws->data(dsName.c_str()));
    npar = ws->pdf("pdfMASS_Bkg_PbPb")->getParameters(*ws->data(dsName.c_str()))->getSize();
  }
 
  if (!nll) {
    delete ws;
    f->Close(); delete f;
    cout << "[ERROR] NLL was not found" << endl; return false;
  } else {
    NLL = nll->getVal();
  }

  value.npar = npar;
  value.nll = NLL;

  delete nll; delete ws;
  f->Close(); delete f;
        
  return true;
};


bool existDir(string dir)
{
  bool exist = false;
  void * dirp = gSystem->OpenDirectory(dir.c_str());
  if (dirp){
    gSystem->FreeDirectory(dirp);
    exist = true;
  }
  return exist;
};


bool readFiles(string dirPath, vector<string> fileNames, map<string, setModels_t>& content, string type)
{
  for (vector<string>::iterator it = fileNames.begin(); it < fileNames.end(); it++) {
    string fileName = *it;

    string binName;   
    string modelName; 
    splitString(fileName, type, binName, modelName); 
    
    model_t modelNLL;
    modelNLL.binName = binName;
    modelNLL.modelName = modelName;
    modelNLL.fileName = fileName;
    modelNLL.cnt=0;
    
    if (extractNLL(dirPath+fileName, modelNLL) && modelName!="") {
      if (content.find(binName) == content.end()) content[binName] = setModels_t();
      content[binName].insert(modelNLL);
    }
  }
  if (content.size()==0) {
    cout << "[ERROR] No NLL values were found in the input files" << endl; return false;
  }
  return true;
};

 
void splitString(string stringOriginal, const string Key, string& stringWithoutKey, string& stringWithKey) 
{
  string tmp  = stringOriginal;   
  
  if (tmp.find(Key+"_")==std::string::npos) {
    stringWithoutKey = "";
    stringWithKey = stringOriginal;
    return;
  }
  

  tmp.erase(0, tmp.find(Key+"_")+(Key+"_").length()); stringWithKey = tmp;
  stringWithKey.erase(stringWithKey.find("_"), stringWithKey.size());
      
  tmp.erase(0, tmp.find("_"));
  stringWithoutKey = stringOriginal;
  stringWithoutKey = stringWithoutKey.erase(stringWithoutKey.find(Key+"_")-1, stringWithoutKey.size())+tmp.erase(tmp.find(".root"), tmp.size());
};


bool findFiles(string dirPath, vector<string>& fileNames) 
{

  DIR *dpdf;
  struct dirent *epdf;

  // Open the working directory
  dpdf = opendir(dirPath.c_str());
  // Search for all the files inside this directory
  if (dpdf != NULL){
    while ( (epdf=readdir(dpdf)) ){
      if (strcmp(epdf->d_name,".")!=0 && strcmp(epdf->d_name,"..")!=0 ) {
        std::cout << "[INFO] Adding file: " << epdf->d_name << std::endl;
        string filePath = epdf->d_name; 
        fileNames.push_back(filePath);
      }
    }
  } else {
    cout << "[ERROR] Working directory was not found!" << endl; return false;
  }
  return true;

};
        
        
///////////////////////////////////////////////////////////////

bool keepLines(string InputFile, vector< int >& lineIndexToKeep, vector< vector< string > > winnerLabels)
{
  vector< string > inputContent;
  bool isPbPb = false;
  if ( InputFile.find("PbPb")!=std::string::npos ) isPbPb = true;
  if(!readInputFile(InputFile, inputContent, -1)){ return false; }

  int i = 0;
  for(vector< string >::iterator iRow=inputContent.begin(); iRow!=inputContent.end(); ++iRow) {
    string row = *iRow;
    if ( i==0 ) { lineIndexToKeep.push_back(i); }
    else {
      bool found = false;
      for(vector< vector< string > >::iterator iLine=winnerLabels.begin(); iLine!=winnerLabels.end(); ++iLine) {
        string kinematic = (*iLine)[0]; if (isPbPb==false) { kinematic.erase(kinematic.find("-0-500;"),string("-0-100;").length()); }
        string modelName = (*iLine)[1];
        if ( row.find(kinematic)!=std::string::npos && row.find(modelName)!=std::string::npos ) { found=true; break; }
      }
      if (found==true) { lineIndexToKeep.push_back(i); }
    }
    i++;
  }
  return true;
};

bool reduceInputFile(string InputFile, string OutputFile, vector< int > lineIndexToKeep)
{

  vector< string > inputContent, outputContent; 
  if(!readInputFile(InputFile, inputContent, -1)){ return false; }
  
  for(vector< int >::iterator index=lineIndexToKeep.begin(); index!=lineIndexToKeep.end(); ++index) {
    int i = *index;
    outputContent.push_back(inputContent[i]);
  }

  if (outputContent.size()>0) {
    ofstream myfile(OutputFile.c_str());
    if (myfile.is_open()){
      for (vector<string>::iterator line = outputContent.begin(); line < outputContent.end(); line++) {
        myfile << *line << endl;
      }
    }
  }

  return true;
};


bool readInputFile(string FileName, vector< string >& content, int nRow)
{
  if (nRow==0) { 
    cout << "[WARNING] Ignoring content of File: " << FileName << endl; return true; 
  }
  if (nRow!=1) { cout << "[INFO] Checking file: " << FileName << endl; }
  ifstream myfile(FileName.c_str());
  if (myfile.is_open()){ 
    string line;
    while ( getline(myfile, line) ){
      if (nRow==0) break; else {nRow=nRow-1;}
      content.push_back(line);
    }
  } else {
    cout << "[WARNING] File: " << FileName << " was not found!" << endl; return false;
  }
  return true;
};


void findSubDir(vector<string>& dirlist, string dirname, string ext)
{
  TSystemDirectory dir(dirname.c_str(), dirname.c_str());
  TList *subdirs = dir.GetListOfFiles();
  if (subdirs) {
    TSystemFile *subdir;
    TIter next(subdirs);
    int countDirs(0);
    while ((subdir=(TSystemFile*)next())) {
      if (subdir->IsDirectory() && string(subdir->GetName())!="." && string(subdir->GetName())!=".." && string(subdir->GetName()).find("LLR")==std::string::npos && string(subdir->GetName()).find("mass")==std::string::npos && string(subdir->GetName()).find("ctau")==std::string::npos) {
        dirlist.push_back(dirname + ext + subdir->GetName() + "/");
        cout << "[INFO] Input subdirectory: " << dirname + ext + subdir->GetName() + "/" << " found!" << endl;
        countDirs++;
      }
    }
    if (countDirs==0){
      dirlist.push_back(dirname + ext);
      cout << "[INFO] Input subdirectory: " << dirname + ext << " found!" << endl;
    }
  }
  delete subdirs;
  return;
};


bool extractWinnerLabel(string fileName, string& winnerLabel) 
{
  TFile *f = new TFile( fileName.c_str() );
  if (!f) {
    cout << "[Error] " << fileName << " not found" << endl; return false;
  }
  RooWorkspace *ws = (RooWorkspace*) f->Get("workspace");
  if (!ws) {
    f->Close(); delete f;
    cout << "[ERROR] Workspace not found in " << fileName << endl; return false;
  }
  
  winnerLabel = Form("%.1f-%.1f;%.1f-%.1f;%.0f-%.0f;",
                     ws->var("rap")->getMin(), ws->var("rap")->getMax(),
                     ws->var("pt")->getMin(), ws->var("pt")->getMax(),
                     ws->var("cent")->getMin()/2.0, ws->var("cent")->getMax()/2.0
                     );

  delete ws;
  f->Close(); delete f;
        
  return true;
};
