#define printChi2Table_cxx

#include "printChi2Table.h"

#include <map>
#include <vector>
#include <iostream>
#include <fstream>

using namespace std;

struct chi2List {
  double chi2CB;
  double chi2CB2;
  double chi2CBG;
};

void printTable(map< string, map<anabin, chi2List> > chi2FullMap, const char* fileName);

void printChi2Table::Loop()
{
  //   In a ROOT session, you can do:
  //      Root > .L printChi2Table.C
  //      Root > printChi2Table t
  //      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
  
  if (fChain == 0) return;
  
  Long64_t nentries = fChain->GetEntriesFast();
  
  map< string, map<anabin, chi2List> > chi2FullMap;
  map<anabin, chi2List> chi2MapPP;
  map<anabin, chi2List> chi2MapPbPb;
  chi2FullMap["PP"] = chi2MapPP;
  chi2FullMap["PbPb"] = chi2MapPbPb;
  
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++)
  {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    
    map<anabin, chi2List> chi2Map = chi2FullMap[collSystem];
    
    chi2List list;
    
    anabin thebin(ymin, ymax, ptmin, ptmax, centmin, centmax);
    if (chi2Map.count(thebin)) list = chi2Map[thebin];
    
    if (!strcmp(jpsiName,"SingleCrystalBall")) list.chi2CB = normchi2;
    else if (!strcmp(jpsiName,"DoubleCrystalBall")) list.chi2CB2 = normchi2;
    else if (!strcmp(jpsiName,"GaussianAndCrystalBall")) list.chi2CBG = normchi2;
    else cout << "[ERROR] Unrecognized signal name" << endl;
      
    chi2Map[thebin] = list;
    chi2FullMap[collSystem] = chi2Map;
  }

  printTable(chi2FullMap,texName);
             
}

void printTable(map< string, map<anabin, chi2List> > chi2FullMap, const char* fileName)
{
  //string texName = "testTable.txt";
  ofstream file(fileName);
  file.precision(3);
  file << "\\begin{tabular}{|cccc|c|c|c|";
  file << "}" << endl;
  file << "\\hline" << endl;
  file << "$|y|$ & \\pt & Centrality & System";
  file << " & Single CB" << " & Double CB" << " & Gauss + CB" << " \\\\" << endl;
  file << "\\hline" << endl;
  
  int sysN = 2;
  const char* sysname[2] = {"PP", "PbPb"};
  
  for ( int i = 0 ; i < sysN ; i++ )
  {
    map<anabin, chi2List> chi2Map = chi2FullMap[sysname[i]];
   
    map<anabin, chi2List>::const_iterator itm;
    for (itm=chi2Map.begin(); itm!=chi2Map.end(); itm++)
    {
      anabin thebin = itm->first;
      chi2List list = itm->second;
      
      file << thebin.rapbin().low() << "-" << thebin.rapbin().high() << " & "
      << thebin.ptbin().low() << "-" << thebin.ptbin().high() << " & "
      << thebin.centbin().low()/2. << "-" << thebin.centbin().high()/2. << " & "
      << sysname[i] << " & ";
      
      file << list.chi2CB << " & " << list.chi2CB2 << " & " << list.chi2CBG;
      
      file << " \\\\" << endl;
    }
    file << "\\hline" << endl;
  }
  
  file << "\\end{tabular}" << endl;
  file.close();
  cout << "Closed " << fileName << endl;
}
