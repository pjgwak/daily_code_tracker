#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/LikelihoodInterval.h"
#include "RooStats/LikelihoodIntervalPlot.h"
#include "RooStats/FeldmanCousins.h"
#include "RooStats/ProfileLikelihoodTestStat.h"
#include "TLatex.h"
#include "TFile.h"

#include <iostream>
#include <fstream>
#include <algorithm>

int nworkers = 0;                        // number of worker for ProofLite (default use all available cores)
int nPointsToScan = 10;//20;
double additionalToysFac = 2.;//0.5;

pair<double,double> runLimit_RaaNS_Workspace_FC(const char *filename="TRIAL.root", const char *poiname="raa3", const char *pdfname="joint", const char *wsname="wcombo", const char* dataname="dOS_DATA", const char* ACTag="", double CI = 0.95, double poimin=0, double poimax=10);

// #include "StandardHypoTestInvDemo.C"

using namespace std;
using namespace RooStats;

pair<double,double> runLimit_RaaNS_Workspace_FC(const char *filename, const char *poiname, const char *pdfname, const char *wsname, const char* dataname, const char* ACTag, double CI, double poimin, double poimax)
{
   // RooMsgService::instance().addStream(RooFit::MsgLevel::DEBUG); 
   RooRealVar *theVar; RooDataSet *data; RooAbsPdf *pdf;
   TFile *f = new TFile(filename) ;
  
   // Open text file to write results
   ofstream of("results.txt");

   // Retrieve workspace from file
   RooWorkspace* ws = (RooWorkspace*) f->Get(wsname);
   RooStats::ModelConfig *sbHypo = (RooStats::ModelConfig*) ws->obj("SbHypo");
   RooStats::ModelConfig *bHypo = (RooStats::ModelConfig*) ws->obj("BHypo");
   cout << "observables: " << sbHypo->GetObservables()->getSize() << endl;
   theVar = ws->var(poiname);
   // if (!theVar) theVar = ws->function(poiname);
   pdf = ws->pdf(pdfname);
   data =(RooDataSet *) ws->data(dataname);

   // Print structure of composite p.d.f.
   pdf->Print("t") ;

   of << "Scanning " << filename << " using the POI " << poiname << ", the pdf " << pdfname << " and the workspace " << wsname << " with CL " << CI << endl;
   ws->loadSnapshot("SbHypo_poiAndNuisance_snapshot");
   of << poiname << " " << theVar->getVal() << endl;

   FeldmanCousins fc(*data,*sbHypo);
   fc.SetConfidenceLevel(CI);
   fc.AdditionalNToysFactor(additionalToysFac); // improve sampling that defines confidence belt
   //  fc.UseAdaptiveSampling(true); // speed it up a bit, but don't use for expectd limits
   RooDataSet poivals("poivals","poivals",RooArgSet(*theVar));
   double oldval = theVar->getVal();
   double dvar = (poimax-poimin)/nPointsToScan;
   cout << "TEST " << poimin << " " << poimax << endl;
   for (int i=0; i<nPointsToScan; i++) {
      theVar->setVal(poimin+i*dvar);
      poivals.add(RooArgSet(*theVar));
   }
   theVar->setVal(oldval);
   // fc.SetNBins(nPointsToScan); // set how many points per parameter of interest to scan
   // fc.SetPOIPointsToTest(poivals);
   // if (!sbHypo->GetNuisanceParameters()) fc.SetParameterPointsToTest(poivals);
   fc.SetPOIPointsToTest(poivals);
   fc.CreateConfBelt(true); // save the information in the belt for plotting

   ToyMCSampler*  toymcsampler = (ToyMCSampler*) fc.GetTestStatSampler();
   ProfileLikelihoodTestStat* testStat = dynamic_cast<ProfileLikelihoodTestStat*>(toymcsampler->GetTestStatistic());

   if(!sbHypo->GetPdf()->canBeExtended()){
      if(data->numEntries()==1)
         fc.FluctuateNumDataEntries(false);
      else
         cout <<"Not sure what to do about this model" <<endl;
   }

   ProofConfig pc(*ws, nworkers, "",false);
   toymcsampler->SetProofConfig(&pc); // enable proof

   if(sbHypo->GetGlobalObservables()){
      cout << "will use global observables for unconditional ensemble"<<endl;
      sbHypo->GetGlobalObservables()->Print();
      toymcsampler->SetGlobalObservables(*sbHypo->GetGlobalObservables());
   }


   // Now get the interval
   PointSetInterval* interval = fc.GetInterval();
   ConfidenceBelt* belt = fc.GetConfidenceBelt();


   pair<double,double> lims;
   lims.first = interval->LowerLimit(*theVar);
   if (lims.first == poimin) lims.first = -1e30;
   lims.second = interval->UpperLimit(*theVar);
   if (lims.second == poimax) lims.second = 1e30;

   cout << "The limits are [" << lims.first << ", " << lims.second << "]" << endl;
  
  return lims;
}
