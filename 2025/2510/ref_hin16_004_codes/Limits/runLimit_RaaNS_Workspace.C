#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/LikelihoodInterval.h"
#include "RooStats/LikelihoodIntervalPlot.h"
#include "TLatex.h"
#include <utility>

using namespace RooStats;
using namespace std;

pair<double,double> runLimit_RaaNS_Workspace(const char *filename="TRIAL.root", const char *poiname="raa3", const char *pdfname="joint", const char *wsname="wcombo", const char* dataname="dOS_DATA", const char* ACTag="", double CI = 0.95, int calculatorType = 2, int testStatType = 2,  bool useCLs = false, double poival=-1, int rndseed=-1);

#include "StandardHypoTestInvDemo.C"


pair<double,double> runLimit_RaaNS_Workspace(const char *filename, const char *poiname, const char *pdfname, const char *wsname, const char* dataname, const char* ACTag, double CI, int calculatorType, int testStatType,  bool useCLs, double poival, int rndseed)
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

   double poimin,poimax;
   if (calculatorType==0 && rndseed>0) {
      poimin=poival;
      poimax=poival;
   } else {
      // get a first estimate of where to look

      // pdf->fitTo(*data, Offset(kTRUE), Extended(kTRUE), NumCPU(2), Range("MassWindow"), Save(), Minimizer("Minuit2","Migrad"));
      // pdf->fitTo(*data, Extended(kTRUE), Offset(kTRUE));
      ProfileLikelihoodCalculator pl(*data,*sbHypo,1.-CI);
      pl.SetConfidenceLevel(CI); 
      int ci = 100*CI;
      LikelihoodInterval* interval = pl.GetInterval();
      LikelihoodIntervalPlot plot(interval);
      TCanvas c4; c4.cd(); 
      plot.SetRange(0.,1.,0.05,3.);
      // plot.Draw();
      TLatex latexCI;
      latexCI.SetNDC();
      latexCI.SetTextSize(0.035);
      latexCI.DrawLatex(0.5,1.-0.05*2,Form("%s %d %s C.I.",poiname,ci,"%"));
      latexCI.DrawLatex(0.5,1.-0.05*3,Form("Upper limit: %f",interval->UpperLimit(*theVar)));
      latexCI.DrawLatex(0.5,1.-0.05*4,Form("Lower limit: %f",interval->LowerLimit(*theVar)));
      TString intrvlName = theVar->GetTitle();
      // print out the iterval on the Parameter of Interest
      cout <<endl<< CI <<"\% interval on " <<theVar->GetName()<<" is : ["<<
         interval->LowerLimit(*theVar) << ", "<<
         interval->UpperLimit(*theVar) << "] "<<endl;
      pair<double, double> CnfdncIntrvl;
      CnfdncIntrvl.first  = interval->LowerLimit(*theVar);
      CnfdncIntrvl.second = interval->UpperLimit(*theVar);
      c4.SaveAs("ULtest.pdf");

      // pair<double, double> CnfdncIntrvl;
      // CnfdncIntrvl.first  = theVar->getVal()+theVar->getErrorLo();
      // CnfdncIntrvl.second = theVar->getVal()+theVar->getErrorHi();


      of << 100.*CI << "% CL confidence interval from a likelihood scan on " << poiname << endl;
      of << CnfdncIntrvl.first << " " << CnfdncIntrvl.second << endl;

      // now let's move to the hypo test inverter
      /////////////////////////////////////////////////////////////
      // Now get the POI for convenience
      // you may want to adjust the range of your POI
      ////////////////////////////////////////////////////////////
      RooRealVar* firstPOI = (RooRealVar*) sbHypo->GetParametersOfInterest()->first();
      poimin = 0.5*max((double) 0.,CnfdncIntrvl.first);
      poimax = 1.5*CnfdncIntrvl.second;
      cout << "Will scan between " << poimin << " and " << poimax << endl;
      // firstPOI->setMin(poimin);
      // firstPOI->setMax(poimax);
      cout << "POI defined between" << firstPOI->getMin() << " and " << firstPOI->getMax() << endl;
   }

   // actually run the hypo test inverter

   // define calc type and test stat type
   // type = 0 Freq calculator 
   // type = 1 Hybrid calculator
   // type = 2 Asymptotic calculator  
   // type = 3 Asymptotic calculator using nominal Asimov data sets (not using fitted parameter values but nominal ones)
   //
   // testStatType = 0 LEP
   //              = 1 Tevatron 
   //              = 2 Profile Likelihood two sided
   //              = 3 Profile Likelihood one sided (i.e. = 0 if mu < mu_hat)
   //              = 4 Profile Likelihood signed ( pll = -pll if mu < mu_hat) 
   //              = 5 Max Likelihood Estimate as test statistic
   //              = 6 Number of observed event as test statistic
   //

   // 0,3,true for frequentist CLs; 2,3,true for asymptotic CLs; 0,2,false for FC
//   int calculatorType = 2;
//   int testStatType = 2;
//   bool useCLs = false;
   // root> StandardHypoTestInvDemo("fileName","workspace name","S+B modelconfig name","B model name","data set name",calculator type, test statistic type, use CLS,
   //                                number of points, xmin, xmax, number of toys, use number counting)


   TStopwatch tw; 
   tw.Start();

   int npoints=20;
   int ntoys=1000;

   if (calculatorType==0 && rndseed>0) {
      npoints=1;
      ntoys=20;
   }

   HypoTestInverterResult* r = StandardHypoTestInvDemo(f->GetName(),
         wsname,
         sbHypo,
         bHypo,
         dataname,
         calculatorType,
         testStatType, 
         useCLs,
         npoints,
         poimin,
         poimax,
         ntoys,
         false,
         0,
         CI,
         ACTag,
         rndseed);

   pair<double,double> lims;
   lims.first = r->LowerLimit();
   lims.second = r->UpperLimit();

   std::streambuf *coutbuf = std::cout.rdbuf();
   std::cout.rdbuf(of.rdbuf());
   tw.Print();
   std::cout.rdbuf(coutbuf);

   of << "Feldman-Cousins interval from " 
   << "StandardHypoTestInvDemo(" << f->GetName() << ","
         << wsname << ","
         <<sbHypo<<","
         <<bHypo<<","
         <<dataname<<","
         <<calculatorType<<","
         <<testStatType<<","
         <<useCLs<<","  
         <<20<<","   
         <<poimin<<","  
         <<poimax<<"," 
         <<1000<<","
         <<false<<","
         <<0<<");" << endl;
   of << lims.first << " " << lims.second << endl;
  
  return lims;
}
