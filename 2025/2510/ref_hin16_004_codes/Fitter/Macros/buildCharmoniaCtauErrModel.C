#ifndef buildCharmoniaCtauErrModel_C
#define buildCharmoniaCtauErrModel_C

#include "Utilities/initClasses.h"
#include "RooStats/SPlot.h"

using namespace RooStats;

void setCtauErrDefaultParameters(map<string, string> &parIni, bool isPbPb, double numEntries);
bool histToPdf(RooWorkspace& ws, TH1D* hist, string pdfName, vector<double> rangeErr);
bool makeCtauErrPdf(RooWorkspace& ws, vector<TH1D*>& outputHist, vector<string> objectColl, vector<string> rangeColl, vector<string> cutColl, string dsName, struct KinCuts cut, double binWidth);
bool createCtauErrTemplateUsingMatrix(RooWorkspace& ws, string dsName, string pdfType, struct KinCuts cut,  map<string,string> parIni, bool incJpsi, bool incPsi2S, double binWidth);
bool createCtauErrTemplateUsingSPLOT(RooWorkspace& ws, string dsName, string pdfType, struct KinCuts cut, bool incJpsi, bool incPsi2S, double binWidth);
bool addCtauErrModel(RooWorkspace& ws, string object, string pdfType, map<string,string> parIni, bool isPbPb);
TH1* rebinhist(TH1 *hist, double xmin=1e99, double xmax=-1e99);
void MatrixInverse2x2(TMatrixD& Inverse);


bool buildCharmoniaCtauErrModel(RooWorkspace& ws, map<string, string>  parIni, 
                                struct KinCuts cut,          // Variable containing all kinematic cuts
                                string dsName,               // Name of current input dataset
                                bool incJpsi,                // Include Jpsi model
                                bool incPsi2S,               // Include Psi(2S) model
                                double  binWidth,            // Bin width
                                double  numEntries = 300000. // Number of entries in the dataset
                                )
{


  bool isMC = false;
  if (dsName.find("MC")!=std::string::npos) { isMC = true; }
  bool isPbPb = false;
  if (dsName.find("PbPb")!=std::string::npos) { isPbPb = true; }

  // If the initial parameters are empty, set defaul parameter values
  setCtauErrDefaultParameters(parIni, isPbPb, numEntries);

  // C r e a t e   m o d e l 

  string pdfType = "pdfCTAUERR";
  if(!createCtauErrTemplateUsingSPLOT(ws, dsName, pdfType, cut, incJpsi, incPsi2S, binWidth)) { cout << "[ERROR] Creating the Ctau Error Templates using sPLOT failed" << endl; return false; }
  //if(!createCtauErrTemplateUsingMatrix(ws, dsName, pdfType, cut, parIni, incJpsi, incPsi2S, binWidth)) { cout << "[ERROR] Creating the Ctau Error Templates failed" << endl; return false; }
  if (incJpsi)  {  if(!addCtauErrModel(ws, "Jpsi", pdfType, parIni, isPbPb))  { return false; }  }
  if (incPsi2S) {  if(!addCtauErrModel(ws, "Psi2S", pdfType, parIni, isPbPb)) { return false; }  }
  if (!isMC)    {  if(!addCtauErrModel(ws, "Bkg", pdfType, parIni, isPbPb))   { return false; }  }

  // Total PDF
  
  RooArgList pdfList;
  if (incJpsi)  { pdfList.add( *ws.pdf(Form("%sTot_Jpsi_%s", pdfType.c_str(), (isPbPb?"PbPb":"PP"))) );  }
  if (incPsi2S) { pdfList.add( *ws.pdf(Form("%sTot_Psi2S_%s", pdfType.c_str(), (isPbPb?"PbPb":"PP"))) ); }
  if (!isMC)    { pdfList.add( *ws.pdf(Form("%sTot_Bkg_%s", pdfType.c_str(), (isPbPb?"PbPb":"PP"))) );   }
  if (!incJpsi && !incPsi2S && isMC) {
    cout << "[ERROR] User did not include any model for MC, please fix your input settings!" << endl; return false;
  }
  string pdfName = Form("%s_Tot_%s", pdfType.c_str(), (isPbPb?"PbPb":"PP"));
  RooAbsPdf *themodel = new RooAddPdf(pdfName.c_str(), pdfName.c_str(), pdfList);
  ws.import(*themodel);
  ws.pdf(pdfName.c_str())->setNormRange("CtauErrWindow");
  
  // save the initial values of the model we've just created
  RooArgSet* params = (RooArgSet*) themodel->getParameters(RooArgSet(*ws.var("ctauErr")));
  ws.saveSnapshot((pdfName+"_parIni").c_str(),*params,kTRUE);

  return true;
};


bool addCtauErrModel(RooWorkspace& ws, string object, string pdfType, map<string,string> parIni, bool isPbPb)
{
  if (ws.pdf(Form("%sTot_%s_%s", pdfType.c_str(), object.c_str(), (isPbPb?"PbPb":"PP")))) { 
    cout << Form("[ERROR] The %s Extended Ctau Error Model has already been implemented!", object.c_str()) << endl;
    return false; 
  }

  if (!ws.var(Form("N_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")))){ ws.factory( parIni[Form("N_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str() ); } 
  // create the Extended PDF
  ws.factory(Form("RooExtendPdf::%s(%s,%s)", Form("%sTot_%s_%s", pdfType.c_str(), object.c_str(), (isPbPb?"PbPb":"PP")),
                  Form("%s_%s_%s", pdfType.c_str(), object.c_str(), (isPbPb?"PbPb":"PP")),
                  Form("N_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))
                  ));

  return true;
};


bool createCtauErrTemplateUsingSPLOT(RooWorkspace& ws, string dsName, string pdfType, struct KinCuts cut, bool incJpsi, bool incPsi2S, double binWidth)
{
  string hType = pdfType;
  hType.replace(hType.find("pdf"), string("pdf").length(), "h");

  bool isPbPb = false;
  if (dsName.find("PbPb")!=std::string::npos) { isPbPb = true; }
  if (dsName.find("MC")!=std::string::npos)   { cout << "[ERROR] sPlot can only be run on data!" << endl; return false;  }  // Only accept data

  string pdfMassName = Form("pdfMASS_Tot_%s", (isPbPb?"PbPb":"PP"));
  RooArgList yieldList;
  if (incJpsi)  { yieldList.add( *ws.var(Form("N_Jpsi_%s", (isPbPb?"PbPb":"PP"))) );  }
  if (incPsi2S) { yieldList.add( *ws.var(Form("N_Psi2S_%s", (isPbPb?"PbPb":"PP"))) ); }
  yieldList.add( *ws.var(Form("N_Bkg_%s", (isPbPb?"PbPb":"PP"))) ); // Always add background
  RooDataSet* data = (RooDataSet*)ws.data(dsName.c_str())->Clone("TMP_DATA");
  RooAbsPdf* pdf = clone(*ws.pdf(pdfMassName.c_str()));

  RooStats::SPlot sData = RooStats::SPlot("sData","An SPlot", *data, pdf, yieldList);
  ws.import(*data, Rename((dsName+"_SPLOT").c_str()));
  if (incJpsi) {
    cout <<  "[INFO] Jpsi yield -> Mass Fit: " << ws.var(Form("N_Jpsi_%s", (isPbPb?"PbPb":"PP")))->getVal() << " , sWeights: " << sData.GetYieldFromSWeight(Form("N_Jpsi_%s", (isPbPb?"PbPb":"PP"))) << std::endl;
  }
  if (incPsi2S) {
    cout <<  "[INFO] Psi2S yield -> Mass Fit: " << ws.var(Form("N_Psi2S_%s", (isPbPb?"PbPb":"PP")))->getVal() << " , sWeights: " << sData.GetYieldFromSWeight(Form("N_Psi2S_%s", (isPbPb?"PbPb":"PP"))) << std::endl;
  }
  cout <<  "[INFO] Bkg yield -> Mass Fit: " << ws.var(Form("N_Bkg_%s", (isPbPb?"PbPb":"PP")))->getVal() << " , sWeights: " << sData.GetYieldFromSWeight(Form("N_Bkg_%s", (isPbPb?"PbPb":"PP"))) << std::endl;
  
  // create weighted data sets
  double ctauErrMax = cut.dMuon.ctauErr.Max, ctauErrMin = cut.dMuon.ctauErr.Min;
  vector<double> rangeErr; rangeErr.push_back(cut.dMuon.ctauErr.Min); rangeErr.push_back(cut.dMuon.ctauErr.Max);
  int nBins = min(int( round((ctauErrMax - ctauErrMin)/binWidth) ), 1000);

  TH1D* hTot = (TH1D*)ws.data(dsName.c_str())->createHistogram(Form("%s_Tot_%s", hType.c_str(), (isPbPb?"PbPb":"PP")), *ws.var("ctauErr"), Binning(nBins, ctauErrMin, ctauErrMax));
  if ( !histToPdf(ws, hTot, Form("%sTot_Tot_%s", pdfType.c_str(), (isPbPb?"PbPb":"PP")), rangeErr)) { return false; }
  hTot->Delete();

  RooDataSet* dataw_Bkg  = new RooDataSet("TMP_BKG_DATA","TMP_BKG_DATA", (RooDataSet*)ws.data((dsName+"_SPLOT").c_str()), RooArgSet(*ws.var("ctauErr"), *ws.var(Form("N_Bkg_%s_sw", (isPbPb?"PbPb":"PP")))), 0, Form("N_Bkg_%s_sw", (isPbPb?"PbPb":"PP")));
  TH1D* hBkg = (TH1D*)dataw_Bkg->createHistogram(Form("%s_Bkg_%s_sWEIGHT", hType.c_str(), (isPbPb?"PbPb":"PP")), *ws.var("ctauErr"), Binning(nBins, ctauErrMin, ctauErrMax));
  if ( !histToPdf(ws, hBkg, Form("%s_Bkg_%s", pdfType.c_str(), (isPbPb?"PbPb":"PP")), rangeErr)) { return false; }
  delete hBkg; delete dataw_Bkg;
  
  if (incJpsi) {
    RooDataSet* dataw_Jpsi  = new RooDataSet("TMP_JPSI_DATA","TMP_JPSI_DATA", (RooDataSet*)ws.data((dsName+"_SPLOT").c_str()), RooArgSet(*ws.var("ctauErr"), *ws.var(Form("N_Jpsi_%s_sw", (isPbPb?"PbPb":"PP")))), 0, Form("N_Jpsi_%s_sw", (isPbPb?"PbPb":"PP")));
    TH1D* hJpsi = (TH1D*)dataw_Jpsi->createHistogram(Form("%s_Jpsi_%s_sWEIGHT", hType.c_str(), (isPbPb?"PbPb":"PP")), *ws.var("ctauErr"), Binning(nBins, ctauErrMin, ctauErrMax));
    if ( !histToPdf(ws, hJpsi, Form("%s_Jpsi_%s", pdfType.c_str(), (isPbPb?"PbPb":"PP")), rangeErr)) { return false; }
    delete hJpsi; delete dataw_Jpsi;
  }
  if (incPsi2S) {
    RooDataSet* dataw_Psi2S  = new RooDataSet("TMP_PSI2S_DATA","TMP_PSI2S_DATA", (RooDataSet*)ws.data((dsName+"_SPLOT").c_str()), RooArgSet(*ws.var("ctauErr"), *ws.var(Form("N_Psi2S_%s_sw", (isPbPb?"PbPb":"PP")))), 0, Form("N_Psi2S_%s_sw", (isPbPb?"PbPb":"PP")));
    TH1D* hPsi2S = (TH1D*)dataw_Psi2S->createHistogram(Form("%s_Psi2S_%s_sWEIGHT", hType.c_str(), (isPbPb?"PbPb":"PP")), *ws.var("ctauErr"), Binning(nBins, ctauErrMin, ctauErrMax));
    if ( !histToPdf(ws, hPsi2S, Form("%s_Psi2S_%s", pdfType.c_str(), (isPbPb?"PbPb":"PP")), rangeErr)) { return false; }
    delete hPsi2S; delete dataw_Psi2S;
  }

  return true;
};
  

bool createCtauErrTemplateUsingMatrix(RooWorkspace& ws, string dsName, string pdfType, struct KinCuts cut,  map<string,string> parIni, bool incJpsi, bool incPsi2S, double binWidth)
{
  string hType = pdfType;
  hType.replace(hType.find("pdf"), string("pdf").length(), "h");

  bool isPbPb = false;
  if (dsName.find("PbPb")!=std::string::npos) { isPbPb = true; }
  bool isMC = false;
  if (dsName.find("MC")!=std::string::npos) { isMC = true; }

  double ctauErrMax = cut.dMuon.ctauErr.Max, ctauErrMin = cut.dMuon.ctauErr.Min;
  int nBins = min(int( round((ctauErrMax - ctauErrMin)/binWidth) ), 1000);
  TH1D* hTot = (TH1D*)ws.data(dsName.c_str())->createHistogram(Form("%s_Tot_%s", hType.c_str(), (isPbPb?"PbPb":"PP")), *ws.var("ctauErr"), Binning(nBins, ctauErrMin, ctauErrMax));
  vector<double> rangeErr; rangeErr.push_back(cut.dMuon.ctauErr.Min); rangeErr.push_back(cut.dMuon.ctauErr.Max);

  // Set the range of the ctau Error
  if (isMC) {
    TH1D* hSig = (TH1D*)hTot->Clone("SIGTMP");
    if (incJpsi != incPsi2S) {
      if ( !histToPdf(ws, hSig, Form("%s_%s_%s", pdfType.c_str(), (incJpsi?"Jpsi":"Psi2S"), (isPbPb?"PbPb":"PP")), rangeErr)) { return false; }
    }
    delete hSig;
  }
  else {
    TH1D* hBkg = (TH1D*)hTot->Clone("BKGTMP");
    if (incJpsi) {
      vector<string> rangeColl; vector<string> cutColl; vector<string> objectColl;
      rangeColl.push_back("JpsiWindow"); cutColl.push_back(parIni["JpsiMassRange_Cut"]);  objectColl.push_back("Jpsi");
      rangeColl.push_back(parIni["BkgMassRange_JPSI_Label"]); cutColl.push_back(parIni["BkgMassRange_JPSI_Cut"]);   objectColl.push_back("Bkg");
      vector<TH1D*> outputHist_JPSI;
      if ( !makeCtauErrPdf(ws, outputHist_JPSI, objectColl, rangeColl, cutColl, dsName, cut, binWidth) ) { return false; }
      hBkg->Add(outputHist_JPSI.at(0), -1.0);
      if ( !histToPdf(ws, outputHist_JPSI.at(0), Form("%s_Jpsi_%s", pdfType.c_str(), (isPbPb?"PbPb":"PP")), rangeErr)) { return false; }
      if (!incPsi2S) {
        if ( !histToPdf(ws, hBkg/*outputHist_JPSI.at(1)*/, Form("%s_Bkg_%s", pdfType.c_str(), (isPbPb?"PbPb":"PP")), rangeErr)) { return false; }
      }
      delete outputHist_JPSI.at(0); delete outputHist_JPSI.at(1);
    }
    if (incPsi2S) {
      vector<string> rangeColl; vector<string> cutColl; vector<string> objectColl;
      rangeColl.push_back("Psi2SWindow"); cutColl.push_back(parIni["Psi2SMassRange_Cut"]);  objectColl.push_back("Psi2S");
      rangeColl.push_back(parIni["BkgMassRange_PSI2S_Label"]); cutColl.push_back(parIni["BkgMassRange_PSI2S_Cut"]);   objectColl.push_back("Bkg");
      vector<TH1D*> outputHist_PSI2S;
      if ( !makeCtauErrPdf(ws, outputHist_PSI2S, objectColl, rangeColl, cutColl, dsName, cut, binWidth) ) { return false; }
      hBkg->Add(outputHist_PSI2S.at(0), -1.0);
      if ( !histToPdf(ws, outputHist_PSI2S.at(0), Form("%s_Psi2S_%s", pdfType.c_str(), (isPbPb?"PbPb":"PP")), rangeErr)) { return false; }
      if (!incJpsi) {
        if ( !histToPdf(ws, outputHist_PSI2S.at(1), Form("%s_Bkg_%s", pdfType.c_str(), (isPbPb?"PbPb":"PP")), rangeErr)) { return false; }
      }
      delete outputHist_PSI2S.at(0); delete outputHist_PSI2S.at(1);
    }
    if (incPsi2S == incJpsi) {
      if ( !histToPdf(ws, hBkg, Form("%s_Bkg_%s", pdfType.c_str(), (isPbPb?"PbPb":"PP")), rangeErr)) { return false; }
    }
    delete hBkg;
  }
  if ( !histToPdf(ws, hTot, Form("%sTot_Tot_%s", pdfType.c_str(), (isPbPb?"PbPb":"PP")), rangeErr)) { return false; }
  hTot->Delete();

  return true;
};


bool makeCtauErrPdf(RooWorkspace& ws, vector<TH1D*>& outputHist, vector<string> objectColl, vector<string> rangeColl, vector<string> cutColl, string dsName, struct KinCuts cut, double binWidth)
{
  if (objectColl.size()!=rangeColl.size() || objectColl.size()<1) {
    cout << "[ERROR] Object and range collections are wrong, CtauErrPdf maker failed!" << endl; return false;
  }

  double ctauErrMax = cut.dMuon.ctauErr.Max, ctauErrMin = cut.dMuon.ctauErr.Min;
  int nBins = min(int( round((ctauErrMax - ctauErrMin)/binWidth) ), 1000);

  bool isPbPb = false;
  if (dsName.find("PbPb")!=std::string::npos) { isPbPb = true; }

  // Extract the PDFs and the Number of events for each physics object
  vector<Double_t>    N;
  vector<RooAbsPdf*>  pdfMass;
  for (unsigned int i=0; i<objectColl.size(); i++) {
    if (!ws.pdf(Form("pdfMASS_%s_%s", objectColl.at(i).c_str(), (isPbPb?"PbPb":"PP")))) {
      cout << Form("[ERROR] Mass PDF for %s is missing!", objectColl.at(i).c_str()) << endl; return false;
    } else {
      pdfMass.push_back(ws.pdf(Form("pdfMASS_%s_%s", objectColl.at(i).c_str(), (isPbPb?"PbPb":"PP"))));
    }
    if (ws.var(Form("N_%s_%s", objectColl.at(i).c_str(), (isPbPb?"PbPb":"PP")))) { 
      N.push_back(ws.var(Form("N_%s_%s", objectColl.at(i).c_str(), (isPbPb?"PbPb":"PP")))->getValV());
    } else if (ws.function(Form("N_%s_%s", objectColl.at(i).c_str(), (isPbPb?"PbPb":"PP")))) {
      N.push_back(ws.function(Form("N_%s_%s", objectColl.at(i).c_str(), (isPbPb?"PbPb":"PP")))->getValV());
    } else {
     cout << Form("[ERROR] N parameter for %s is missing!", objectColl.at(i).c_str()) << endl; return false;
    } 
  }

  // Extract the ctau error histograms for each mass range 
  vector<TH1*> histCtauErr;
  for (unsigned int j=0; j<rangeColl.size(); j++) {
    RooDataSet* data = (RooDataSet*)ws.data(dsName.c_str())->reduce(*ws.var("ctauErr"), cutColl.at(j).c_str());
    histCtauErr.push_back( (TH1*)data->createHistogram(Form("hCtauErr_Range_%s_%s", rangeColl.at(j).c_str(), (isPbPb?"PbPb":"PP")), *ws.var("ctauErr"), Binning(nBins, ctauErrMin, ctauErrMax)) );
    histCtauErr.at(j)->Scale(1/histCtauErr.at(j)->GetEntries());
    delete data;
  }

  // Get the mass yield coefficients by integrating the mass PDFs
  TMatrixD CoefficientMatrix = TMatrixD(objectColl.size(), rangeColl.size());
  for (unsigned int j=0; j<rangeColl.size(); j++) {
    Double_t Norm = 0;
    for (unsigned int i=0; i<objectColl.size(); i++) {
      CoefficientMatrix(j, i) = pdfMass.at(i)->createIntegral(*ws.var("invMass"), NormSet(*ws.var("invMass")), Range(rangeColl.at(j).c_str()))->getValV();
      CoefficientMatrix(j, i) = N.at(i) * CoefficientMatrix(j, i);
      Norm = Norm + CoefficientMatrix(j, i);
    }
    if (Norm>0.0) {
      double testNorm = 0.0;
      for (unsigned int i=0; i<objectColl.size(); i++) {
        CoefficientMatrix(j, i) = CoefficientMatrix(j, i) / Norm;
        testNorm += CoefficientMatrix(j, i);
      }
      if (abs(testNorm-1.0)>0.000000000001) { cout << "[ERROR] Sum of alphas is: " << Form("%.10f",testNorm) << " !" << endl; return false; }
    } else {
      cout << "[ERROR] Normalization factor for Range: " << rangeColl.at(j) << " is invalid: " << Form("%.4f",Norm) << endl; return false;
    }
  }

  // Get the inverse of the coefficients
  CoefficientMatrix.Print();
  Double_t Det = CoefficientMatrix.Determinant();
  TMatrixD InvMatrix = CoefficientMatrix;
  if (Det>0) {
    InvMatrix.Invert();
    TMatrixD TestMatrix = InvMatrix * CoefficientMatrix;
    double testDet = TestMatrix.Determinant();
    if (fabs(testDet-1.0)>0.00000000001) { cout << Form("[ERROR] Determinant of invMatrix * CoefficienctMatrix is: %.10f", testDet) << endl; return false; }
  } else {
    cout << "[ERROR] Determinant of Coefficient Matrix is invalid: " << Form("%.4f",Det) << endl; return false;
  }

  // Get the ctau error histrograms for each physics object
  for (unsigned int j=0; j<objectColl.size(); j++) {
    outputHist.push_back(new TH1D(Form("h_%s_%s", objectColl.at(j).c_str(), (isPbPb?"PbPb":"PP")), "", nBins, ctauErrMin, ctauErrMax));
    for (unsigned int i=0; i<rangeColl.size(); i++) {
      outputHist.at(j)->Add(histCtauErr.at(i), InvMatrix(j, i));
    }
    outputHist.at(j)->Scale(N.at(j)*pdfMass.at(j)->createIntegral(*ws.var("invMass"), NormSet(*ws.var("invMass")), Range("FullWindow"))->getValV());
  }

  for (unsigned int i=0; i<rangeColl.size(); i++) {
    histCtauErr.at(i)->Delete();
  }
   
  return true;
  
};


bool histToPdf(RooWorkspace& ws, TH1D* hist, string pdfName, vector<double> rangeErr)
{
  if (ws.pdf(pdfName.c_str())) { 
    cout << Form("[INFO] The %s Template has already been created!", pdfName.c_str()) << endl;
    return true; 
  }
  // Cleaning the input histogram
  // 1) Remove the Under and Overflow bins
  hist->ClearUnderflowAndOverflow();
  // 2) Set negative bin content to zero
  for (int i=0; i<=hist->GetNbinsX(); i++) { if (hist->GetBinContent(i)<0) { hist->SetBinContent(i, 0.0000000001); } }
  // 2) Reduce the range of histogram and rebin it
  TH1* hClean = rebinhist(hist, rangeErr[0], rangeErr[1]);

  cout << Form("[INFO] Implementing %s Template", pdfName.c_str()) << endl;

  string dataName = pdfName;
  dataName.replace(dataName.find("pdf"), string("pdf").length(), "dh");
  RooDataHist* dataHist = new RooDataHist(dataName.c_str(), "", *ws.var("ctauErr"), hClean);
  if (dataHist==NULL) { cout << "[ERROR] DataHist used to create " << pdfName << " failed!" << endl; return false; } 
  if (dataHist->sumEntries()==0) { cout << "[ERROR] DataHist used to create " << pdfName << " is empty!" << endl; return false; } 
  if (fabs(dataHist->sumEntries()-hClean->GetSumOfWeights())>0.001) { cout << "[ERROR] DataHist used to create " << pdfName << "  " << " is invalid!  " << endl; return false; } 
  ws.import(*dataHist);

  bool isPbPb=false;
  if (pdfName.find("PbPb")!=std::string::npos) isPbPb=true;
  RooHistPdf* pdf = new RooHistPdf(pdfName.c_str(), pdfName.c_str(), *ws.var("ctauErr"), *((RooDataHist*)ws.data(dataName.c_str())));
  //RooKeysPdf* pdf = new RooKeysPdf(pdfName.c_str(), pdfName.c_str(), *ws.var("ctauErr"), *((RooDataSet*)ws.data(dataName.c_str())),RooKeysPdf::NoMirror, isPbPb?0.4:0.4);
  if (pdf==NULL) { cout << "[ERROR] RooKeysPDF " << pdfName << " is NULL!" << endl; return false; } 
  ws.import(*pdf);

  delete pdf;
  delete dataHist;
 
  return true;
};


TH1* rebinhist(TH1 *hist, double xmin, double xmax)
{
  TH1 *hcopy = (TH1*) hist->Clone("hcopy");

  // range of the new hist
  int imin = hcopy->FindBin(xmin);
  if (imin>=hcopy->GetNbinsX()) imin=1;
  int imax = hcopy->FindBin(0.999999*xmax);
  if (imax<=1) imax=hcopy->GetNbinsX();

  vector<double> newbins;
  newbins.push_back(hcopy->GetBinLowEdge(imin));
  for (int i=imin; i<=imax; i++) {
    if (hcopy->GetBinContent(i)>0.1) {
      newbins.push_back(hcopy->GetBinLowEdge(i)+hcopy->GetBinWidth(i));
    } else {
      int nrebin=2;
      for (i++; i<=imax; i++) {
        if (hcopy->GetBinContent(i)>0.00000001) {
          newbins.push_back(hcopy->GetBinLowEdge(i)+hcopy->GetBinWidth(i));
          hcopy->SetBinContent(i,hcopy->GetBinContent(i)/nrebin);
          break;
        }
        nrebin++;
      }
    }
  }

  if (xmin < newbins[1]) newbins[0] = xmin;
  if (xmax > newbins[newbins.size()-2]) newbins[newbins.size()-1] = xmax;

  TH1 *ans = hcopy->Rebin(newbins.size()-1,"hnew",newbins.data());

  delete hcopy;
  return ans;
};


void setCtauErrDefaultParameters(map<string, string> &parIni, bool isPbPb, double numEntries)
{

  cout << "[INFO] Setting user undefined initial parameters to their default values" << endl;

  // DEFAULT RANGE OF NUMBER OF EVENTS
  if (parIni.count(Form("N_Jpsi_%s", (isPbPb?"PbPb":"PP")))==0 || parIni[Form("N_Jpsi_%s", (isPbPb?"PbPb":"PP"))]=="") { 
    parIni[Form("N_Jpsi_%s", (isPbPb?"PbPb":"PP"))]  = Form("%s[%.12f,%.12f,%.12f]", Form("N_Jpsi_%s", (isPbPb?"PbPb":"PP")), numEntries, 0.0, numEntries*2.0);
  }
  if (parIni.count(Form("N_Psi2S_%s", (isPbPb?"PbPb":"PP")))==0 || parIni[Form("N_Psi2S_%s", (isPbPb?"PbPb":"PP"))]=="") { 
    parIni[Form("N_Psi2S_%s", (isPbPb?"PbPb":"PP"))]  = Form("%s[%.12f,%.12f,%.12f]", Form("N_Psi2S_%s", (isPbPb?"PbPb":"PP")), numEntries, 0.0, numEntries*2.0);
  }
  if (parIni.count(Form("N_Bkg_%s", (isPbPb?"PbPb":"PP")))==0 || parIni[Form("N_Bkg_%s", (isPbPb?"PbPb":"PP"))]=="") { 
    parIni[Form("N_Bkg_%s", (isPbPb?"PbPb":"PP"))]  = Form("%s[%.12f,%.12f,%.12f]", Form("N_Bkg_%s", (isPbPb?"PbPb":"PP")), numEntries, 0.0, numEntries*2.0);
  }

  return;
};


void MatrixInverse2x2(TMatrixD& Inverse)
{
  Double_t A = Inverse(0,0)  ,  B = Inverse(1,0);
  Double_t C = Inverse(0,1)  ,  D = Inverse(1,1);

  Double_t det = ( ( A*D ) - ( B*C ) );

  Double_t Ainv = ( D/det )       ,  Binv = -1.0*( B/det );
  Double_t Cinv = -1.0*( C/det )  ,  Dinv = ( A/det );

  Inverse(0,0) = Ainv;   Inverse(1,0) = Binv;
  Inverse(0,1) = Cinv;   Inverse(1,1) = Dinv;

  return;
};


#endif // #ifndef buildCharmoniaCtauErrModel_C
