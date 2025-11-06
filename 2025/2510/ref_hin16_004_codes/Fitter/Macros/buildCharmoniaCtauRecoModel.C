#ifndef buildCharmoniaCtauRecoModel_C
#define buildCharmoniaCtauRecoModel_C

#include "Utilities/initClasses.h"

bool createCtauRecoTemplate(RooWorkspace& ws, string dsName, string pdfType, struct KinCuts cut, bool incJpsi, bool incPsi2S, double binWidth);
void setCtauRecoDefaultParameters(map<string, string> &parIni, bool isPbPb, double numEntries);
bool ctauHistToPdf(RooWorkspace& ws, TH1D* hist, string pdfName, string dsName, vector<double> range, bool useDataSet=true);
TH1* rebinctauhist(TH1 *hist, double xmin=1e99, double xmax=-1e99);


bool buildCharmoniaCtauRecoModel(RooWorkspace& ws, map<string, string>  parIni,
                                 struct KinCuts cut,          // Variable containing all kinematic cuts
                                 string dsName,               // Name of current input dataset
                                 bool incJpsi,                // Include Jpsi model
                                 bool incPsi2S,                // Include Psi(2S) model
                                 double  binWidth,            // Bin width
                                 double  numEntries = 300000. // Number of entries in the dataset
                                 )
{
  bool isPbPb = false;
  if (dsName.find("PbPb")!=std::string::npos) { isPbPb = true; }

  // If the initial parameters are empty, set defaul parameter values
  setCtauRecoDefaultParameters(parIni, isPbPb, numEntries);

  // C r e a t e   m o d e l
  // Total PDF
  string pdfType = "pdfCTAU";
  string pdfName = Form("%s_%s_%s", pdfType.c_str(), incJpsi?"JpsiNoPR":"Psi2SNoPR", (isPbPb?"PbPb":"PP"));
  
  if(!createCtauRecoTemplate(ws, dsName, pdfType, cut, incJpsi, incPsi2S, binWidth)) { cout << "[ERROR] Creating the Ctau Reco Template failed" << endl; return false; }

  if (incJpsi && incPsi2S) {
    cout << "[ERROR] User included Jpsi and Psi2S together, but can only template them separately, please fix your input settings!" << endl; return false;
  }
  if (!incJpsi && !incPsi2S) {
    cout << "[ERROR] User did not include any model, please fix your input settings!" << endl; return false;
  }
  
//  setFixedVarsToContantVars(ws);

  // save the initial values of the model we've just created
  RooArgSet* params = (RooArgSet*) ((ws.pdf(pdfName.c_str()))->getParameters(RooArgSet(*ws.var("ctau"))));
  
  ws.saveSnapshot((pdfName+"_parIni").c_str(),*params,kTRUE);
  
  return true;
};


bool createCtauRecoTemplate(RooWorkspace& ws, string dsName, string pdfType, struct KinCuts cut, bool incJpsi, bool incPsi2S, double binWidth)
{
  string hType = pdfType;
  hType.replace(hType.find("pdf"), string("pdf").length(), "h");
  
  bool isPbPb = false;
  if (dsName.find("PbPb")!=std::string::npos) { isPbPb = true; }
  if (dsName.find("MC")==std::string::npos)   { std::cout << "[ERROR] Ctau Reco templates can not be build using data!" << std::endl; return false;  }  // Only accept MC
  
  // create weighted data sets
  double ctauRecoMax = cut.dMuon.ctauTrue.Max, ctauRecoMin = cut.dMuon.ctauTrue.Min;
  vector<double> rangeErr; rangeErr.push_back(cut.dMuon.ctauTrue.Min); rangeErr.push_back(cut.dMuon.ctauTrue.Max);
  int nBins = min(int( round((ctauRecoMax - ctauRecoMin)/binWidth) ), 1000);

  TH1D* hTot = (TH1D*)ws.data(dsName.c_str())->createHistogram(Form("%s_%s_%s", hType.c_str(), incJpsi?"JpsiNoPR":"Psi2SNoPR", (isPbPb?"PbPb":"PP")), *ws.var("ctau"), Binning(nBins, ctauRecoMin, ctauRecoMax));
  RooDataHist* dataHist = new RooDataHist("tmp", "", *ws.var("ctau"), hTot); delete dataHist; // KEEP THIS LINE, IT IS A FIX FOR ROOKEYS
  if ( !ctauHistToPdf(ws, hTot, Form("%s_%s_%s", pdfType.c_str(), incJpsi?"JpsiNoPR":"Psi2SNoPR", (isPbPb?"PbPb":"PP")), dsName, rangeErr)) { return false; }
  hTot->Delete();
  
  return true;
};


bool ctauHistToPdf(RooWorkspace& ws, TH1D* hist, string pdfName, string dsName, vector<double> range, bool useDataSet)
{
  if (ws.pdf(pdfName.c_str())) {
    cout << Form("[INFO] The %s Template has already been created!", pdfName.c_str()) << endl;
    return true;
  }
  
  cout << Form("[INFO] Implementing %s Template", pdfName.c_str()) << endl;
  
  if (useDataSet) {
    cout << Form("[INFO] Implementing %s Template using RooKeysPdf", pdfName.c_str()) << endl;
    if (ws.data(dsName.c_str())==NULL) { cout << "[ERROR] DataSet " << dsName << " for Bkg Template was not found!" << endl; return false; }
    RooKeysPdf* pdf = new RooKeysPdf(pdfName.c_str(), pdfName.c_str(), *ws.var("ctau"), *((RooDataSet*)ws.data(dsName.c_str())), RooKeysPdf::MirrorAsymBoth);
    if (pdf==NULL) { cout << "[ERROR] RooKeysPdf " << pdfName << " is NULL!" << endl; return false; }
    pdf->setNormRange("CtauRecoWindow");
    ws.import(*pdf);
    delete pdf;
  }
  else {
    // Cleaning the input histogram
    // 1) Remove the Under and Overflow bins
    hist->ClearUnderflowAndOverflow();
    // 2) Set negative bin content to zero
    for (int i=0; i<=hist->GetNbinsX(); i++) { if (hist->GetBinContent(i)<0) { hist->SetBinContent(i, 0.0000000001); } }
    // 3) Reduce the range of histogram and rebin it
    TH1* hClean = rebinctauBkghist(ws, hist, range[0], range[1]);
    // 4) Create the RooDataHist for the Template
    cout << Form("[INFO] Implementing %s Template using RooDataHist", pdfName.c_str()) << endl;
    string dataName = pdfName;
    dataName.replace(dataName.find("pdf"), string("pdf").length(), "dh");
    RooDataHist* dataHist = new RooDataHist(dataName.c_str(), "", *ws.var("ctau"), hClean);
    if (dataHist==NULL) { cout << "[ERROR] DataHist used to create " << pdfName << " failed!" << endl; return false; }
    if (dataHist->sumEntries()==0) { cout << "[ERROR] DataHist used to create " << pdfName << " is empty!" << endl; return false; }
    if (fabs(dataHist->sumEntries()-hClean->GetSumOfWeights())>0.001) { cout << "[ERROR] DataHist used to create " << pdfName << "  " << " is invalid!  " << endl; return false; }
    ws.import(*dataHist);
    // 4) Make the Template
    RooHistPdf* pdf = new RooHistPdf(pdfName.c_str(), pdfName.c_str(), *ws.var("ctau"), *((RooDataHist*)ws.data(dataName.c_str())));
    if (pdf==NULL) { cout << "[ERROR] RooHistPdf " << pdfName << " is NULL!" << endl; return false; }
    //RooKeysPdf* pdf = new RooKeysPdf(pdfName.c_str(), pdfName.c_str(), *ws.var("ctau"), *((RooDataSet*)ws.data(dataName.c_str())),RooKeysPdf::NoMirror, 0.4)
    pdf->setNormRange("CtauRecoWindow");
    ws.import(*pdf);
    delete pdf;
    delete dataHist;
  }
  
  return true;
};


TH1* rebinctauhist(TH1 *hist, double xmin, double xmax)
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
    if (hcopy->GetBinContent(i)>0.000000000001) {
      newbins.push_back(hcopy->GetBinLowEdge(i)+hcopy->GetBinWidth(i));
    } else {
      int nrebin=2;
      for (i++; i<=imax; i++) {
        if (hcopy->GetBinContent(i)>0.000000000001) {
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


void setCtauRecoDefaultParameters(map<string, string> &parIni, bool isPbPb, double numEntries)
{
  
  cout << "[INFO] Setting user undefined initial parameters to their default values" << endl;
  
  // DEFAULT RANGE OF NUMBER OF EVENTS
  if (parIni.count(Form("N_Jpsi_%s", (isPbPb?"PbPb":"PP")))==0 || parIni[Form("N_Jpsi_%s", (isPbPb?"PbPb":"PP"))]=="") { 
    parIni[Form("N_Jpsi_%s", (isPbPb?"PbPb":"PP"))]  = Form("%s[%.12f,%.12f,%.12f]", Form("N_Jpsi_%s", (isPbPb?"PbPb":"PP")), numEntries, 0.0, numEntries*2.0);
  }
  if (parIni.count(Form("N_Psi2S_%s", (isPbPb?"PbPb":"PP")))==0 || parIni[Form("N_Psi2S_%s", (isPbPb?"PbPb":"PP"))]=="") { 
    parIni[Form("N_Psi2S_%s", (isPbPb?"PbPb":"PP"))]  = Form("%s[%.12f,%.12f,%.12f]", Form("N_Psi2S_%s", (isPbPb?"PbPb":"PP")), numEntries, 0.0, numEntries*2.0);
  }
  
};


#endif // #ifndef buildCharmoniaCtauRecoModel_C
