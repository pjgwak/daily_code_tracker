// -- CLASS DESCRIPTION [PDF] --
// Gaussian distribution
// 

#include "myRooGaussian.h"
#include "TMath.h"

//ClassImp(myRooJohnsonSU)
  

//_____________________________________________________________________________
myRooGaussian::myRooGaussian(const char *name, const char *title,
			   RooAbsReal& _x,
			   RooAbsReal& _mu,
			   RooAbsReal& _sigma) :
RooAbsPdf(name,title), 
  x("x","Observable",this,_x),
  mu("mu","Mean",this,_mu),
  sigma("sigma","Width",this,_sigma)
 {;}

//_____________________________________________________________________________
myRooGaussian::myRooGaussian(const myRooGaussian& other, const char* name) :
  RooAbsPdf(other,name), 
  x("x",this,other.x),
  mu("mu",this,other.mu),
  sigma("sigma",this,other.sigma)
 {;}

//_____________________________________________________________________________
Double_t myRooGaussian::evaluate() const
{ 
	//if(sigma==0.) return 0.;
	if(sigma<=0.) return 0.;
	double value = TMath::Gaus(x,mu,sigma,false);
	return value;
} 

//_____________________________________________________________________________
Double_t myRooGaussian::analyticalIntegral(Int_t code, const char* rangeName) const
{
	assert(code==1);
	
	Double_t x_min = x.min(rangeName);
	Double_t x_max = x.max(rangeName);
	if (x_max<x_min) return 0.;
	
	//if(sigma==0.) return 0.;
	
	double sqr2s = TMath::Sqrt2()*sigma;
	double norm = TMath::Sqrt(TMath::Pi()/2.)*sigma*(TMath::Erf((mu-x_min)/sqr2s)-TMath::Erf((mu-x_max)/sqr2s));
	return norm;
}



































