#pragma once

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
 


class myRooGaussian : public RooAbsPdf {
public:
	myRooGaussian() {} ;
	myRooGaussian(const char *name, const char *title,
						RooAbsReal& _x,
						RooAbsReal& _mu,
						RooAbsReal& _sigma);
  
	myRooGaussian(const myRooGaussian& other, const char* name=0) ;
	
	virtual TObject* clone(const char* newname) const
	{
		return new myRooGaussian(*this,newname);
	}

	inline virtual ~myRooGaussian() {;}

	Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const 
	{
		if (matchArgs(allVars,analVars,x)) return 1;
		return 0;
	}
	
	Double_t analyticalIntegral(Int_t code, const char* rangeName) const ;
	
protected:
  RooRealProxy x;
  RooRealProxy mu;
  RooRealProxy sigma;
 
  Double_t evaluate() const ;

private:
  //ClassDef(myRooGaussian,1)
};












 
