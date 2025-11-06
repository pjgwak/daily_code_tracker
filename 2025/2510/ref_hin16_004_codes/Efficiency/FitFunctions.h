/////////////////////////////
// Different Fit functions
// Heavy Ion Dilepton PiNG
// Prashant Shukla and Vineet Kumar
// Use as follows:
// TF1 *fitFunction = new TF1("fitFunction", RBWPol2, 0, 200, 6);
// histMass->Fit("fitFunction","Q0","MER", 20, 120);


Double_t CrystalBall(Double_t *x, Double_t *par)
{
 Double_t N = par[0];
 // Double_t width = par[1];
 Double_t width = 1;

 Double_t mean = par[1];
 Double_t sigma = par[2];
 Double_t alpha = par[3];
 Double_t n = par[4];

 Double_t A = pow(n/fabs(alpha),n)*exp(-pow(alpha,2)/2);
 Double_t B = n/fabs(alpha) - fabs(alpha);

 if ((x[0]-mean)/sigma>-alpha)
  return N*width*TMath::Gaus(x[0],mean,sigma,1);
 else
  return N/(sqrt(TMath::TwoPi())*sigma)*width*A*pow(B-(x[0]-mean)/sigma,-n);
}



//crystall ball
Double_t CrystalBall_1(Double_t *x,Double_t *par) {
  // par[2]->Mean, par[3]->Sigma, par[4]->Yield
  // par[0]-> alpha, par[1] -> n 
  Double_t t = (x[0]-par[2])/par[3];
  if (par[0] < 0) t = -t;
  Double_t absAlpha = fabs((Double_t)par[0]);
  if (t >= -absAlpha) {
    return par[4]*exp(-0.5*t*t);
  }
  else {
    Double_t a =  TMath::Power(par[1]/absAlpha,par[1])*exp(-0.5*absAlpha*absAlpha);
    Double_t b = par[1]/absAlpha - absAlpha;
    return par[4]*(a/TMath::Power(b - t, par[1]));
  }
}




// Polynomials

Double_t Pol1(Double_t *x, Double_t *par)
{
  Double_t pol = par[0] + par[1]*x[0];
  return pol;
}

Double_t Pol11(Double_t *x, Double_t *par)
{
  Double_t pol = par[0] + par[1]*x[0] - par[2]*x[0]*x[0];
  return pol;
}

Double_t Pol2(Double_t *x, Double_t *par)
{
  Double_t pol = par[0] + par[1]*x[0] + par[2]*x[0]*x[0];
  return pol;
}

Double_t Pol3(Double_t *x, Double_t *par)
{
  Double_t pol = par[0] + par[1]*x[0] + par[2]*x[0]*x[0] + par[3]*x[0]*x[0]*x[0];
  return pol;
}

Double_t Pol4(Double_t *x, Double_t *par)
{
  Double_t pol = par[0] + par[1]*x[0] + par[2]*x[0]*x[0] + par[3]*x[0]*x[0]*x[0] + par[4]*x[0]*x[0]*x[0]*x[0];
  return pol;
}


// CB + Pol2
Double_t CBPol2(Double_t *x, Double_t *par)
{
  return CrystalBall(x,par) + Pol2(x,&par[5]);
}


Double_t CB3Pol2(Double_t *x, Double_t *par)
{
  Double_t binwidth = 0.1;
  
  Double_t par1[12], parM[10], parN[10];
  
  par1[0] = binwidth*par[0];
  par1[1] = par[1];
  par1[2] = par[2];
  par1[3] = par[3];
  par1[4] = par[4];
  
  
  parM[0] = binwidth*par[5];
  parM[1] = 1.059*par[1];
  parM[2] = 1.059*par[2];
  parM[3] = par[3];
  parM[4] = par[4];
  
  parN[0] = binwidth*par[6];
  parN[1] = 1.094*par[1];
  parN[2] = 1.094*par[2];
  parN[3] = par[3];
  parN[4] = par[4];
  
  return (CrystalBall(x,par1) + CrystalBall(x, parM) + CrystalBall(x, parN) + Pol2(x, &par[7]));
}

// Exponential

Double_t Exp(Double_t *x, Double_t *par)
{
  Double_t exp = par[0] * TMath::Exp(par[1]*x[0]);
  return exp;
}

Double_t ExpPol2(Double_t *x, Double_t *par)
{
  //Double_t exp = par[0] * TMath::Exp(par[1]*x[0]);
 // Double_t exp = TMath::Exp(par[0] + par[1]*x[0]);
 Double_t exp = TMath::Exp(par[0] + par[1]*x[0] + par[2]*x[0]*x[0]);
  return exp;
  //return Exp(x,par) + Pol2(x,&par[3]);
}

Double_t RBWC(Double_t x, Double_t par1, Double_t par2)
{
  Double_t rbw = x*par1*par2*0.5/TMath::Pi()/TMath::Max(1.e-10,(x*x-par1*par1) * (x*x-par1*par1) + par1*par1*par2*par2);
 return rbw;
}

Double_t GausC(Double_t x, Double_t par1, Double_t par2)
{
  Double_t invsq2pi = 0.3989422804014;   // 1/sqrt(2 pi)
  Double_t arg = 0;
  if (par2) arg = (x - par1)/par2;
  return TMath::Exp(-0.5*arg*arg)*invsq2pi/par2;
}


// Relativitic BreitWignor + BKG

Double_t RBW(Double_t *x, Double_t *par)
{
  Double_t rbw = par[0]*x[0]*par[1]*par[2]*0.5/TMath::Pi()/TMath::Max(1.e-10,(x[0]*x[0]-par[1]*par[1]) * (x[0]*x[0]-par[1]*par[1])+par[1]*par[1]*par[2]*par[2]);
 
  //  Double_t rbw = par[0]*x[0]*par[1]*par[2]/TMath::Max(1.e-10,(x[0]*x[0]-par[1]*par[1]) * (x[0]*x[0]-par[1]*par[1])+par[1]*par[1]*par[2]*par[2]);

 return rbw;

}

Double_t RBWPol1(Double_t *x, Double_t *par)
{
  return RBW(x,par) + Pol1(x,&par[3]);
}

Double_t RBWPol2(Double_t *x, Double_t *par)
{
  return RBW(x,par) + Pol2(x,&par[3]);
}


Double_t RBWPol3(Double_t *x, Double_t *par)
{
  return RBW(x,par) + Pol3(x,&par[3]);
}

Double_t RBWPol4(Double_t *x, Double_t *par)
{
  return RBW(x,par) + Pol4(x,&par[3]);
}

Double_t RBWExp(Double_t *x, Double_t *par)
{
  return RBW(x,par) + Exp(x,&par[3]);

} 


// Non Relativitic BreitWignor + BKG

Double_t BW(Double_t *x, Double_t *par)
{
  Double_t rbw = par[0]*par[2]*0.5/TMath::Pi()/TMath::Max(1.e-10,(x[0]*x[0]-par[1]*par[1]) * (x[0]*x[0]-par[1]*par[1]) + 0.25*par[2]*par[2]);
  return rbw;
}

Double_t BWPol1(Double_t *x, Double_t *par)
{
  return BW(x,par) + Pol1(x,&par[3]);
}

Double_t BWPol2(Double_t *x, Double_t *par)
{
  return BW(x,par) + Pol2(x,&par[3]);
}


Double_t BWPol3(Double_t *x, Double_t *par)
{
  return BW(x,par) + Pol3(x,&par[3]);
}

Double_t BWPol4(Double_t *x, Double_t *par)
{
  return BW(x,par) + Pol4(x,&par[3]);
}

Double_t BWExp(Double_t *x, Double_t *par)
{
  return BW(x,par) + Exp(x,&par[3]);

} 



// Gaussians + BKG

Double_t GausT(Double_t *x, Double_t *par)
{
  Double_t invsq2pi = 0.3989422804014;   // 1/sqrt(2 pi)
  Double_t arg = 0;
  if (par[2]) arg = (x[0] - par[1])/par[2];
  Double_t gau = par[0]*TMath::Exp(-0.5*arg*arg);
  return gau*invsq2pi/par[2];
}

Double_t Gaus(Double_t *x, Double_t *par)
{
  Double_t gaus = par[0]*TMath::Gaus(x[0],par[1],par[2],1);
  return gaus;
}


Double_t GausPol1(Double_t *x, Double_t *par) 
{
  return Gaus(x,par) + Pol1(x,&par[3]);
}

Double_t GausPol2(Double_t *x, Double_t *par) 
{
  return Gaus(x,par) + Pol2(x,&par[3]);
}

Double_t GausPol3(Double_t *x, Double_t *par) 
{
  return Gaus(x,par) + Pol3(x,&par[3]);
}

Double_t GausPol4(Double_t *x, Double_t *par) 
{
  return Gaus(x,par) + Pol4(x,&par[3]);
}

Double_t GausExp(Double_t *x, Double_t *par) 
{
  return Gaus(x,par) + Exp(x,&par[3]);
}


// Three Gaussians 

Double_t GausThree(Double_t *x, Double_t *par) 
{
  double delM1 = 0.5629; // Y(2s) - Y(1s)
  double delM2 = 0.8949; // Y(3s) - Y(1s)
  
  Double_t par2S[3];
  Double_t par3S[3];
  
  par2S[0] = par[3];
  par2S[1] = par[1] + delM1;
  par2S[2] = par[4];

  par3S[0] = par[5];
  par3S[1] = par[1] + delM2;
  par3S[2] = par[6];

  return Gaus(x,par) + Gaus(x, par2S) + Gaus(x,par3S);
}

Double_t CBThree(Double_t *x, Double_t *par) 
{
  return CrystalBall(x,par) + CrystalBall(x, &par[6]) + CrystalBall(x, &par[12]);
}




// Three Gaussians + Pol2 
Double_t GausThreePol2(Double_t *x, Double_t *par) 
{
  return GausThree(x,par) + Pol2(x,&par[7]);
}


Double_t GausTwo(Double_t *x, Double_t *par) 
{
  double delM1 = 0.5629; // Y(2s) - Y(1s)
  double delM2 = 0.8949; // Y(3s) - Y(1s)
  
  Double_t par2S[3];
  
  par2S[0] = par[3];
  par2S[1] = par[1] + (delM1+delM2)/2.0;
  par2S[2] = par[4];

  return Gaus(x,par) + Gaus(x, par2S);
}



// Three Gaussians + Pol2 
Double_t GausTwoPol2(Double_t *x, Double_t *par) 
{
  return GausTwo(x,par) + Pol2(x,&par[5]);
}



Double_t RBWGaus(Double_t *x, Double_t *par) {
  ////// RBW convoluted with gaussian ////// 
  //Fit parameters:
  //par[0]=Normalization
  //par[1]=Mean of RBW
  //par[2]=Width of RBW
  //par[3]=Width (sigma) of convoluted Gaussian function
  
  // Range of convolution integral
  Double_t np = 100.0;      // number of convolution steps
  Double_t xlow = x[0] - 5.0 * par[3];
  Double_t xupp = x[0] + 5.0 * par[3];
  Double_t step = (xupp-xlow) / np;

  // Convolution integral of Landau and Gaussian by sum
  // Variables
  Double_t xx, fland, sum = 0.0;

  for(Double_t i=1.0; i<=np/2; i++) {
    xx = xlow + (i-.5) * step;
    fland = RBWC(xx, par[1], par[2]);
    sum += fland * GausC(x[0], xx, par[3]);
    
    xx = xupp - (i-.5) * step;
    fland = RBWC(xx, par[1], par[2]);
    sum += fland * GausC(x[0], xx, par[3]);
  }
  
  return (par[0] * step * sum);
}


Double_t RBWGausPol1(Double_t *x, Double_t *par)
{
  return RBWGaus(x,par) + Pol1(x,&par[4]);
}


Double_t RBWGausPol2(Double_t *x, Double_t *par)
{
  return RBWGaus(x,par) + Pol2(x,&par[4]);
}


Double_t RBWGausPol3(Double_t *x, Double_t *par)
{
  return RBWGaus(x,par) + Pol3(x,&par[4]);
}


Double_t RBWGausExp(Double_t *x, Double_t *par)
{
  return RBWGaus(x,par) + Exp(x,&par[4]);
}

Double_t abpT(Double_t *x, Double_t *par)
{
  Double_t apT = par[0] + par[1]/x[0];
  return apT;
}
//a+b/pt 
