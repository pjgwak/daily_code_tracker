#ifndef rappGraph_h
#define rappGraph_h

#include "TGraphErrors.h"
#include <fstream>
#include <vector>
#include <iostream>

using namespace std;

TGraphErrors* rapp_cent_mid_276(const char* prependPath=".") {
   ifstream file(Form("%s/Input/chi_Rapp_cent_276.dat",prependPath));
   double npart, chimin, chimax, dummy;
   vector<double> vnpart, vchi, vdchi;
   while (file.good()) {
      file >> npart >> chimin >> chimax >> dummy >> dummy;
      if (!file.good()) break;
      vnpart.push_back(npart);
      vchi.push_back((chimax+chimin)/2.);
      vdchi.push_back((chimax-chimin)/2.);
   }
   int n = vnpart.size();
   return new TGraphErrors(n,vnpart.data(),vchi.data(),vector<double>(n).data(),vdchi.data());
}
TGraphErrors* rapp_cent_fwd_276(const char* prependPath=".") {
   ifstream file(Form("%s/Input/chi_Rapp_cent_276.dat",prependPath));
   double npart, chimin, chimax, dummy;
   vector<double> vnpart, vchi, vdchi;
   while (file.good()) {
      file >> npart >> dummy >> dummy >> chimin >> chimax;
      if (!file.good()) break;
      vnpart.push_back(npart);
      vchi.push_back((chimax+chimin)/2.);
      vdchi.push_back((chimax-chimin)/2.);
   }
   int n = vnpart.size();
   return new TGraphErrors(n,vnpart.data(),vchi.data(),vector<double>(n).data(),vdchi.data());
}
TGraphErrors* rapp_cent_mid_502(const char* prependPath=".") {
   ifstream file(Form("%s/Input/chi_Rapp_cent_502.dat",prependPath));
   double npart, chimin, chimax, dummy;
   vector<double> vnpart, vchi, vdchi;
   while (file.good()) {
      file >> npart >> chimin >> chimax >> dummy >> dummy;
      if (!file.good()) break;
      vnpart.push_back(npart);
      vchi.push_back((chimax+chimin)/2.);
      vdchi.push_back((chimax-chimin)/2.);
   }
   int n = vnpart.size();
   return new TGraphErrors(n,vnpart.data(),vchi.data(),vector<double>(n).data(),vdchi.data());
}
TGraphErrors* rapp_cent_fwd_502(const char* prependPath=".") {
   ifstream file(Form("%s/Input/chi_Rapp_cent_502.dat",prependPath));
   double npart, chimin, chimax, dummy;
   vector<double> vnpart, vchi, vdchi;
   while (file.good()) {
      file >> npart >> dummy >> dummy >> chimin >> chimax;
      if (!file.good()) break;
      vnpart.push_back(npart);
      vchi.push_back((chimax+chimin)/2.);
      vdchi.push_back((chimax-chimin)/2.);
   }
   int n = vnpart.size();
   return new TGraphErrors(n,vnpart.data(),vchi.data(),vector<double>(n).data(),vdchi.data());
}
TGraphErrors* rapp_pt_mid_502(const char* prependPath=".") {
   ifstream file(Form("%s/Input/chi_Rapp_pt_502_mid.dat",prependPath));
   double pt, raa2min, raa2max, raa1min, raa1max;
   vector<double> vpt, vchi, vdchi;
   while (file.good()) {
      file >> pt >> raa2min >> raa2max >> raa1min >> raa1max;
      if (!file.good()) break;
      vpt.push_back(pt);
      double chimin = raa2min / raa1max;
      double chimax = raa2max / raa1min;
      vchi.push_back((chimax+chimin)/2.);
      vdchi.push_back((chimax-chimin)/2.);
   }
   int n = vpt.size();
   return new TGraphErrors(n,vpt.data(),vchi.data(),vector<double>(n).data(),vdchi.data());
}
TGraphErrors* rapp_pt_fwd_502(const char* prependPath=".") {
   ifstream file(Form("%s/Input/chi_Rapp_pt_502_fwd.dat",prependPath));
   double pt, raa2min, raa2max, raa1min, raa1max;
   vector<double> vpt, vchi, vdchi;
   while (file.good()) {
      file >> pt >> raa2min >> raa2max >> raa1min >> raa1max;
      if (!file.good()) break;
      vpt.push_back(pt);
      double chimin = raa2min / raa1max;
      double chimax = raa2max / raa1min;
      vchi.push_back((chimax+chimin)/2.);
      vdchi.push_back((chimax-chimin)/2.);
   }
   int n = vpt.size();
   return new TGraphErrors(n,vpt.data(),vchi.data(),vector<double>(n).data(),vdchi.data());
}

#endif // ifndef rappGraph_h
