#ifndef bin_h
#define bin_h

#include <utility>
#include <tuple>
#include <iostream>
#include <set>

using namespace std;

// a simple template class to store a bin and overload the equality operator

// define a few common uses of the template class
template <typename T> class bin : public pair<T,T> {
   public:
      bin(T a, T b) : pair<T,T>(b,a) {};
      T low() const {return this->second;}
      T high() const {return this->first;}
};
typedef bin<double> binD;
typedef bin<float> binF;
typedef bin<int> binI;

// associate three such bins to make an analysis bin
class anabin : public tuple<binF,binF,binI> {
   public:
      anabin(float rapmin, float rapmax, float ptmin, float ptmax, int centmin, int centmax) :
         tuple<binF,binF,binI> (binF(rapmin,rapmax), binF(ptmin, ptmax), binI(centmin, centmax)) {};
      binF rapbin() const {return get<0>(*this);};
      binF ptbin() const {return get<1>(*this);};
      binI centbin() const {return get<2>(*this);};
      void setrapbin(binF rapbin) {get<0>(*this) = rapbin;};
      void setptbin(binF ptbin) {get<1>(*this) = ptbin;};
      void setcentbin(binI centbin) {get<2>(*this) = centbin;};
      void print() const {
         cout << "rap=[" << get<0>(*this).low() << "," << get<0>(*this).high() <<
            "], pt=[" << get<1>(*this).low() << "," << get<1>(*this).high() <<
            "], cent=[" << get<2>(*this).low() << "," << get<2>(*this).high() << "]" << endl;
      }
};

set<anabin> allbins16004() {
  set<anabin> ans;

  // HIN-16-004
  ans.insert(anabin(0,1.6,6.5,30,0,200));
  ans.insert(anabin(0,1.6,6.5,30,0,20));
  ans.insert(anabin(0,1.6,6.5,30,20,40));
  ans.insert(anabin(0,1.6,6.5,30,40,60));
  ans.insert(anabin(0,1.6,6.5,30,60,80));
  ans.insert(anabin(0,1.6,6.5,30,80,100));
  ans.insert(anabin(0,1.6,6.5,30,100,200));
  ans.insert(anabin(0,1.6,6.5,9,0,200));
  ans.insert(anabin(0,1.6,9,12,0,200));
  ans.insert(anabin(0,1.6,12,15,0,200));
  ans.insert(anabin(0,1.6,15,20,0,200));
  ans.insert(anabin(0,1.6,20,30,0,200));
  ans.insert(anabin(1.6,2.4,3,30,0,200));
  ans.insert(anabin(1.6,2.4,3,30,0,40));
  ans.insert(anabin(1.6,2.4,3,30,40,80));
  ans.insert(anabin(1.6,2.4,3,30,80,200));
  ans.insert(anabin(1.6,2.4,3,6.5,0,200));
  ans.insert(anabin(1.6,2.4,6.5,12,0,200));
  ans.insert(anabin(1.6,2.4,12,30,0,200));
  // all integrated
  ans.insert(anabin(0,2.4,6.5,50,0,200));
  
  return ans;
}

set<anabin> allbins() {
   set<anabin> ans;

   // centrality dependence
   // thinner binning (for jpsi RAA only)
   // 1 rapidity interval
   ans.insert(anabin(0,2.4,6.5,50,0,10));
   ans.insert(anabin(0,2.4,6.5,50,10,20));
   ans.insert(anabin(0,2.4,6.5,50,20,30));
   ans.insert(anabin(0,2.4,6.5,50,30,40));
   ans.insert(anabin(0,2.4,6.5,50,40,50));
   ans.insert(anabin(0,2.4,6.5,50,50,60));
   ans.insert(anabin(0,2.4,6.5,50,60,70));
   ans.insert(anabin(0,2.4,6.5,50,70,80));
   ans.insert(anabin(0,2.4,6.5,50,80,90));
   ans.insert(anabin(0,2.4,6.5,50,90,100));
   ans.insert(anabin(0,2.4,6.5,50,100,120));
   ans.insert(anabin(0,2.4,6.5,50,120,140));
   ans.insert(anabin(0,2.4,6.5,50,140,200));

   // fwd rapidity and low pt
   ans.insert(anabin(1.8,2.4,3.,6.5,0,20));
   ans.insert(anabin(1.8,2.4,3.,6.5,20,40));
   ans.insert(anabin(1.8,2.4,3.,6.5,40,60));
   ans.insert(anabin(1.8,2.4,3.,6.5,60,80));
   ans.insert(anabin(1.8,2.4,3.,6.5,80,100));
   ans.insert(anabin(1.8,2.4,3.,6.5,100,200));
   // integrated bin in pp
   ans.insert(anabin(1.8,2.4,3.,6.5,0,200));

   // 4 rapidity intervals
   // |y|<0.6
   ans.insert(anabin(0,0.6,6.5,50,0,20));
   ans.insert(anabin(0,0.6,6.5,50,20,40));
   ans.insert(anabin(0,0.6,6.5,50,40,60));
   ans.insert(anabin(0,0.6,6.5,50,60,80));
   ans.insert(anabin(0,0.6,6.5,50,80,100));
   ans.insert(anabin(0,0.6,6.5,50,100,200));
   // 0.6<|y|<1.2
   ans.insert(anabin(0.6,1.2,6.5,50,0,20));
   ans.insert(anabin(0.6,1.2,6.5,50,20,40));
   ans.insert(anabin(0.6,1.2,6.5,50,40,60));
   ans.insert(anabin(0.6,1.2,6.5,50,60,80));
   ans.insert(anabin(0.6,1.2,6.5,50,80,100));
   ans.insert(anabin(0.6,1.2,6.5,50,100,200));
   // 1.2<|y|<1.8
   ans.insert(anabin(1.2,1.8,6.5,50,0,20));
   ans.insert(anabin(1.2,1.8,6.5,50,20,40));
   ans.insert(anabin(1.2,1.8,6.5,50,40,60));
   ans.insert(anabin(1.2,1.8,6.5,50,60,80));
   ans.insert(anabin(1.2,1.8,6.5,50,80,100));
   ans.insert(anabin(1.2,1.8,6.5,50,100,200));
   // 1.8<|y|<2.4
   ans.insert(anabin(1.8,2.4,6.5,50,0,20));
   ans.insert(anabin(1.8,2.4,6.5,50,20,40));
   ans.insert(anabin(1.8,2.4,6.5,50,40,60));
   ans.insert(anabin(1.8,2.4,6.5,50,60,80));
   ans.insert(anabin(1.8,2.4,6.5,50,80,100));
   ans.insert(anabin(1.8,2.4,6.5,50,100,200));
   // integrated bins in pp
   ans.insert(anabin(0,0.6,6.5,50,0,200));
   ans.insert(anabin(0.6,1.2,6.5,50,0,200));
   ans.insert(anabin(1.2,1.8,6.5,50,0,200));
   ans.insert(anabin(1.8,2.4,6.5,50,0,200));
   ans.insert(anabin(0.0,2.4,6.5,50,0,200));


   // pt dependence
   // thinner binning (for jpsi raa only)
   // 1 rapidity interval
   ans.insert(anabin(0,2.4,6.5,7.5,0,200));
   ans.insert(anabin(0,2.4,7.5,8.5,0,200));
   ans.insert(anabin(0,2.4,8.5,9.5,0,200));
   ans.insert(anabin(0,2.4,9.5,11,0,200));
   ans.insert(anabin(0,2.4,11,13,0,200));
   ans.insert(anabin(0,2.4,13,15,0,200));
   ans.insert(anabin(0,2.4,15,17.5,0,200));
   ans.insert(anabin(0,2.4,17.5,20,0,200));
   ans.insert(anabin(0,2.4,20,25,0,200));
   ans.insert(anabin(0,2.4,25,30,0,200));
   ans.insert(anabin(0,2.4,30,35,0,200));
   ans.insert(anabin(0,2.4,35,50,0,200));

   // 4 rapidity intervals
   // |y|<0.6
   ans.insert(anabin(0,0.6,6.5,8.5,0,200));
   ans.insert(anabin(0,0.6,8.5,9.5,0,200));
   ans.insert(anabin(0,0.6,9.5,11,0,200));
   ans.insert(anabin(0,0.6,11,15,0,200));
   ans.insert(anabin(0,0.6,15,20,0,200));
   ans.insert(anabin(0,0.6,20,30,0,200));
//   ans.insert(anabin(0,0.6,30,50,0,200));
   // 0.6<|y|<1.2
   ans.insert(anabin(0.6,1.2,6.5,8.5,0,200));
   ans.insert(anabin(0.6,1.2,8.5,9.5,0,200));
   ans.insert(anabin(0.6,1.2,9.5,11,0,200));
   ans.insert(anabin(0.6,1.2,11,15,0,200));
   ans.insert(anabin(0.6,1.2,15,20,0,200));
   ans.insert(anabin(0.6,1.2,20,30,0,200));
//   ans.insert(anabin(0.6,1.2,30,50,0,200));
   // 0.6<|y|<1.8
   ans.insert(anabin(1.2,1.8,6.5,7.5,0,200));
   ans.insert(anabin(1.2,1.8,7.5,8.5,0,200));
   ans.insert(anabin(1.2,1.8,8.5,9.5,0,200));
   ans.insert(anabin(1.2,1.8,9.5,11,0,200));
   ans.insert(anabin(1.2,1.8,11,15,0,200));
   ans.insert(anabin(1.2,1.8,15,20,0,200));
   ans.insert(anabin(1.2,1.8,20,30,0,200));
//   ans.insert(anabin(1.2,1.8,30,50,0,200));
   // 1.8<|y|<2.4
   ans.insert(anabin(1.8,2.4,3.0,4.5,0,200));
   ans.insert(anabin(1.8,2.4,4.5,5.5,0,200));
   ans.insert(anabin(1.8,2.4,5.5,6.5,0,200));
   ans.insert(anabin(1.8,2.4,6.5,7.5,0,200));
   ans.insert(anabin(1.8,2.4,7.5,8.5,0,200));
   ans.insert(anabin(1.8,2.4,8.5,9.5,0,200));
   ans.insert(anabin(1.8,2.4,9.5,11,0,200));
   ans.insert(anabin(1.8,2.4,11,15,0,200));
   ans.insert(anabin(1.8,2.4,15,20,0,200));
   ans.insert(anabin(1.8,2.4,20,30,0,200));
//   ans.insert(anabin(1.8,2.4,30,50,0,200));

   // 3 centrality intervals
   // 0-10%
   ans.insert(anabin(0,2.4,6.5,7.5,0,20));
   ans.insert(anabin(0,2.4,7.5,8.5,0,20));
   ans.insert(anabin(0,2.4,8.5,9.5,0,20));
   ans.insert(anabin(0,2.4,9.5,11,0,20));
   ans.insert(anabin(0,2.4,11,13,0,20));
   ans.insert(anabin(0,2.4,13,15,0,20));
   ans.insert(anabin(0,2.4,15,20,0,20));
   ans.insert(anabin(0,2.4,20,30,0,20));
//   ans.insert(anabin(0,2.4,30,50,0,20));
   // 10-30%
   ans.insert(anabin(0,2.4,6.5,7.5,20,60));
   ans.insert(anabin(0,2.4,7.5,8.5,20,60));
   ans.insert(anabin(0,2.4,8.5,9.5,20,60));
   ans.insert(anabin(0,2.4,9.5,11,20,60));
   ans.insert(anabin(0,2.4,11,13,20,60));
   ans.insert(anabin(0,2.4,13,15,20,60));
   ans.insert(anabin(0,2.4,15,20,20,60));
   ans.insert(anabin(0,2.4,20,30,20,60));
//   ans.insert(anabin(0,2.4,30,50,20,60));
   // 30-100%
   ans.insert(anabin(0,2.4,6.5,7.5,60,200));
   ans.insert(anabin(0,2.4,7.5,8.5,60,200));
   ans.insert(anabin(0,2.4,8.5,9.5,60,200));
   ans.insert(anabin(0,2.4,9.5,11,60,200));
   ans.insert(anabin(0,2.4,11,13,60,200));
   ans.insert(anabin(0,2.4,13,15,60,200));
   ans.insert(anabin(0,2.4,15,20,60,200));
   ans.insert(anabin(0,2.4,20,30,60,200));
//   ans.insert(anabin(0,2.4,30,50,60,200));
   // 0-100% (for pp)
   ans.insert(anabin(0,2.4,6.5,7.5,0,200));
   ans.insert(anabin(0,2.4,7.5,8.5,0,200));
   ans.insert(anabin(0,2.4,8.5,9.5,0,200));
   ans.insert(anabin(0,2.4,9.5,11,0,200));
   ans.insert(anabin(0,2.4,11,13,0,200));
   ans.insert(anabin(0,2.4,13,15,0,200));
   ans.insert(anabin(0,2.4,15,17.5,0,200));
   ans.insert(anabin(0,2.4,17.5,20,0,200));
   ans.insert(anabin(0,2.4,20,25,0,200));
   ans.insert(anabin(0,2.4,25,30,0,200));
   ans.insert(anabin(0,2.4,30,35,0,200));
   ans.insert(anabin(0,2.4,35,50,0,200));

   // rapidity dependence
   ans.insert(anabin(0.0,0.4,6.5,50,0,200));
   ans.insert(anabin(0.4,0.8,6.5,50,0,200));
   ans.insert(anabin(0.8,1.2,6.5,50,0,200));
   ans.insert(anabin(1.2,1.6,6.5,50,0,200));
   ans.insert(anabin(1.6,2.0,6.5,50,0,200));
   ans.insert(anabin(2.0,2.4,6.5,50,0,200));

   return ans;
};

#endif // #ifndef bin_h
