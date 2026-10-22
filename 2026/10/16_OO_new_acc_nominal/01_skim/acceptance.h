#ifndef OO_NEW_ACC_ACCEPTANCE_H
#define OO_NEW_ACC_ACCEPTANCE_H

#include "reference/pp_hx_selection.h"

namespace oo_new_acc
{
// Use the exact imported pp_HX selection for both OO converters.
inline bool passMuon(float pt, float eta)
{
  return ::passAcc(pt, eta);
}

constexpr const char *description =
    "pp_HX passAcc: finite pT and eta; |eta|<2.4; pT>3.5 for |eta|<1.6, "
    "otherwise pT>1.5; both muons";
constexpr const char *sourceSha256 = "f8863289063ef5d69a3e0b4b14edf0cb017d1d08b08dcc0b37d3362c76309b14";
}

#endif
