#ifndef JPSI_SELECTION_H
#define JPSI_SELECTION_H

#include <cmath>

// single muon acceptance cut
inline bool passAcc(float pt, float eta) {
  const double absEta = std::abs(static_cast<double>(eta));

  return std::isfinite(pt) && std::isfinite(eta) && absEta < 2.4 && pt > (absEta < 1.6 ? 3.5 : 1.5);
}

// single muon selection cut - soft muon
inline bool passQuality(int selectionType, int nTrk, int nPix, float dxy, float dz) {
  // selectionType & 8: turn on bit 3, tracker muon
  return (selectionType & 8) == 8 
    && nTrk > 5 
    && nPix > 0 
    && std::abs(dxy) < 0.3 
    && std::abs(dz) < 20.0;
};

// HLT for dimuon, muons 
// HLT, 2ULL is HLT bit "1" ("1" is second HLT),HLT_PPRefL1DoubleMu0_v*
inline bool passTrigger(unsigned long long qqTrig, unsigned long long muplTrig, unsigned long long mumiTrig) {
  constexpr unsigned long long mask = 2ULL;
  return (qqTrig & mask) == mask
    && (muplTrig & mask) == mask
    && (mumiTrig & mask) == mask;
}

// Reference parent-ID classification: NP keeps B-parent IDs; PR keeps the rest.
inline bool passReferenceTruth(bool isPR, int momId) {
    const int id = std::abs(momId);

    const bool isBParent =
        id == 5 ||
        (id >= 500    && id < 600) ||
        (id >= 5000   && id < 6000) ||
        (id >= 50000  && id < 60000) ||
        (id >= 500000 && id < 600000);

    return isPR ? !isBParent : isBParent;
  }
#endif