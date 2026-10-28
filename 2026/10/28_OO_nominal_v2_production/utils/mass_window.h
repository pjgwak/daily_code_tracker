#pragma once
// Dimuon mass fit window (2026-10-27 study hook). Default 2.6-3.5 GeV; per-bin settings mass_lo / mass_hi
// (config/prefit_settings.txt) are exported as OO_MASS_LO / OO_MASS_HI by study/production.py for every stage.
#include "TSystem.h"
#include <cstdlib>
inline double oo_mass_lo() { return gSystem->Getenv("OO_MASS_LO") ? std::atof(gSystem->Getenv("OO_MASS_LO")) : 2.6; }
inline double oo_mass_hi() { return gSystem->Getenv("OO_MASS_HI") ? std::atof(gSystem->Getenv("OO_MASS_HI")) : 3.5; }
