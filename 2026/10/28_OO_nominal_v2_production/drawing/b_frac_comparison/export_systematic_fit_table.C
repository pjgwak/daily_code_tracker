#include "../latest_fit_results.h"
#include "TSystem.h"
#include <iostream>
#include <sstream>

// From this folder:
// root -n -l -b -q 'export_systematic_fit_table.C()'       // forward and mid
// root -n -l -b -q 'export_systematic_fit_table.C("mid")'  // mid only
// Tab-separated text for Excel; bin intervals stay text rather than dates.
// bFrac/Error = current systematic-variation fit value/statistical error.
// raw_yield/raw_error = Nsig/statistical error, without bin-width normalization.
// These errors are NOT the systematic shift abs(Bfit - Bnominal).
void export_systematic_fit_table(const char *region = "all",
                                 const char *rootsDir = "", const char *outputDir = "")
{
    const TString selected(region);
    if (selected != "all" && selected != "fwd" && selected != "mid")
        throw std::runtime_error("Region must be all, fwd, or mid");
    TString macro = __FILE__;
    if (!gSystem->IsAbsoluteFileName(macro)) macro = TString(gSystem->WorkingDirectory()) + "/" + macro;
    const TString directory = gSystem->DirName(macro);
    const TString roots = rootsDir[0] ? TString(rootsDir) : directory + "/../../roots";
    const TString output = outputDir[0] ? TString(outputDir) : directory;
    if (gSystem->AccessPathName(output) && gSystem->mkdir(output, true) != 0)
        throw std::runtime_error("Cannot create output directory");

    for (int mid = 0; mid < 2; ++mid) {
        const char *name = mid ? "mid" : "fwd";
        if (selected != "all" && selected != name) continue;
        const TString tag = mid ? "y0.00_1.60" : "y1.60_2.40";
        std::ostringstream table;
        table << "bin\tbFrac\tError\traw_yield\traw_error\n" << std::setprecision(17);
        for (int i = 0; i < (mid ? 8 : 12); ++i) {
            const auto &bin = mid ? nominal_bfraction::midBins()[i] : nominal_bfraction::bins()[i];
            const TString path = Form("%s/%s/fit2d/fit2d_result_%s_pt%.2f_%.2f.root",
                                      roots.Data(), tag.Data(), tag.Data(), bin.low, bin.high);
            const auto fit = oo_drawing::readFit(path);
            table << Form("(%g, %g)", bin.low, bin.high) << '\t'
                  << fit.b << '\t' << fit.bError << '\t' << fit.yield << '\t' << fit.yieldError << '\n';
        }
        const TString path = output + Form("/systematic_fit_%s.txt", name);
        std::ofstream file(path.Data());
        file << table.str();
        file.close();
        if (!file) throw std::runtime_error("Cannot write table: " + std::string(path.Data()));
        std::cout << "\n" << name << ": pT bins (GeV/c), " << tag << '\n'
                  << table.str() << "Saved: " << path << '\n';
    }
}
