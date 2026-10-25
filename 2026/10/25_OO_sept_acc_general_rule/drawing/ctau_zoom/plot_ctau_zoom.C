#include "../../fit2d.C"

// From this folder: root -n -l -b -q 'plot_ctau_zoom.C()'
// Optional display range: root -n -l -b -q 'plot_ctau_zoom.C(-2,3)'
// Exports two lifetime PDFs; fit2d requires >=99% of each saved PR/NP signal.
// Preserve fit selection and normalization; refine bins to at least
// max(400, 2*sqrt(N_eff)) visible bins (rounded to 50) inside the zoom.
void plot_ctau_zoom(double ctauLow = -2.0, double ctauHigh = 3.0)
{
    TString here = __FILE__;
    if (!gSystem->IsAbsoluteFileName(here)) here = TString(gSystem->WorkingDirectory()) + "/" + here;
    const TString output = gSystem->DirName(here);
    const TString previous = gSystem->WorkingDirectory();
    if (!gSystem->ChangeDirectory(output + "/../.."))
        throw std::runtime_error("Cannot find the analysis directory");
    try {
        for (int pt = 1; pt <= 2; ++pt)
            fit2d(pt, pt+1, 1.6, 2.4, false, true, true, "", output, ctauLow, ctauHigh);
    } catch (...) {
        gSystem->ChangeDirectory(previous);
        throw;
    }
    gSystem->ChangeDirectory(previous);
}
