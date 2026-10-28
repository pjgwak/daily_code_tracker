#ifndef NOMINAL_BFRACTION_H
#define NOMINAL_BFRACTION_H
#include <array>
#include <algorithm>
#include <cmath>
#include <stdexcept>

// Reference (2026-10-28): nominal ver2 = AN-26-048 nominal ctau template fit,
// note_draft/AN_26_048/figs/images/ctau_fit/fit_results_roodecay_{fwdy,midy}_nominal/prompt_fraction_results.root (gNonpromptFraction_vs_pT).
// Earlier reference (25 folder): ../../nominal_bfraction_fwdy.xlsx and the 07/27 mid OO 2D values.
// Used only for reporting and validation; the fit does not include this header.
namespace nominal_bfraction {
constexpr double targetPercent = 2.5;
constexpr double forwardAllowedPercent = 3.0; // User-reviewed forward tolerance.
constexpr double preferredPercent = 2.0;
constexpr double bestSigma = 0.138;
constexpr double maginotSigma = 0.414;
constexpr double lastSigma = 1.0;
// KST - Tag mid references so forward-only sigma criteria cannot leak into mid 8-9/9-10.
struct Bin { double low, high, value, error; bool mid = false; };
inline const std::array<Bin, 12> &bins() {
    static const std::array<Bin, 12> values = {{
        {1, 2, 0.09575, 0.01310},
        {2, 3, 0.11930, 0.00879},
        {3, 4, 0.15179, 0.00846},
        {4, 5, 0.18095, 0.00843},
        {5, 6, 0.19068, 0.00772},
        {6, 7, 0.21654, 0.00955},
        {7, 8, 0.24432, 0.01001},
        {8, 9, 0.25202, 0.01163},
        {9, 10, 0.28555, 0.01563},
        {10, 12, 0.34426, 0.01446},
        {12, 14, 0.41924, 0.02164},
        {14, 20, 0.41358, 0.02118},
    }};
    return values;
}
// Mid reference: ../../../07/27_OO_production_2D_fit_with_ctua2D/drawing/b_frac/bFraction_OO_mid_pp.C
// val_oo_y16 / stat_oo_y16 (OO 2D), validation only; forward references unchanged.
inline const std::array<Bin, 8> &midBins() {
    static const std::array<Bin, 8> values = {{
        {7, 8, 0.24703, 0.00745, true},
        {8, 9, 0.26467, 0.00790, true},
        {9, 10, 0.29272, 0.00950, true},
        {10, 12, 0.31715, 0.00898, true},
        {12, 14, 0.37757, 0.01253, true},
        {14, 16, 0.40657, 0.01693, true},
        {16, 20, 0.48169, 0.01816, true},
        {20, 30, 0.54026, 0.02372, true},
    }};
    return values;
}
inline const Bin &find(double low, double high, double yLow, double yHigh) {
    if (std::abs(yLow-1.6)<1e-6 && std::abs(yHigh-2.4)<1e-6)
        for (const auto &bin : bins())
            if (std::abs(bin.low-low)<1e-6 && std::abs(bin.high-high)<1e-6) return bin;
    // KST - Resolve mid nominal values independently of the forward workbook table.
    if (std::abs(yLow)<1e-6 && std::abs(yHigh-1.6)<1e-6)
        for (const auto &bin : midBins())
            if (std::abs(bin.low-low)<1e-6 && std::abs(bin.high-high)<1e-6) return bin;
    throw std::runtime_error("No nominal b-fraction reference for this bin");
}
inline double differencePercent(double fitted, const Bin &bin) {
    return 100.0 * (fitted/bin.value - 1.0);
}
inline double differenceSigma(double fitted, const Bin &bin) {
    return std::abs(fitted-bin.value)/bin.error;
}
inline double targetForBin(const Bin &bin) {
    return 100.0 * bestSigma * bin.error / bin.value;
}
// The last tier is explicit: validation prints both percent and sigma criteria.
inline double allowedTargetForBin(const Bin &bin) {
    return std::max(bin.mid ? targetPercent : forwardAllowedPercent, 100.0 * lastSigma * bin.error / bin.value);
}
inline const char *tier(double fitted, const Bin &bin) {
    const double sigma = differenceSigma(fitted, bin);
    const double percent = std::abs(differencePercent(fitted, bin));
    if (!std::isfinite(sigma) || !std::isfinite(percent)) return "OUTSIDE";
    if (sigma <= bestSigma) return "BEST";
    if (sigma <= maginotSigma) return "MAGINOT";
    if (percent < preferredPercent) return "PERCENT_2";
    if (percent <= targetPercent) return "PERCENT_2P5";
    if (!bin.mid && percent <= forwardAllowedPercent) return "PERCENT_3";
    if (sigma <= lastSigma) return "LAST_1SIGMA";
    return "OUTSIDE";
}
inline bool agrees(double fitted, const Bin &bin, double target = -1.0) {
    if (target == -1.0) target = allowedTargetForBin(bin);
    return std::isfinite(fitted) && std::isfinite(target) && target > 0
        && std::abs(differencePercent(fitted, bin)) <= target + 1e-12;
}
}
#endif
