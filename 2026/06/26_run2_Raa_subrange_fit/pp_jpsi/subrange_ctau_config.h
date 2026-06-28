#ifndef PP_JPSI_SUBRANGE_CTAU_CONFIG_H
#define PP_JPSI_SUBRANGE_CTAU_CONFIG_H

#include <algorithm>
#include <cmath>
#include <vector>

namespace pp_jpsi_subrange
{
// Per-bin manual override for the final ctau fit window.
// The requested ctLow/ctHigh are snapped to the overlapping adaptive ctau
// slices. The actual binned range is written to the output ROOT file as
// actualManualCtauFitLow/High.
struct ManualCtauFitRange
{
	float ptLow = 0.0f;
	float ptHigh = 0.0f;
	float yLow = 0.0f;
	float yHigh = 0.0f;
	double ctLow = 0.0;
	double ctHigh = 0.0;
	bool enabled = false;
	const char *note = "";
};

// Per-bin override for the ctau binning used by the final signal-yield fit.
// This table is exact-bin matched so local pull fixes cannot leak into
// neighboring pT/rapidity bins.
struct CtauBinningBand
{
	double ctLow = 0.0;
	double ctHigh = 0.0;
	double targetWidth = 0.0;
};

struct ManualCtauBinningOverride
{
	float ptLow = 0.0f;
	float ptHigh = 0.0f;
	float yLow = 0.0f;
	float yHigh = 0.0f;
	bool enabled = false;
	const char *note = "";
	std::vector<CtauBinningBand> bands;
};

inline const std::vector<ManualCtauFitRange> &manual_ctau_fit_ranges()
{
	static const std::vector<ManualCtauFitRange> ranges = {
		// Format:
		//   {ptLow, ptHigh, yLow, yHigh, ctLow, ctHigh, enabled, note}
		{3.5f, 6.5f, 1.6f, 2.4f, -1.00, 2.00, true, "trim high tail after monotonic p-value scan"},
		{6.5f, 9.0f, 1.6f, 2.4f, -0.50, 2.25, true, "exclude abnormal 2.5-3.0 tail drop"},
		{9.0f, 12.0f, 1.6f, 2.4f, -0.25, 2.25, true, "restore positive tail for peak-preserving binning"},
		{12.0f, 40.0f, 1.6f, 2.4f, -0.25, 2.25, true, "exclude abnormal 2.5-3.0 tail drop"},
		{3.5f, 40.0f, 1.6f, 2.4f, -1.00, 2.25, true, "restore negative side for peak-preserving binning"},
		{6.5f, 9.0f, 0.0f, 1.6f, -0.50, 2.25, true, "exclude abnormal 2.5-3.0 tail drop"},
		{9.0f, 12.0f, 0.0f, 1.6f, -0.50, 2.25, true, "restore positive tail for peak-preserving binning"},
		{12.0f, 15.0f, 0.0f, 1.6f, -0.15, 2.25, true, "exclude negative-side pull and abnormal 2.5-3.0 tail drop"},
		{15.0f, 20.0f, 0.0f, 1.6f, -0.12, 5.00, true, "exclude negative-side pull"},
		{20.0f, 25.0f, 0.0f, 1.6f, -0.25, 3.00, true, ""},
		{25.0f, 40.0f, 0.0f, 1.6f, -0.10, 2.20, true, "trim high-tail pull"},
		{6.5f, 40.0f, 0.0f, 1.6f, -0.10, 2.25, true, "exclude negative-side pull and abnormal 2.5-3.0 tail drop"},
	};
	return ranges;
}

inline const char *default_ctau_binning_note()
{
	return "monotonic coarser ctau binning; abnormal 2.5-3.0 tail drops are excluded by fit range";
}

inline const std::vector<CtauBinningBand> &default_ctau_binning_bands()
{
	static const std::vector<CtauBinningBand> bands = {
		// ctLow, ctHigh, targetWidth. Width targets are nondecreasing away
		// from zero; the final slice builder also enforces this after snapping.
		{-0.060,  0.060, 0.0200},
		{-0.150, -0.060, 0.0250},
		{ 0.060,  0.150, 0.0250},
		{-0.300, -0.150, 0.0300},
		{ 0.150,  0.300, 0.0300},
		{-0.600, -0.300, 0.0500},
		{ 0.300,  0.600, 0.0500},
		{-1.050, -0.600, 0.0750},
		{ 0.600,  1.050, 0.0750},
		{-2.250, -1.050, 0.1000},
		{ 1.050,  2.250, 0.1000},
		{ 2.250,  5.000, 0.2000},
	};
	return bands;
}

inline bool has_default_ctau_binning()
{
	return !default_ctau_binning_bands().empty();
}

inline const std::vector<ManualCtauBinningOverride> &manual_ctau_binning_overrides()
{
	static const std::vector<ManualCtauBinningOverride> overrides = {
		// Exact-bin exceptions for low-p-value bins. These keep enough bins around
		// the prompt peak while merging only the problematic transition/tail bands.
		{25.0f, 40.0f, 0.0f, 1.6f, true, "native peak-preserving binning", {}},
		{20.0f, 25.0f, 0.0f, 1.6f, true, "native peak-preserving binning", {}},
		{6.5f, 40.0f, 0.0f, 1.6f, true, "native peak-preserving binning", {}},
		{12.0f, 15.0f, 0.0f, 1.6f, true, "peak-preserving ctau binning", {{-0.060, 0.060, 0.0050}, {-0.150, -0.060, 0.0100}, {0.060, 0.150, 0.0100}, {-0.300, -0.150, 0.0200}, {0.150, 0.300, 0.0200}, {-0.600, -0.300, 0.0500}, {0.300, 0.600, 0.0500}, {-1.050, -0.600, 0.0750}, {0.600, 1.050, 0.0750}, {-2.250, -1.050, 0.1000}, {1.050, 2.250, 0.1000}, {2.250, 5.000, 0.2000}}},
		{15.0f, 20.0f, 0.0f, 1.6f, true, "peak-preserving ctau binning", {{-0.060, 0.060, 0.0050}, {-0.150, -0.060, 0.0100}, {0.060, 0.150, 0.0100}, {-0.300, -0.150, 0.0200}, {0.150, 0.300, 0.0200}, {-0.600, -0.300, 0.0500}, {0.300, 0.600, 0.0500}, {-1.050, -0.600, 0.0750}, {0.600, 1.050, 0.0750}, {-2.250, -1.050, 0.1000}, {1.050, 2.250, 0.1000}, {2.250, 5.000, 0.2000}}},
		{9.0f, 12.0f, 0.0f, 1.6f, true, "core plus negative transition and positive shoulder", {{-0.060, 0.060, 0.0200}, {-0.300, -0.150, 0.0250}, {0.150, 0.300, 0.0250}}},
		{9.0f, 12.0f, 1.6f, 2.4f, true, "core plus transition bands", {{-0.060, 0.060, 0.0200}, {-0.300, -0.150, 0.0400}, {0.150, 0.300, 0.0400}}},
		{3.5f, 40.0f, 1.6f, 2.4f, true, "wide transition bands with restored negative side", {{-0.060, 0.060, 0.0200}, {-0.300, -0.150, 0.1500}, {0.150, 0.300, 0.1500}}},
	};
	return overrides;
}

inline bool same_bin(double a, double b)
{
	return std::abs(a - b) < 1e-5;
}

inline const ManualCtauFitRange *find_manual_ctau_fit_range(float ptLow, float ptHigh, float yLow, float yHigh)
{
	for (const auto &range : manual_ctau_fit_ranges())
	{
		if (!range.enabled || !(range.ctLow < range.ctHigh))
			continue;
		if (same_bin(range.ptLow, ptLow) && same_bin(range.ptHigh, ptHigh) &&
			same_bin(range.yLow, yLow) && same_bin(range.yHigh, yHigh))
			return &range;
	}
	return nullptr;
}

inline const ManualCtauBinningOverride *find_manual_ctau_binning_override(float ptLow, float ptHigh, float yLow, float yHigh)
{
	for (const auto &entry : manual_ctau_binning_overrides())
	{
		if (!entry.enabled)
			continue;
		if (same_bin(entry.ptLow, ptLow) && same_bin(entry.ptHigh, ptHigh) &&
			same_bin(entry.yLow, yLow) && same_bin(entry.yHigh, yHigh))
			return &entry;
	}
	return nullptr;
}

inline double target_ctau_bin_width(const ManualCtauBinningOverride *binningOverride, double center, double currentWidth)
{
	double targetWidth = currentWidth;
	const auto &bands = binningOverride ? binningOverride->bands : default_ctau_binning_bands();
	for (const auto &band : bands)
	{
		if (!(band.ctLow < band.ctHigh) || !(band.targetWidth > 0.0))
			continue;
		if (center >= band.ctLow - 1e-9 && center < band.ctHigh - 1e-9)
			targetWidth = std::max(targetWidth, band.targetWidth);
	}
	return targetWidth;
}

} // namespace pp_jpsi_subrange

#endif
