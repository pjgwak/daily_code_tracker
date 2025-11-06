//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////


#ifndef EVENTUTILS_H_
#define EVENTUTILS_H_

#include "TLorentzVector.h"
#include "initOniaTree.C"

namespace HI {
  enum TRIGGERBIT {
    HLT_HIL1DoubleMu0_v1 = 0,
    HLT_HIL1DoubleMu0_2HF_v1 = 1,
    HLT_HIL1DoubleMu0_2HF0_v1 = 2,
    HLT_HIL1DoubleMu10_v1 = 3,
    HLT_HIL2DoubleMu0_NHitQ_v2 = 4,
    HLT_HIL2DoubleMu0_NHitQ_2HF_v1 = 5,
    HLT_HIL2DoubleMu0_NHitQ_2HF0_v1 = 6,
    HLT_HIL1DoubleMu0_2HF_Cent30100_v1 = 7,
    HLT_HIL1DoubleMu0_2HF0_Cent30100_v1 = 8,
    HLT_HIL2DoubleMu0_2HF_Cent30100_NHitQ_v1 = 9,
    HLT_HIL1DoubleMu0_Cent30_v1 = 10,
    HLT_HIL2DoubleMu0_2HF0_Cent30100_NHitQ_v1 = 11,
    HLT_HIL2DoubleMu0_Cent30_NHitQ_v1 = 12,
    HLT_HIL2DoubleMu0_Cent30_OS_NHitQ_v1 = 13,
    HLT_HIL3DoubleMu0_Cent30_v1 = 14,
    HLT_HIL3DoubleMu0_Cent30_OS_m2p5to4p5_v1 = 15,
    HLT_HIL3DoubleMu0_Cent30_OS_m7to14_v1 = 16,
    HLT_HIL3DoubleMu0_OS_m2p5to4p5_v1 = 17,
    HLT_HIL3DoubleMu0_OS_m7to14_v1 = 18,
    HLT_HIL2Mu3_NHitQ10_2HF_v1 = 19,
    HLT_HIL2Mu3_NHitQ10_2HF0_v1 = 20,
    HLT_HIL3Mu3_NHitQ15_2HF_v1 = 21,
    HLT_HIL3Mu3_NHitQ15_2HF0_v1 = 22,
    HLT_HIL2Mu5_NHitQ10_2HF_v1 = 23,
    HLT_HIL2Mu5_NHitQ10_2HF0_v1 = 24,
    HLT_HIL3Mu5_NHitQ15_2HF_v1 = 25,
    HLT_HIL3Mu5_NHitQ15_2HF0_v1 = 26,
    HLT_HIL2Mu7_NHitQ10_2HF_v1 = 27,
    HLT_HIL2Mu7_NHitQ10_2HF0_v1 = 28,
    HLT_HIL3Mu7_NHitQ15_2HF_v1 = 29,
    HLT_HIL3Mu7_NHitQ15_2HF0_v1 = 30,
    HLT_HIL2Mu15_v2 = 31,
    HLT_HIL2Mu15_2HF_v1 = 32,
    HLT_HIL2Mu15_2HF0_v1 = 33,
    HLT_HIL3Mu15_v1 = 34,
    HLT_HIL3Mu15_2HF_v1 = 35,
    HLT_HIL3Mu15_2HF0_v1 = 36,
    HLT_HIL2Mu20_v1 = 37,
    HLT_HIL2Mu20_2HF_v1 = 38,
    HLT_HIL2Mu20_2HF0_v1 = 39,
    HLT_HIL3Mu20_v1 = 40,
    HLT_HIL3Mu20_2HF_v1 = 41,
    HLT_HIL3Mu20_2HF0_v1 = 42
  };
  
  Double_t findNcoll(int hiBin) {
    const int nbins = 200;
    const Double_t Ncoll[nbins] = {1976.95, 1944.02, 1927.29, 1891.9, 1845.3, 1807.2, 1760.45, 1729.18, 1674.8, 1630.3, 1590.52, 1561.72, 1516.1, 1486.5, 1444.68, 1410.88, 1376.4, 1347.32, 1309.71, 1279.98, 1255.31, 1219.89, 1195.13, 1165.96, 1138.92, 1113.37, 1082.26, 1062.42, 1030.6, 1009.96, 980.229, 955.443, 936.501, 915.97, 892.063, 871.289, 847.364, 825.127, 806.584, 789.163, 765.42, 751.187, 733.001, 708.31, 690.972, 677.711, 660.682, 640.431, 623.839, 607.456, 593.307, 576.364, 560.967, 548.909, 530.475, 519.575, 505.105, 490.027, 478.133, 462.372, 451.115, 442.642, 425.76, 416.364, 405.154, 392.688, 380.565, 371.167, 360.28, 348.239, 340.587, 328.746, 320.268, 311.752, 300.742, 292.172, 281.361, 274.249, 267.025, 258.625, 249.931, 240.497, 235.423, 228.63, 219.854, 214.004, 205.425, 199.114, 193.618, 185.644, 180.923, 174.289, 169.641, 161.016, 157.398, 152.151, 147.425, 140.933, 135.924, 132.365, 127.017, 122.127, 117.817, 113.076, 109.055, 105.16, 101.323, 98.098, 95.0548, 90.729, 87.6495, 84.0899, 80.2237, 77.2201, 74.8848, 71.3554, 68.7745, 65.9911, 63.4136, 61.3859, 58.1903, 56.4155, 53.8486, 52.0196, 49.2921, 47.0735, 45.4345, 43.8434, 41.7181, 39.8988, 38.2262, 36.4435, 34.8984, 33.4664, 31.8056, 30.351, 29.2074, 27.6924, 26.7754, 25.4965, 24.2802, 22.9651, 22.0059, 21.0915, 19.9129, 19.1041, 18.1487, 17.3218, 16.5957, 15.5323, 14.8035, 14.2514, 13.3782, 12.8667, 12.2891, 11.61, 11.0026, 10.3747, 9.90294, 9.42648, 8.85324, 8.50121, 7.89834, 7.65197, 7.22768, 6.7755, 6.34855, 5.98336, 5.76555, 5.38056, 5.11024, 4.7748, 4.59117, 4.23247, 4.00814, 3.79607, 3.68702, 3.3767, 3.16309, 2.98282, 2.8095, 2.65875, 2.50561, 2.32516, 2.16357, 2.03235, 1.84061, 1.72628, 1.62305, 1.48916, 1.38784, 1.28366, 1.24693, 1.18552, 1.16085, 1.12596, 1.09298, 1.07402, 1.06105, 1.02954};
    return Ncoll[hiBin];
  };
  
  Double_t findNcollAverage(int hiBinLow, int hiBinHigh) {
    // if the bin corresponds to a bin we already know, give the official number from https://twiki.cern.ch/twiki/bin/view/CMS/GlauberTables
    
    //16-025 bins
    //For centrality dependence in rapidity intervals:
    if (hiBinLow==0&&hiBinHigh==20) return 23.22;
    if (hiBinLow==20&&hiBinHigh==40) return 14.35;
    if (hiBinLow==40&&hiBinHigh==60) return 8.660;
    if (hiBinLow==60&&hiBinHigh==80) return 4.978;
    if (hiBinLow==80&&hiBinHigh==100) return 2.660;
    if (hiBinLow==100&&hiBinHigh==200) return 0.4395;
    
    // For pt dependence in centrality intervals
    if (hiBinLow==20&&hiBinHigh==60) return 11.51;
    if (hiBinLow==60&&hiBinHigh==200) return 1.405;
    
    // For centrality depencence in inclusive rapidity
    if (hiBinLow==0&&hiBinHigh==10) return 1819;
    if (hiBinLow==10&&hiBinHigh==20) return 1432;
    if (hiBinLow==20&&hiBinHigh==30) return 1127;
    if (hiBinLow==30&&hiBinHigh==40) return 882;
    if (hiBinLow==40&&hiBinHigh==50) return 685.9;
    if (hiBinLow==50&&hiBinHigh==60) return 526.5;
    if (hiBinLow==60&&hiBinHigh==70) return 399.4;
    if (hiBinLow==70&&hiBinHigh==80) return 297.6;
    if (hiBinLow==80&&hiBinHigh==90) return 217.2;
    if (hiBinLow==90&&hiBinHigh==100) return 155.5;
    if (hiBinLow==100&&hiBinHigh==120) return 90.72;
    if (hiBinLow==120&&hiBinHigh==140) return 40.1;
    if (hiBinLow==140&&hiBinHigh==200) return 7.663;
    
    //16-004 bins (fwd-y, cause mid-y is the same as 16-025)
    if (hiBinLow==0&&hiBinHigh==40) return 1315;
    if (hiBinLow==40&&hiBinHigh==80) return 477.3;
    if (hiBinLow==80&&hiBinHigh==200) return 56.67;
    
    //Inclusive value
    if (hiBinLow==0&&hiBinHigh==200) return 392.5;
    
    // else... compute the average ncoll ourselves
    cout << "[WARNIGN]: Official Ncoll number not found for range [" << hiBinLow << "," << hiBinHigh << "], returning -1" << endl;
    Double_t w=0;
    const int nbins = 200;
    const Double_t Ncoll[nbins] = {1976.95, 1944.02, 1927.29, 1891.9, 1845.3, 1807.2, 1760.45, 1729.18, 1674.8, 1630.3, 1590.52, 1561.72, 1516.1, 1486.5, 1444.68, 1410.88, 1376.4, 1347.32, 1309.71, 1279.98, 1255.31, 1219.89, 1195.13, 1165.96, 1138.92, 1113.37, 1082.26, 1062.42, 1030.6, 1009.96, 980.229, 955.443, 936.501, 915.97, 892.063, 871.289, 847.364, 825.127, 806.584, 789.163, 765.42, 751.187, 733.001, 708.31, 690.972, 677.711, 660.682, 640.431, 623.839, 607.456, 593.307, 576.364, 560.967, 548.909, 530.475, 519.575, 505.105, 490.027, 478.133, 462.372, 451.115, 442.642, 425.76, 416.364, 405.154, 392.688, 380.565, 371.167, 360.28, 348.239, 340.587, 328.746, 320.268, 311.752, 300.742, 292.172, 281.361, 274.249, 267.025, 258.625, 249.931, 240.497, 235.423, 228.63, 219.854, 214.004, 205.425, 199.114, 193.618, 185.644, 180.923, 174.289, 169.641, 161.016, 157.398, 152.151, 147.425, 140.933, 135.924, 132.365, 127.017, 122.127, 117.817, 113.076, 109.055, 105.16, 101.323, 98.098, 95.0548, 90.729, 87.6495, 84.0899, 80.2237, 77.2201, 74.8848, 71.3554, 68.7745, 65.9911, 63.4136, 61.3859, 58.1903, 56.4155, 53.8486, 52.0196, 49.2921, 47.0735, 45.4345, 43.8434, 41.7181, 39.8988, 38.2262, 36.4435, 34.8984, 33.4664, 31.8056, 30.351, 29.2074, 27.6924, 26.7754, 25.4965, 24.2802, 22.9651, 22.0059, 21.0915, 19.9129, 19.1041, 18.1487, 17.3218, 16.5957, 15.5323, 14.8035, 14.2514, 13.3782, 12.8667, 12.2891, 11.61, 11.0026, 10.3747, 9.90294, 9.42648, 8.85324, 8.50121, 7.89834, 7.65197, 7.22768, 6.7755, 6.34855, 5.98336, 5.76555, 5.38056, 5.11024, 4.7748, 4.59117, 4.23247, 4.00814, 3.79607, 3.68702, 3.3767, 3.16309, 2.98282, 2.8095, 2.65875, 2.50561, 2.32516, 2.16357, 2.03235, 1.84061, 1.72628, 1.62305, 1.48916, 1.38784, 1.28366, 1.24693, 1.18552, 1.16085, 1.12596, 1.09298, 1.07402, 1.06105, 1.02954};
    for(int i=hiBinLow; i<hiBinHigh; i++)  w+=Ncoll[i]/(hiBinHigh-hiBinLow);
    return w;
  };
  
  Double_t findNpartSyst_low(int hiBinLow, int hiBinHigh) {
    // if the bin corresponds to a bin we already know, give the official number from https://twiki.cern.ch/twiki/bin/view/CMS/GlauberTables
    
    //16-025 bins
    //For centrality dependence in rapidity intervals:
    if (hiBinLow==0&&hiBinHigh==20) return 2.6;
    if (hiBinLow==20&&hiBinHigh==40) return 3.8;
    if (hiBinLow==40&&hiBinHigh==60) return 4.1;
    if (hiBinLow==60&&hiBinHigh==80) return 4.0;
    if (hiBinLow==80&&hiBinHigh==100) return 3.7;
    if (hiBinLow==100&&hiBinHigh==200) return 1.0;
    
    // For pt dependence in centrality intervals
    if (hiBinLow==20&&hiBinHigh==60) return 3.9;
    if (hiBinLow==60&&hiBinHigh==200) return 1.2;
    
    // For centrality depencence in inclusive rapidity
    if (hiBinLow==0&&hiBinHigh==10) return 2.0;
    if (hiBinLow==10&&hiBinHigh==20) return 3.2;
    if (hiBinLow==20&&hiBinHigh==30) return 3.7;
    if (hiBinLow==30&&hiBinHigh==40) return 3.9;
    if (hiBinLow==40&&hiBinHigh==50) return 4.1;
    if (hiBinLow==50&&hiBinHigh==60) return 4.0;
    if (hiBinLow==60&&hiBinHigh==70) return 4.0;
    if (hiBinLow==70&&hiBinHigh==80) return 4.0;
    if (hiBinLow==80&&hiBinHigh==90) return 3.8;
    if (hiBinLow==90&&hiBinHigh==100) return 3.6;
    if (hiBinLow==100&&hiBinHigh==120) return 3.1;
    if (hiBinLow==120&&hiBinHigh==140) return 2.4;
    if (hiBinLow==140&&hiBinHigh==200) return 0.6;
    
    //16-004 bins (fwd-y, cause mid-y is the same as 16-025)
    if (hiBinLow==0&&hiBinHigh==40) return 3.1;
    if (hiBinLow==40&&hiBinHigh==80) return 4.0;
    if (hiBinLow==80&&hiBinHigh==200) return 1.1;
    
    //Inclusive value
    if (hiBinLow==0&&hiBinHigh==200) return 2.6;
    
    // else... return a dummy number
    cout << "[WARNIGN]: Official Npart low unc not found for range [" << hiBinLow << "," << hiBinHigh << "], returning -5" << endl;
    return -5;
  };
  
  Double_t findNpartSyst_high(int hiBinLow, int hiBinHigh) {
    // if the bin corresponds to a bin we already know, give the official number from https://twiki.cern.ch/twiki/bin/view/CMS/GlauberTables
    
    //16-025 bins
    //For centrality dependence in rapidity intervals:
    if (hiBinLow==0&&hiBinHigh==20) return 2.4;
    if (hiBinLow==20&&hiBinHigh==40) return 3.6;
    if (hiBinLow==40&&hiBinHigh==60) return 4.0;
    if (hiBinLow==60&&hiBinHigh==80) return 4.0;
    if (hiBinLow==80&&hiBinHigh==100) return 3.7;
    if (hiBinLow==100&&hiBinHigh==200) return 1.8;
    
    // For pt dependence in centrality intervals
    if (hiBinLow==20&&hiBinHigh==60) return 3.7;
    if (hiBinLow==60&&hiBinHigh==200) return 2.4;
    
    // For centrality depencence in inclusive rapidity
    if (hiBinLow==0&&hiBinHigh==10) return 1.8;
    if (hiBinLow==10&&hiBinHigh==20) return 3.0;
    if (hiBinLow==20&&hiBinHigh==30) return 3.5;
    if (hiBinLow==30&&hiBinHigh==40) return 3.8;
    if (hiBinLow==40&&hiBinHigh==50) return 3.9;
    if (hiBinLow==50&&hiBinHigh==60) return 4.0;
    if (hiBinLow==60&&hiBinHigh==70) return 4.0;
    if (hiBinLow==70&&hiBinHigh==80) return 4.0;
    if (hiBinLow==80&&hiBinHigh==90) return 3.8;
    if (hiBinLow==90&&hiBinHigh==100) return 3.7;
    if (hiBinLow==100&&hiBinHigh==120) return 3.2;
    if (hiBinLow==120&&hiBinHigh==140) return 2.6;
    if (hiBinLow==140&&hiBinHigh==200) return 1.0;
    
    //16-004 bins (fwd-y, cause mid-y is the same as 16-025)
    if (hiBinLow==0&&hiBinHigh==40) return 2.9;
    if (hiBinLow==40&&hiBinHigh==80) return 4.0;
    if (hiBinLow==80&&hiBinHigh==200) return 2.1;
    
    //Inclusive value
    if (hiBinLow==0&&hiBinHigh==200) return 2.6;
    
    // else... return a dummy number
    cout << "[WARNIGN]: Official Npart high unc not found for range [" << hiBinLow << "," << hiBinHigh << "], returning -5" << endl;
    return -5;
  };
  
  // Outdated method
//  Double_t findNpartSyst(int hiBinLow, int hiBinHigh) {
//    // if the bin corresponds to a bin we already know, give the official number from https://twiki.cern.ch/twiki/pub/CMS/HiCentrality2016/Centrality_ST_May03v1.pdf
//    cout << "[WARNING] : the Npart uncertainties should be asymetric, please update them" << endl;
//    if (hiBinLow==0&&hiBinHigh==20) return 2.3;
//    if (hiBinLow==20&&hiBinHigh==40) return 2.1;
//    if (hiBinLow==40&&hiBinHigh==60) return 3.6;
//    if (hiBinLow==60&&hiBinHigh==80) return 3.7;
//    if (hiBinLow==80&&hiBinHigh==100) return 2.8;
//    if (hiBinLow==100&&hiBinHigh==200) return 1.4;
//    if (hiBinLow==0&&hiBinHigh==40) return 2.2;
//    if (hiBinLow==40&&hiBinHigh==80) return 2.5;
//    if (hiBinLow==80&&hiBinHigh==200) return 3.3;
//    
//    // else... return a dummy number
//    return 5;
//  };
  
  float findNpart(int hiBin) {
    const int nbins = 200;
    const float Npart[nbins] = {401.99, 398.783, 396.936, 392.71, 387.901, 383.593, 377.914, 374.546, 367.507, 361.252, 356.05, 352.43, 345.701, 341.584, 335.148, 330.581, 325.135, 320.777, 315.074, 310.679, 306.687, 301.189, 296.769, 291.795, 287.516, 283.163, 277.818, 274.293, 269.29, 265.911, 260.574, 256.586, 252.732, 249.194, 245.011, 241.292, 236.715, 232.55, 229.322, 225.328, 221.263, 218.604, 214.728, 210.554, 206.878, 203.924, 200.84, 196.572, 193.288, 189.969, 186.894, 183.232, 180.24, 177.36, 174.008, 171.222, 168.296, 165.319, 162.013, 158.495, 156.05, 154.218, 150.559, 148.455, 145.471, 142.496, 139.715, 137.395, 134.469, 131.926, 129.817, 127.045, 124.467, 122.427, 119.698, 117.607, 114.543, 112.662, 110.696, 108.294, 105.777, 103.544, 101.736, 99.943, 97.4951, 95.4291, 93.2148, 91.2133, 89.5108, 87.2103, 85.7498, 83.5134, 81.9687, 79.7456, 78.1684, 76.4873, 74.7635, 72.761, 71.0948, 69.6102, 67.7806, 66.2215, 64.5813, 63.0269, 61.4325, 59.8065, 58.2423, 57.2432, 55.8296, 54.2171, 52.8809, 51.3254, 49.9902, 48.6927, 47.5565, 46.136, 44.8382, 43.6345, 42.3964, 41.4211, 39.9681, 39.178, 37.9341, 36.9268, 35.5626, 34.5382, 33.6912, 32.8156, 31.6695, 30.6552, 29.7015, 28.8655, 27.9609, 27.0857, 26.105, 25.3163, 24.4872, 23.6394, 23.0484, 22.2774, 21.4877, 20.5556, 19.9736, 19.3296, 18.5628, 17.916, 17.2928, 16.6546, 16.1131, 15.4013, 14.8264, 14.3973, 13.7262, 13.2853, 12.8253, 12.2874, 11.7558, 11.2723, 10.8829, 10.4652, 9.96477, 9.6368, 9.09316, 8.84175, 8.48084, 8.05694, 7.64559, 7.29709, 7.07981, 6.70294, 6.45736, 6.10284, 5.91788, 5.5441, 5.33311, 5.06641, 4.96415, 4.6286, 4.38214, 4.2076, 4.01099, 3.81054, 3.63854, 3.43403, 3.23244, 3.08666, 2.86953, 2.74334, 2.62787, 2.48354, 2.38115, 2.26822, 2.23137, 2.1665, 2.14264, 2.10636, 2.07358, 2.05422, 2.04126, 2.00954};
    return Npart[hiBin];
  };
  
  float findNpartAverage(int hiBinLow, int hiBinHigh) {
    // if the bin corresponds to a bin we already know, give the official number from https://twiki.cern.ch/twiki/bin/view/CMS/GlauberTables
    
    //16-025 bins
    //For centrality dependence in rapidity intervals:
    if (hiBinLow==0&&hiBinHigh==20) return 358.8;
    if (hiBinLow==20&&hiBinHigh==40) return 264.2;
    if (hiBinLow==40&&hiBinHigh==60) return 189.2;
    if (hiBinLow==60&&hiBinHigh==80) return 131.4;
    if (hiBinLow==80&&hiBinHigh==100) return 86.95;
    if (hiBinLow==100&&hiBinHigh==200) return 21.87;
    
    // For pt dependence in centrality intervals
    if (hiBinLow==20&&hiBinHigh==60) return 226.7;
    if (hiBinLow==60&&hiBinHigh==200) return 46.81;
    
    // For centrality depencence in inclusive rapidity
    if (hiBinLow==0&&hiBinHigh==10) return 384.3;
    if (hiBinLow==10&&hiBinHigh==20) return 333.3;
    if (hiBinLow==20&&hiBinHigh==30) return 285.4;
    if (hiBinLow==30&&hiBinHigh==40) return 242.9;
    if (hiBinLow==40&&hiBinHigh==50) return 205.7;
    if (hiBinLow==50&&hiBinHigh==60) return 172.7;
    if (hiBinLow==60&&hiBinHigh==70) return 144.1;
    if (hiBinLow==70&&hiBinHigh==80) return 118.7;
    if (hiBinLow==80&&hiBinHigh==90) return 96.51;
    if (hiBinLow==90&&hiBinHigh==100) return 77.39;
    if (hiBinLow==100&&hiBinHigh==120) return 53.86;
    if (hiBinLow==120&&hiBinHigh==140) return 30.57;
    if (hiBinLow==140&&hiBinHigh==200) return 8.297;
    
    //16-004 bins (fwd-y, cause mid-y is the same as 16-025)
    if (hiBinLow==0&&hiBinHigh==40) return 311.5;
    if (hiBinLow==40&&hiBinHigh==80) return 160.3;
    if (hiBinLow==80&&hiBinHigh==200) return 32.71;
    
    //Inclusive value
    if (hiBinLow==0&&hiBinHigh==200) return 114.0;
    
    // else... compute the average npart ourselves
    
    cout << "[WARNIGN]: Official Npart number not found for range [" << hiBinLow << "," << hiBinHigh << "], so we just compute it" << endl;
    float w=0;
    const int nbins = 200;
    const float Npart[nbins] = {401.99, 398.783, 396.936, 392.71, 387.901, 383.593, 377.914, 374.546, 367.507, 361.252, 356.05, 352.43, 345.701, 341.584, 335.148, 330.581, 325.135, 320.777, 315.074, 310.679, 306.687, 301.189, 296.769, 291.795, 287.516, 283.163, 277.818, 274.293, 269.29, 265.911, 260.574, 256.586, 252.732, 249.194, 245.011, 241.292, 236.715, 232.55, 229.322, 225.328, 221.263, 218.604, 214.728, 210.554, 206.878, 203.924, 200.84, 196.572, 193.288, 189.969, 186.894, 183.232, 180.24, 177.36, 174.008, 171.222, 168.296, 165.319, 162.013, 158.495, 156.05, 154.218, 150.559, 148.455, 145.471, 142.496, 139.715, 137.395, 134.469, 131.926, 129.817, 127.045, 124.467, 122.427, 119.698, 117.607, 114.543, 112.662, 110.696, 108.294, 105.777, 103.544, 101.736, 99.943, 97.4951, 95.4291, 93.2148, 91.2133, 89.5108, 87.2103, 85.7498, 83.5134, 81.9687, 79.7456, 78.1684, 76.4873, 74.7635, 72.761, 71.0948, 69.6102, 67.7806, 66.2215, 64.5813, 63.0269, 61.4325, 59.8065, 58.2423, 57.2432, 55.8296, 54.2171, 52.8809, 51.3254, 49.9902, 48.6927, 47.5565, 46.136, 44.8382, 43.6345, 42.3964, 41.4211, 39.9681, 39.178, 37.9341, 36.9268, 35.5626, 34.5382, 33.6912, 32.8156, 31.6695, 30.6552, 29.7015, 28.8655, 27.9609, 27.0857, 26.105, 25.3163, 24.4872, 23.6394, 23.0484, 22.2774, 21.4877, 20.5556, 19.9736, 19.3296, 18.5628, 17.916, 17.2928, 16.6546, 16.1131, 15.4013, 14.8264, 14.3973, 13.7262, 13.2853, 12.8253, 12.2874, 11.7558, 11.2723, 10.8829, 10.4652, 9.96477, 9.6368, 9.09316, 8.84175, 8.48084, 8.05694, 7.64559, 7.29709, 7.07981, 6.70294, 6.45736, 6.10284, 5.91788, 5.5441, 5.33311, 5.06641, 4.96415, 4.6286, 4.38214, 4.2076, 4.01099, 3.81054, 3.63854, 3.43403, 3.23244, 3.08666, 2.86953, 2.74334, 2.62787, 2.48354, 2.38115, 2.26822, 2.23137, 2.1665, 2.14264, 2.10636, 2.07358, 2.05422, 2.04126, 2.00954};
    for(int i=hiBinLow; i<hiBinHigh; i++)  w+=Npart[i]/(hiBinHigh-hiBinLow);
    return w;
  };
  
  float findTaaAverage(int hiBinLow, int hiBinHigh) {
    // if the bin corresponds to a bin we already know, give the official number from https://twiki.cern.ch/twiki/bin/view/CMS/GlauberTables
    
    //16-025 bins
    //For centrality dependence in rapidity intervals:
    if (hiBinLow==0&&hiBinHigh==20) return 23.22;
    if (hiBinLow==20&&hiBinHigh==40) return 14.35;
    if (hiBinLow==40&&hiBinHigh==60) return 8.660;
    if (hiBinLow==60&&hiBinHigh==80) return 4.978;
    if (hiBinLow==80&&hiBinHigh==100) return 2.660;
    if (hiBinLow==100&&hiBinHigh==200) return 0.4395;
    
    // For pt dependence in centrality intervals
    if (hiBinLow==20&&hiBinHigh==60) return 11.51;
    if (hiBinLow==60&&hiBinHigh==200) return 1.405;
    
    // For centrality depencence in inclusive rapidity
    if (hiBinLow==0&&hiBinHigh==10) return 25.98;
    if (hiBinLow==10&&hiBinHigh==20) return 20.46;
    if (hiBinLow==20&&hiBinHigh==30) return 16.11;
    if (hiBinLow==30&&hiBinHigh==40) return 12.60;
    if (hiBinLow==40&&hiBinHigh==50) return 9.799;
    if (hiBinLow==50&&hiBinHigh==60) return 7.522;
    if (hiBinLow==60&&hiBinHigh==70) return 5.706;
    if (hiBinLow==70&&hiBinHigh==80) return 4.251;
    if (hiBinLow==80&&hiBinHigh==90) return 3.103;
    if (hiBinLow==90&&hiBinHigh==100) return 2.217;
    if (hiBinLow==100&&hiBinHigh==120) return 1.296;
    if (hiBinLow==120&&hiBinHigh==140) return 0.5729;
    if (hiBinLow==140&&hiBinHigh==200) return 0.1095;
    
    //16-004 bins (fwd-y, cause mid-y is the same as 16-025)
    if (hiBinLow==0&&hiBinHigh==40) return 18.79;
    if (hiBinLow==40&&hiBinHigh==80) return 6.819;
    if (hiBinLow==80&&hiBinHigh==200) return 0.8096;
    
    //Inclusive value
    if (hiBinLow==0&&hiBinHigh==200) return 5.607;
    
    // if we don't know the answer... return -1, the user should worry and ask the global observables group!!
    cout << "[WARNIGN]: Official Taa number not found for range [" << hiBinLow << "," << hiBinHigh << "], returning -1" << endl;
    return -1;
  };
  
  float findTaaAverage_err_low(int hiBinLow, int hiBinHigh) {
    // if the bin corresponds to a bin we already know, give the official number from https://twiki.cern.ch/twiki/bin/view/CMS/GlauberTables
    
    //16-025 bins
    //For centrality dependence in rapidity intervals:
    if (hiBinLow==0&&hiBinHigh==20) return 0.689;
    if (hiBinLow==20&&hiBinHigh==40) return 0.455;
    if (hiBinLow==40&&hiBinHigh==60) return 0.332;
    if (hiBinLow==60&&hiBinHigh==80) return 0.242;
    if (hiBinLow==80&&hiBinHigh==100) return 0.176;
    if (hiBinLow==100&&hiBinHigh==200) return 0.032;
    
    // For pt dependence in centrality intervals
    if (hiBinLow==20&&hiBinHigh==60) return 0.388;
    if (hiBinLow==60&&hiBinHigh==200) return 0.061;
    
    // For centrality depencence in inclusive rapidity
    if (hiBinLow==0&&hiBinHigh==10) return 0.766;
    if (hiBinLow==10&&hiBinHigh==20) return 0.605;
    if (hiBinLow==20&&hiBinHigh==30) return 0.496;
    if (hiBinLow==30&&hiBinHigh==40) return 0.427;
    if (hiBinLow==40&&hiBinHigh==50) return 0.374;
    if (hiBinLow==50&&hiBinHigh==60) return 0.316;
    if (hiBinLow==60&&hiBinHigh==70) return 0.270;
    if (hiBinLow==70&&hiBinHigh==80) return 0.236;
    if (hiBinLow==80&&hiBinHigh==90) return 0.190;
    if (hiBinLow==90&&hiBinHigh==100) return 0.157;
    if (hiBinLow==100&&hiBinHigh==120) return 0.115;
    if (hiBinLow==120&&hiBinHigh==140) return 0.064;
    if (hiBinLow==140&&hiBinHigh==200) return 0.011;
    
    //16-004 bins (fwd-y, cause mid-y is the same as 16-025)
    if (hiBinLow==0&&hiBinHigh==40) return 0.563;
    if (hiBinLow==40&&hiBinHigh==80) return 0.283;
    if (hiBinLow==80&&hiBinHigh==200) return 0.046;
    
    //Inclusive value
    if (hiBinLow==0&&hiBinHigh==200) return 0.191;
    
    // if we don't know the answer... return -1, the user should worry and ask the global observables group!!
     cout << "[WARNIGN]: Official Taa lower uncertainty not found for range [" << hiBinLow << "," << hiBinHigh << "], returning -1" << endl;
    return -1;
  };
  
  float findTaaAverage_err_high(int hiBinLow, int hiBinHigh) {
    // if the bin corresponds to a bin we already know, give the official number from https://twiki.cern.ch/twiki/bin/view/CMS/GlauberTables
    
    //16-025 bins
    //For centrality dependence in rapidity intervals:
    if (hiBinLow==0&&hiBinHigh==20) return 0.431;
    if (hiBinLow==20&&hiBinHigh==40) return 0.329;
    if (hiBinLow==40&&hiBinHigh==60) return 0.288;
    if (hiBinLow==60&&hiBinHigh==80) return 0.236;
    if (hiBinLow==80&&hiBinHigh==100) return 0.180;
    if (hiBinLow==100&&hiBinHigh==200) return 0.049;
    
    // For pt dependence in centrality intervals
    if (hiBinLow==20&&hiBinHigh==60) return 0.304;
    if (hiBinLow==60&&hiBinHigh==200) return 0.094;
    
    // For centrality depencence in inclusive rapidity
    if (hiBinLow==0&&hiBinHigh==10) return 0.469;
    if (hiBinLow==10&&hiBinHigh==20) return 0.382;
    if (hiBinLow==20&&hiBinHigh==30) return 0.352;
    if (hiBinLow==30&&hiBinHigh==40) return 0.321;
    if (hiBinLow==40&&hiBinHigh==50) return 0.311;
    if (hiBinLow==50&&hiBinHigh==60) return 0.294;
    if (hiBinLow==60&&hiBinHigh==70) return 0.265;
    if (hiBinLow==70&&hiBinHigh==80) return 0.230;
    if (hiBinLow==80&&hiBinHigh==90) return 0.192;
    if (hiBinLow==90&&hiBinHigh==100) return 0.162;
    if (hiBinLow==100&&hiBinHigh==120) return 0.120;
    if (hiBinLow==120&&hiBinHigh==140) return 0.071;
    if (hiBinLow==140&&hiBinHigh==200) return 0.018;
    
    //16-004 bins (fwd-y, cause mid-y is the same as 16-025)
    if (hiBinLow==0&&hiBinHigh==40) return 0.366;
    if (hiBinLow==40&&hiBinHigh==80) return 0.260;
    if (hiBinLow==80&&hiBinHigh==200) return 0.071;
    
    //Inclusive value
    if (hiBinLow==0&&hiBinHigh==200) return 0.155;
    
    // if we don't know the answer... return -1, the user should worry and ask the global observables group!!
     cout << "[WARNIGN]: Official Taa lower uncertainty not found for range [" << hiBinLow << "," << hiBinHigh << "], returning -1" << endl;
    return -1;
  };
  
  //Outdated method:
//  float findTaaAverage_err(int hiBinLow, int hiBinHigh) {
//    // take it from the official table https://twiki.cern.ch/twiki/pub/CMS/HiCentrality2016/AN-15-080_temp_20161206.pdf
//    if (hiBinLow==0&&hiBinHigh==200) return 0.501;
//    if (hiBinLow==0&&hiBinHigh==20) return 0.392;
//    if (hiBinLow==20&&hiBinHigh==40) return 0.287;
//    if (hiBinLow==40&&hiBinHigh==60) return 0.24;
//    if (hiBinLow==60&&hiBinHigh==80) return 0.201;
//    if (hiBinLow==80&&hiBinHigh==100) return 0.155;
//    if (hiBinLow==100&&hiBinHigh==200) return 0.0642;
//    if (hiBinLow==0&&hiBinHigh==60) return 0.274;
//    if (hiBinLow==60&&hiBinHigh==200) return 0.166;
//    if (hiBinLow==0&&hiBinHigh==10) return 0.431;
//    if (hiBinLow==10&&hiBinHigh==20) return 0.351;
//    if (hiBinLow==20&&hiBinHigh==30) return 0.301;
//    if (hiBinLow==30&&hiBinHigh==40) return 0.269;
//    if (hiBinLow==40&&hiBinHigh==50) return 0.246;
//    if (hiBinLow==50&&hiBinHigh==60) return 0.227;
//    if (hiBinLow==60&&hiBinHigh==70) return 0.208;
//    if (hiBinLow==70&&hiBinHigh==80) return 0.188;
//    if (hiBinLow==80&&hiBinHigh==90) return 0.165;
//    if (hiBinLow==90&&hiBinHigh==100) return 0.141;
//    if (hiBinLow==100&&hiBinHigh==120) return 0.106;
//    if (hiBinLow==120&&hiBinHigh==140) return 0.0629;
//    if (hiBinLow==140&&hiBinHigh==200) return 0.0197;
//    if (hiBinLow==20&&hiBinHigh==60) return 0.274;
//    if (hiBinLow==60&&hiBinHigh==100) return 0.188;
//    if (hiBinLow==100&&hiBinHigh==140) return 0.0893;
//    if (hiBinLow==140&&hiBinHigh==180) return 0.0245;
//    // if we don't know the answer... return -1, the user should worry and ask the global observables group!!
//    return -1;
//  };
};

namespace PP {
  enum TRIGGERBIT {
    HLT_HIL1DoubleMu0_v1 = 0,
    HLT_HIL1DoubleMu10_v1 = 1,
    HLT_HIL2DoubleMu0_NHitQ_v1 = 2,
    HLT_HIL3DoubleMu0_OS_m2p5to4p5_v1 = 3,
    HLT_HIL3DoubleMu0_OS_m7to14_v1 = 4,
    HLT_HIL2Mu3_NHitQ10_v1 = 5,
    HLT_HIL3Mu3_NHitQ15_v1 = 6,
    HLT_HIL2Mu5_NHitQ10_v1 = 7,
    HLT_HIL3Mu5_NHitQ15_v1 = 8,
    HLT_HIL2Mu7_NHitQ10_v1 = 9,
    HLT_HIL3Mu7_NHitQ15_v1 = 10,
    HLT_HIL2Mu15_v1 = 11,
    HLT_HIL3Mu15_v1 = 12,
    HLT_HIL2Mu20_v1 = 13,
    HLT_HIL3Mu20_v1 = 14
  };
};

namespace RecoQQ {
  void iniBranches(TChain* fChain)
  {
    if (fChain->GetBranch("HLTriggers")) fChain->SetBranchStatus("HLTriggers",1);
    if (fChain->GetBranch("Reco_QQ_trig")) fChain->SetBranchStatus("Reco_QQ_trig",1);
    if (fChain->GetBranch("Reco_QQ_VtxProb")) fChain->SetBranchStatus("Reco_QQ_VtxProb",1);
    if (fChain->GetBranch("Reco_QQ_mupl_isGoodMuon")) fChain->SetBranchStatus("Reco_QQ_mupl_isGoodMuon",1);
    if (fChain->GetBranch("Reco_QQ_mumi_isGoodMuon")) fChain->SetBranchStatus("Reco_QQ_mumi_isGoodMuon",1);
    if (fChain->GetBranch("Reco_QQ_mupl_nTrkWMea")) fChain->SetBranchStatus("Reco_QQ_mupl_nTrkWMea",1);
    if (fChain->GetBranch("Reco_QQ_mumi_nTrkWMea")) fChain->SetBranchStatus("Reco_QQ_mumi_nTrkWMea",1);
    if (fChain->GetBranch("Reco_QQ_mupl_nPixWMea")) fChain->SetBranchStatus("Reco_QQ_mupl_nPixWMea",1);
    if (fChain->GetBranch("Reco_QQ_mumi_nPixWMea")) fChain->SetBranchStatus("Reco_QQ_mumi_nPixWMea",1);
    if (fChain->GetBranch("Reco_QQ_mupl_dxy")) fChain->SetBranchStatus("Reco_QQ_mupl_dxy",1);
    if (fChain->GetBranch("Reco_QQ_mumi_dxy")) fChain->SetBranchStatus("Reco_QQ_mumi_dxy",1);
    if (fChain->GetBranch("Reco_QQ_mupl_dz")) fChain->SetBranchStatus("Reco_QQ_mupl_dz",1);
    if (fChain->GetBranch("Reco_QQ_mumi_dz")) fChain->SetBranchStatus("Reco_QQ_mumi_dz",1);
    
    if (fChain->GetBranch("Reco_QQ_mupl_nTrkHits")) fChain->SetBranchStatus("Reco_QQ_mupl_nTrkHits",1);
    if (fChain->GetBranch("Reco_QQ_mumi_nTrkHits")) fChain->SetBranchStatus("Reco_QQ_mumi_nTrkHits",1);
    if (fChain->GetBranch("Reco_QQ_mupl_normChi2_global")) fChain->SetBranchStatus("Reco_QQ_mupl_normChi2_global",1);
    if (fChain->GetBranch("Reco_QQ_mumi_normChi2_global")) fChain->SetBranchStatus("Reco_QQ_mumi_normChi2_global",1);
    if (fChain->GetBranch("Reco_QQ_mupl_normChi2_inner")) fChain->SetBranchStatus("Reco_QQ_mupl_normChi2_inner",1);
    if (fChain->GetBranch("Reco_QQ_mumi_normChi2_inner")) fChain->SetBranchStatus("Reco_QQ_mumi_normChi2_inner",1);
    if (fChain->GetBranch("Reco_QQ_mupl_TrkMuArb")) fChain->SetBranchStatus("Reco_QQ_mupl_TrkMuArb",1);
    if (fChain->GetBranch("Reco_QQ_mumi_TrkMuArb")) fChain->SetBranchStatus("Reco_QQ_mumi_TrkMuArb",1);
    
  }
  
  Bool_t isTriggerMatch (Int_t iRecoQQ, Int_t TriggerBit)
  {
    Bool_t cond = true;
    cond = cond && ( (HLTriggers&((ULong64_t)pow(2, TriggerBit))) == ((ULong64_t)pow(2, TriggerBit)) );
    cond = cond && ( (Reco_QQ_trig[iRecoQQ]&((ULong64_t)pow(2, TriggerBit))) == ((ULong64_t)pow(2, TriggerBit)) );
    return cond;
  };
  
  Bool_t isGlobalMuonInAccept2011 (TLorentzVector* Muon)
  {
    return (fabs(Muon->Eta()) < 2.4 &&
            ((fabs(Muon->Eta()) < 1.0 && Muon->Pt() >= 3.4) ||
             (1.0 <= fabs(Muon->Eta()) && fabs(Muon->Eta()) < 1.5 && Muon->Pt() >= 5.8-2.4*fabs(Muon->Eta())) ||
             (1.5 <= fabs(Muon->Eta()) && Muon->Pt() >= 3.3667-7.0/9.0*fabs(Muon->Eta()))));
  };
  
  Bool_t areMuonsInAcceptance2011 (Int_t iRecoQQ)
  {
    TLorentzVector *RecoQQmupl = (TLorentzVector*) Reco_QQ_mupl_4mom->At(iRecoQQ);
    TLorentzVector *RecoQQmumi = (TLorentzVector*) Reco_QQ_mumi_4mom->At(iRecoQQ);
    return ( isGlobalMuonInAccept2011(RecoQQmupl) && isGlobalMuonInAccept2011(RecoQQmumi) );
  };
  
  Bool_t isGlobalMuonInAccept2015 (TLorentzVector* Muon)
  {
    return (fabs(Muon->Eta()) < 2.4 &&
            ((fabs(Muon->Eta()) < 1.2 && Muon->Pt() >= 3.5) ||
             (1.2 <= fabs(Muon->Eta()) && fabs(Muon->Eta()) < 2.1 && Muon->Pt() >= 5.77-1.89*fabs(Muon->Eta())) ||
             (2.1 <= fabs(Muon->Eta()) && Muon->Pt() >= 1.8)));
  };
  
  Bool_t areMuonsInAcceptance2015 (Int_t iRecoQQ)
  {
    TLorentzVector *RecoQQmupl = (TLorentzVector*) Reco_QQ_mupl_4mom->At(iRecoQQ);
    TLorentzVector *RecoQQmumi = (TLorentzVector*) Reco_QQ_mumi_4mom->At(iRecoQQ);
    return ( isGlobalMuonInAccept2015(RecoQQmupl) && isGlobalMuonInAccept2015(RecoQQmumi) );
  };
  
  Bool_t passQualityCuts2015 (Int_t iRecoQQ)
  {
    Bool_t cond = true;
    
    // cond = cond && (Reco_QQ_mumi_highPurity[iRecoQQ]);
    cond = cond && (Reco_QQ_mumi_isGoodMuon[iRecoQQ]==1);
    cond = cond && (Reco_QQ_mumi_nTrkWMea[iRecoQQ] > 5);
    cond = cond && (Reco_QQ_mumi_nPixWMea[iRecoQQ] > 0);
    cond = cond && (fabs(Reco_QQ_mumi_dxy[iRecoQQ]) < 0.3);
    cond = cond && (fabs(Reco_QQ_mumi_dz[iRecoQQ]) < 20.);
    
    // cond = cond && (Reco_QQ_mupl_highPurity[iRecoQQ]);
    cond = cond && (Reco_QQ_mupl_isGoodMuon[iRecoQQ]==1);
    cond = cond && (Reco_QQ_mupl_nTrkWMea[iRecoQQ] > 5);
    cond = cond && (Reco_QQ_mupl_nPixWMea[iRecoQQ] > 0);
    cond = cond && (fabs(Reco_QQ_mupl_dxy[iRecoQQ]) < 0.3);
    cond = cond && (fabs(Reco_QQ_mupl_dz[iRecoQQ]) < 20.);
    
    cond = cond && (Reco_QQ_VtxProb[iRecoQQ] > 0.01);
    
    return cond;
  };
  
  Bool_t passQualityCuts2011 (Int_t iRecoQQ)
  {
    Bool_t cond = true;
    
    cond = cond && (Reco_QQ_mumi_nTrkHits[iRecoQQ] > 10);
    cond = cond && (Reco_QQ_mumi_normChi2_global[iRecoQQ] < 20.0);
    cond = cond && (Reco_QQ_mumi_normChi2_inner[iRecoQQ] < 4.0);
    cond = cond && (Reco_QQ_mumi_TrkMuArb[iRecoQQ]==1);
    cond = cond && (Reco_QQ_mumi_nPixWMea[iRecoQQ]>0);
    cond = cond && (fabs(Reco_QQ_mumi_dxy[iRecoQQ]) < 3.0);
    cond = cond && (fabs(Reco_QQ_mumi_dz[iRecoQQ]) < 15.0);
    
    cond = cond && (Reco_QQ_mupl_nTrkHits[iRecoQQ] > 10);
    cond = cond && (Reco_QQ_mupl_normChi2_global[iRecoQQ] < 20.0);
    cond = cond && (Reco_QQ_mupl_normChi2_inner[iRecoQQ] < 4.0);
    cond = cond && (Reco_QQ_mupl_TrkMuArb[iRecoQQ]==1);
    cond = cond && (Reco_QQ_mupl_nPixWMea[iRecoQQ]>0);
    cond = cond && (fabs(Reco_QQ_mupl_dxy[iRecoQQ]) < 3.0);
    cond = cond && (fabs(Reco_QQ_mupl_dz[iRecoQQ]) < 15.0);
    
    cond = cond && (Reco_QQ_VtxProb[iRecoQQ] > 0.01);
    
    return cond;
  };
  
}

#endif /* EVENTUTILS_H_ */
