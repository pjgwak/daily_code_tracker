#include "TCanvas.h"
#include "TPad.h"
#include "TH1D.h"
#include "TGraphAsymmErrors.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TColor.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TLine.h"
#include "TStyle.h"

void Plot_MASS_1824_4555_0200()
{
  gStyle->SetEndErrorSize(0);
  
//=========Macro generated from canvas: cMassFig_PbPb/cMassFig
//=========  (Wed Jun 28 13:59:22 2017) by ROOT version6.06/00
   TCanvas *cMassFig_PbPb = new TCanvas("cMassFig_PbPb", "cMassFig",0,22,800,800);
   cMassFig_PbPb->Range(0,0,1,1);
   cMassFig_PbPb->SetFillColor(0);
   cMassFig_PbPb->SetBorderMode(0);
   cMassFig_PbPb->SetBorderSize(2);
   cMassFig_PbPb->SetTickx(1);
   cMassFig_PbPb->SetTicky(1);
   cMassFig_PbPb->SetLeftMargin(0.16);
   cMassFig_PbPb->SetRightMargin(0.04);
   cMassFig_PbPb->SetTopMargin(0.08);
   cMassFig_PbPb->SetBottomMargin(0.12);
   cMassFig_PbPb->SetFrameFillStyle(0);
   cMassFig_PbPb->SetFrameBorderMode(0);
  
// ------------>Primitives in pad: pad1_PbPb
//   TPad *pad1_PbPb = new TPad("pad1_PbPb", "",0,0.23,1,1);
//   pad1_PbPb->Draw();
//   pad1_PbPb->cd();
//   pad1_PbPb->Range(2.42,2.332543,3.545,3.534904);
//   pad1_PbPb->SetFillColor(0);
//   pad1_PbPb->SetBorderMode(0);
//   pad1_PbPb->SetBorderSize(2);
//   pad1_PbPb->SetLogy();
//   pad1_PbPb->SetTickx(1);
//   pad1_PbPb->SetTicky(1);
//   pad1_PbPb->SetLeftMargin(0.16);
//   pad1_PbPb->SetRightMargin(0.04);
//   pad1_PbPb->SetTopMargin(0.08);
//   pad1_PbPb->SetBottomMargin(0.015);
//   pad1_PbPb->SetFrameFillStyle(0);
//   pad1_PbPb->SetFrameBorderMode(0);
//   pad1_PbPb->SetFrameFillStyle(0);
//   pad1_PbPb->SetFrameBorderMode(0);
   
   TH1D *frame_25fab720__5 = new TH1D("frame_25fab720__5","",36,2.6,3.5);
   frame_25fab720__5->SetBinContent(1,1091.736);
   frame_25fab720__5->SetMinimum(224.1707);
   frame_25fab720__5->SetMaximum(2746.091);
   frame_25fab720__5->SetEntries(1);
   frame_25fab720__5->SetDirectory(0);
   frame_25fab720__5->SetStats(0);
   frame_25fab720__5->SetLineStyle(0);
   frame_25fab720__5->SetMarkerStyle(20);
   frame_25fab720__5->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV/c^{2})");
   frame_25fab720__5->GetXaxis()->CenterTitle(true);
   frame_25fab720__5->GetXaxis()->SetTitleOffset(1.25);
   frame_25fab720__5->GetXaxis()->SetLabelFont(42);
   frame_25fab720__5->GetXaxis()->SetLabelOffset(0.015);
   frame_25fab720__5->GetXaxis()->SetLabelSize(0.048);
   frame_25fab720__5->GetXaxis()->SetTitleSize(0.044);
   frame_25fab720__5->GetXaxis()->SetTitleFont(42);
   frame_25fab720__5->GetYaxis()->SetTitle("Events / ( 0.025 GeV/c^{2} )");
   frame_25fab720__5->GetYaxis()->CenterTitle(true);
   frame_25fab720__5->GetYaxis()->SetRangeUser(0,1500);
   frame_25fab720__5->GetYaxis()->SetLabelOffset(0.007);
   frame_25fab720__5->GetYaxis()->SetTitleOffset(1.9);
   frame_25fab720__5->GetYaxis()->SetLabelSize(0.048);
   frame_25fab720__5->GetYaxis()->SetTitleSize(0.044);
   frame_25fab720__5->GetYaxis()->SetTitleFont(42);
   frame_25fab720__5->GetZaxis()->SetLabelFont(42);
   frame_25fab720__5->GetZaxis()->SetLabelOffset(0.007);
   frame_25fab720__5->GetZaxis()->SetLabelSize(0.05);
   frame_25fab720__5->GetZaxis()->SetTitleSize(0.06);
   frame_25fab720__5->GetZaxis()->SetTitleFont(42);
   frame_25fab720__5->Draw("FUNC");
   
   Double_t dOS_fx3004[36] = {
   2.6125,
   2.6375,
   2.6625,
   2.6875,
   2.7125,
   2.7375,
   2.7625,
   2.7875,
   2.8125,
   2.8375,
   2.8625,
   2.8875,
   2.9125,
   2.9375,
   2.9625,
   2.9875,
   3.0125,
   3.0375,
   3.0625,
   3.0875,
   3.1125,
   3.1375,
   3.1625,
   3.1875,
   3.2125,
   3.2375,
   3.2625,
   3.2875,
   3.3125,
   3.3375,
   3.3625,
   3.3875,
   3.4125,
   3.4375,
   3.4625,
   3.4875};
   Double_t dOS_fy3004[36] = {
   407,
   379,
   378,
   398,
   356,
   365,
   371,
   343,
   390,
   377,
   378,
   426,
   390,
   422,
   414,
   463,
   528,
   795,
   969,
   1008,
   933,
   736,
   600,
   453,
   383,
   376,
   330,
   371,
   326,
   342,
   324,
   344,
   297,
   345,
   306,
   288};
   Double_t dOS_felx3004[36] = {
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   Double_t dOS_fely3004[36] = {
   20.17424,
   19.46792,
   19.44222,
   19.94994,
   18.86796,
   19.10497,
   19.26136,
   18.52026,
   19.74842,
   19.41649,
   19.44222,
   20.63977,
   19.74842,
   20.54264,
   20.34699,
   21.51743,
   22.97825,
   28.19574,
   31.12876,
   31.74902,
   30.54505,
   27.12932,
   24.4949,
   21.2838,
   19.57039,
   19.39072,
   18.1659,
   19.26136,
   18.05547,
   18.49324,
   18,
   18.54724,
   17.23369,
   18.57418,
   17.49286,
   16.97056};
   Double_t dOS_fehx3004[36] = {
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   Double_t dOS_fehy3004[36] = {
   20.17424,
   19.46792,
   19.44222,
   19.94994,
   18.86796,
   19.10497,
   19.26136,
   18.52026,
   19.74842,
   19.41649,
   19.44222,
   20.63977,
   19.74842,
   20.54264,
   20.34699,
   21.51743,
   22.97825,
   28.19574,
   31.12876,
   31.74902,
   30.54505,
   27.12932,
   24.4949,
   21.2838,
   19.57039,
   19.39072,
   18.1659,
   19.26136,
   18.05547,
   18.49324,
   18,
   18.54724,
   17.23369,
   18.57418,
   17.49286,
   16.97056};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(36,dOS_fx3004,dOS_fy3004,dOS_felx3004,dOS_fehx3004,dOS_fely3004,dOS_fehy3004);
   grae->SetName("dOS");
   grae->SetTitle("Histogram of dOS_DATA_PbPb_plot__invMass");
   grae->SetFillColor(1);
   grae->SetMarkerStyle(8);
   grae->SetMarkerSize(1.2);
   
   TH1F *Graph_dOS3004 = new TH1F("Graph_dOS3004","Histogram of dOS_DATA_PbPb_plot__invMass",100,2.525,3.575);
   Graph_dOS3004->SetMinimum(194.1575);
   Graph_dOS3004->SetMaximum(1116.621);
   Graph_dOS3004->SetDirectory(0);
   Graph_dOS3004->SetStats(0);
   Graph_dOS3004->SetLineStyle(0);
   Graph_dOS3004->SetMarkerStyle(20);
   Graph_dOS3004->GetXaxis()->SetLabelFont(42);
   Graph_dOS3004->GetXaxis()->SetLabelOffset(0.007);
   Graph_dOS3004->GetXaxis()->SetLabelSize(0.05);
   Graph_dOS3004->GetXaxis()->SetTitleSize(0.06);
   Graph_dOS3004->GetXaxis()->SetTitleOffset(0.9);
   Graph_dOS3004->GetXaxis()->SetTitleFont(42);
   Graph_dOS3004->GetYaxis()->SetLabelFont(42);
   Graph_dOS3004->GetYaxis()->SetLabelOffset(0.007);
   Graph_dOS3004->GetYaxis()->SetLabelSize(0.05);
   Graph_dOS3004->GetYaxis()->SetTitleSize(0.06);
   Graph_dOS3004->GetYaxis()->SetTitleOffset(1.25);
   Graph_dOS3004->GetYaxis()->SetTitleFont(42);
   Graph_dOS3004->GetZaxis()->SetLabelFont(42);
   Graph_dOS3004->GetZaxis()->SetLabelOffset(0.007);
   Graph_dOS3004->GetZaxis()->SetLabelSize(0.05);
   Graph_dOS3004->GetZaxis()->SetTitleSize(0.06);
   Graph_dOS3004->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_dOS3004);
   
   grae->Draw("p");
   
   Double_t BKG_fx6[40] = {
   2.6,
   2.6,
   2.625,
   2.65,
   2.675,
   2.7,
   2.725,
   2.75,
   2.775,
   2.8,
   2.825,
   2.85,
   2.875,
   2.9,
   2.925,
   2.95,
   2.975,
   3,
   3.025,
   3.05,
   3.075,
   3.1,
   3.125,
   3.15,
   3.175,
   3.2,
   3.225,
   3.25,
   3.275,
   3.3,
   3.325,
   3.35,
   3.375,
   3.4,
   3.425,
   3.45,
   3.475,
   3.5,
   3.5,
   3.5};
   Double_t BKG_fy6[40] = {
   0,
   387.6149,
   385.6614,
   383.708,
   381.7545,
   379.8011,
   377.8476,
   375.8942,
   373.9407,
   371.9873,
   370.0338,
   368.0804,
   366.1269,
   364.1735,
   362.22,
   360.2666,
   358.3131,
   356.3597,
   354.4062,
   352.4528,
   350.4993,
   348.5459,
   346.5924,
   344.639,
   342.6855,
   340.7321,
   338.7786,
   336.8252,
   334.8717,
   332.9183,
   330.9648,
   329.0114,
   327.0579,
   325.1045,
   323.151,
   321.1976,
   319.2441,
   317.2907,
   317.2907,
   0};
   TGraph *graph = new TGraph(40,BKG_fx6,BKG_fy6);
   graph->SetName("BKG");
   graph->SetTitle("Projection of pdfCTAUMASS_Tot_PbPb");

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#99ccff");
   graph->SetFillColor(ci);

   ci = TColor::GetColor("#0000ff");
   graph->SetLineColor(ci);
//   graph->SetLineStyle(2);
   graph->SetLineWidth(0); // 3
   graph->SetMarkerStyle(20);
   
   TH1F *Graph_BKG6 = new TH1F("Graph_BKG6","Projection of pdfCTAUMASS_Tot_PbPb",100,2.51,3.59);
   Graph_BKG6->SetMinimum(0.4263764);
   Graph_BKG6->SetMaximum(426.3764);
   Graph_BKG6->SetDirectory(0);
   Graph_BKG6->SetStats(0);
   Graph_BKG6->SetLineStyle(0);
   Graph_BKG6->SetMarkerStyle(20);
   Graph_BKG6->GetXaxis()->SetLabelFont(42);
   Graph_BKG6->GetXaxis()->SetLabelOffset(0.007);
   Graph_BKG6->GetXaxis()->SetLabelSize(0.05);
   Graph_BKG6->GetXaxis()->SetTitleSize(0.06);
   Graph_BKG6->GetXaxis()->SetTitleOffset(0.9);
   Graph_BKG6->GetXaxis()->SetTitleFont(42);
   Graph_BKG6->GetYaxis()->SetLabelFont(42);
   Graph_BKG6->GetYaxis()->SetLabelOffset(0.007);
   Graph_BKG6->GetYaxis()->SetLabelSize(0.05);
   Graph_BKG6->GetYaxis()->SetTitleSize(0.06);
   Graph_BKG6->GetYaxis()->SetTitleOffset(1.25);
   Graph_BKG6->GetYaxis()->SetTitleFont(42);
   Graph_BKG6->GetZaxis()->SetLabelFont(42);
   Graph_BKG6->GetZaxis()->SetLabelOffset(0.007);
   Graph_BKG6->GetZaxis()->SetLabelSize(0.05);
   Graph_BKG6->GetZaxis()->SetTitleSize(0.06);
   Graph_BKG6->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_BKG6);
   
   graph->Draw("f");
   
   Double_t JPSIPR_fx7[230] = {
   2.575,
   2.575,
   2.6,
   2.625,
   2.65,
   2.675,
   2.7,
   2.725,
   2.75,
   2.775,
   2.8,
   2.825,
   2.85,
   2.8625,
   2.875,
   2.8875,
   2.9,
   2.90625,
   2.9125,
   2.91875,
   2.925,
   2.93125,
   2.9375,
   2.940625,
   2.94375,
   2.946875,
   2.95,
   2.953125,
   2.95625,
   2.959375,
   2.9625,
   2.964063,
   2.965625,
   2.967188,
   2.96875,
   2.970312,
   2.971875,
   2.973438,
   2.975,
   2.976562,
   2.978125,
   2.979688,
   2.98125,
   2.982813,
   2.984375,
   2.985937,
   2.9875,
   2.989062,
   2.990625,
   2.992188,
   2.99375,
   2.995312,
   2.996875,
   2.998438,
   3,
   3.001562,
   3.003125,
   3.004688,
   3.00625,
   3.007812,
   3.009375,
   3.010938,
   3.0125,
   3.014063,
   3.015625,
   3.017187,
   3.01875,
   3.020312,
   3.021875,
   3.023438,
   3.025,
   3.026562,
   3.028125,
   3.029688,
   3.03125,
   3.034375,
   3.0375,
   3.04375,
   3.046875,
   3.05,
   3.051562,
   3.053125,
   3.054688,
   3.05625,
   3.057812,
   3.059375,
   3.060938,
   3.0625,
   3.064062,
   3.065625,
   3.067188,
   3.06875,
   3.070312,
   3.071094,
   3.071875,
   3.072656,
   3.073438,
   3.074219,
   3.075,
   3.075781,
   3.076563,
   3.077344,
   3.078125,
   3.078906,
   3.079688,
   3.080469,
   3.08125,
   3.082031,
   3.082813,
   3.083594,
   3.084375,
   3.085156,
   3.085938,
   3.086719,
   3.0875,
   3.088281,
   3.089063,
   3.089844,
   3.090625,
   3.091406,
   3.092188,
   3.092969,
   3.09375,
   3.094531,
   3.095312,
   3.096094,
   3.096875,
   3.097656,
   3.098438,
   3.099219,
   3.1,
   3.100781,
   3.101562,
   3.102344,
   3.103125,
   3.103906,
   3.104688,
   3.10625,
   3.107813,
   3.109375,
   3.110937,
   3.1125,
   3.114062,
   3.115625,
   3.117188,
   3.11875,
   3.120312,
   3.121875,
   3.123438,
   3.125,
   3.128125,
   3.13125,
   3.1375,
   3.140625,
   3.14375,
   3.145312,
   3.146875,
   3.148438,
   3.15,
   3.151562,
   3.153125,
   3.154688,
   3.15625,
   3.157812,
   3.159375,
   3.160937,
   3.1625,
   3.164062,
   3.165625,
   3.167187,
   3.16875,
   3.170312,
   3.171875,
   3.173437,
   3.175,
   3.176562,
   3.178125,
   3.179688,
   3.18125,
   3.182812,
   3.184375,
   3.185938,
   3.1875,
   3.189062,
   3.190625,
   3.192188,
   3.19375,
   3.195312,
   3.196875,
   3.198438,
   3.2,
   3.201563,
   3.203125,
   3.204688,
   3.20625,
   3.207813,
   3.209375,
   3.210938,
   3.2125,
   3.215625,
   3.21875,
   3.221875,
   3.225,
   3.228125,
   3.23125,
   3.234375,
   3.2375,
   3.240625,
   3.24375,
   3.246875,
   3.25,
   3.253125,
   3.25625,
   3.2625,
   3.26875,
   3.275,
   3.28125,
   3.2875,
   3.3,
   3.325,
   3.35,
   3.375,
   3.4,
   3.425,
   3.45,
   3.475,
   3.5,
   3.5,
   3.525,
   3.525};
   Double_t JPSIPR_fy7[230] = {
   0,
   391.674,
   391.674,
   389.8809,
   388.1159,
   386.3862,
   384.7014,
   383.0748,
   381.5254,
   380.0812,
   378.7849,
   377.7045,
   376.9557,
   376.7656,
   376.75,
   376.9694,
   377.5149,
   377.9517,
   378.5296,
   379.2814,
   380.2507,
   381.4961,
   383.099,
   384.0688,
   385.1749,
   386.4401,
   387.8929,
   389.5683,
   391.5107,
   393.7766,
   396.4389,
   397.9477,
   399.594,
   401.395,
   403.3708,
   405.5451,
   407.9459,
   410.6067,
   413.552,
   416.6761,
   419.9526,
   423.3851,
   426.9771,
   430.732,
   434.6528,
   438.7424,
   443.0036,
   447.4387,
   452.05,
   456.8395,
   461.8088,
   466.9594,
   472.2924,
   477.8089,
   483.5095,
   489.3949,
   495.4654,
   501.7214,
   508.1633,
   514.7918,
   521.6076,
   528.6125,
   535.8087,
   543.2002,
   550.7929,
   558.5959,
   566.6226,
   574.8925,
   583.3689,
   592.0122,
   600.8151,
   609.7695,
   618.8663,
   628.0955,
   637.4461,
   656.4623,
   675.8061,
   714.953,
   734.4514,
   753.669,
   763.1134,
   772.415,
   781.5482,
   790.4868,
   799.2045,
   807.6744,
   815.8697,
   823.7637,
   831.3298,
   838.5421,
   845.3749,
   851.8037,
   857.8045,
   860.6372,
   863.3546,
   865.954,
   868.4329,
   870.7888,
   873.0193,
   875.1222,
   877.0955,
   878.9371,
   880.6451,
   882.2178,
   883.6535,
   884.9508,
   886.1081,
   887.1244,
   887.9985,
   888.7295,
   889.3165,
   889.7589,
   890.0562,
   890.208,
   890.2141,
   890.0744,
   889.789,
   889.3581,
   888.7822,
   888.0617,
   887.1973,
   886.1898,
   885.0402,
   883.7496,
   882.3193,
   880.7506,
   879.045,
   877.2042,
   875.23,
   873.1242,
   870.889,
   868.5263,
   866.0384,
   863.4278,
   860.6967,
   857.8479,
   854.8839,
   848.6212,
   841.9315,
   834.8386,
   827.3672,
   819.5428,
   811.3916,
   802.9401,
   794.215,
   785.2431,
   776.0512,
   766.6658,
   757.1129,
   747.4182,
   727.7026,
   707.7097,
   667.5828,
   647.753,
   628.2529,
   618.6611,
   609.1911,
   599.854,
   590.6595,
   581.6169,
   572.7342,
   564.0186,
   555.4764,
   547.1134,
   538.9343,
   530.9433,
   523.1438,
   515.5387,
   508.1302,
   500.92,
   493.9093,
   487.0987,
   480.4886,
   474.0786,
   467.8683,
   461.8566,
   456.0423,
   450.4237,
   444.9989,
   439.7656,
   434.7215,
   429.8638,
   425.1896,
   420.6958,
   416.379,
   412.2357,
   408.2623,
   404.455,
   400.8099,
   397.3229,
   393.99,
   390.8069,
   387.7694,
   384.8731,
   382.1136,
   379.4866,
   376.9877,
   374.6124,
   372.3562,
   368.184,
   364.4362,
   361.0788,
   358.0787,
   355.4041,
   353.0246,
   350.9116,
   349.0378,
   347.3779,
   345.9085,
   344.6079,
   343.4561,
   342.4352,
   341.5286,
   340.0009,
   338.7727,
   337.7652,
   336.9174,
   336.1836,
   334.931,
   332.8051,
   330.8169,
   328.8499,
   326.8855,
   324.9213,
   322.9572,
   320.993,
   319.0289,
   319.0289,
   319.0289,
   0};
   graph = new TGraph(230,JPSIPR_fx7,JPSIPR_fy7);
   graph->SetName("JPSIPR");
   graph->SetTitle("Projection of pdfCTAUMASS_Tot_PbPb");
   graph->SetFillColor(1);
   graph->SetLineColor(2);

   ci = TColor::GetColor("#000099");
   graph->SetLineColor(2);
   graph->SetLineStyle(9);
   graph->SetLineWidth(3);
   graph->SetMarkerStyle(20);
   
   TH1F *Graph_JPSIPR7 = new TH1F("Graph_JPSIPR7","Projection of pdfCTAUMASS_Tot_PbPb",230,2.48,3.62);
   Graph_JPSIPR7->SetMinimum(0.9792355);
   Graph_JPSIPR7->SetMaximum(979.2355);
   Graph_JPSIPR7->SetDirectory(0);
   Graph_JPSIPR7->SetStats(0);
   Graph_JPSIPR7->SetLineStyle(0);
   Graph_JPSIPR7->SetMarkerStyle(20);
   Graph_JPSIPR7->GetXaxis()->SetLabelFont(42);
   Graph_JPSIPR7->GetXaxis()->SetLabelOffset(0.007);
   Graph_JPSIPR7->GetXaxis()->SetLabelSize(0.05);
   Graph_JPSIPR7->GetXaxis()->SetTitleSize(0.06);
   Graph_JPSIPR7->GetXaxis()->SetTitleOffset(0.9);
   Graph_JPSIPR7->GetXaxis()->SetTitleFont(42);
   Graph_JPSIPR7->GetYaxis()->SetLabelFont(42);
   Graph_JPSIPR7->GetYaxis()->SetLabelOffset(0.007);
   Graph_JPSIPR7->GetYaxis()->SetLabelSize(0.05);
   Graph_JPSIPR7->GetYaxis()->SetTitleSize(0.06);
   Graph_JPSIPR7->GetYaxis()->SetTitleOffset(1.25);
   Graph_JPSIPR7->GetYaxis()->SetTitleFont(42);
   Graph_JPSIPR7->GetZaxis()->SetLabelFont(42);
   Graph_JPSIPR7->GetZaxis()->SetLabelOffset(0.007);
   Graph_JPSIPR7->GetZaxis()->SetLabelSize(0.05);
   Graph_JPSIPR7->GetZaxis()->SetTitleSize(0.06);
   Graph_JPSIPR7->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_JPSIPR7);
   
   graph->Draw("l");
   
   Double_t JPSINOPR_fx8[218] = {
   2.575,
   2.575,
   2.6,
   2.625,
   2.65,
   2.675,
   2.7,
   2.725,
   2.75,
   2.775,
   2.8,
   2.825,
   2.85,
   2.8625,
   2.875,
   2.8875,
   2.9,
   2.90625,
   2.9125,
   2.91875,
   2.925,
   2.93125,
   2.9375,
   2.940625,
   2.94375,
   2.946875,
   2.95,
   2.953125,
   2.95625,
   2.959375,
   2.9625,
   2.964063,
   2.965625,
   2.967188,
   2.96875,
   2.970312,
   2.971875,
   2.973438,
   2.975,
   2.976562,
   2.978125,
   2.979688,
   2.98125,
   2.982813,
   2.984375,
   2.985937,
   2.9875,
   2.989062,
   2.990625,
   2.992188,
   2.99375,
   2.995312,
   2.996875,
   2.998438,
   3,
   3.001562,
   3.003125,
   3.004688,
   3.00625,
   3.007812,
   3.009375,
   3.010938,
   3.0125,
   3.014063,
   3.015625,
   3.017187,
   3.01875,
   3.020312,
   3.021875,
   3.023438,
   3.025,
   3.026562,
   3.028125,
   3.03125,
   3.034375,
   3.0375,
   3.04375,
   3.046875,
   3.05,
   3.051562,
   3.053125,
   3.054688,
   3.05625,
   3.057812,
   3.059375,
   3.060938,
   3.0625,
   3.064062,
   3.065625,
   3.067188,
   3.06875,
   3.070312,
   3.071875,
   3.073438,
   3.075,
   3.075781,
   3.076563,
   3.077344,
   3.078125,
   3.078906,
   3.079688,
   3.080469,
   3.08125,
   3.082031,
   3.082813,
   3.083594,
   3.084375,
   3.085156,
   3.085938,
   3.086719,
   3.0875,
   3.088281,
   3.089063,
   3.089844,
   3.090625,
   3.091406,
   3.092188,
   3.092969,
   3.09375,
   3.094531,
   3.095312,
   3.096094,
   3.096875,
   3.097656,
   3.098438,
   3.099219,
   3.1,
   3.101562,
   3.103125,
   3.104688,
   3.10625,
   3.107813,
   3.109375,
   3.110937,
   3.1125,
   3.114062,
   3.115625,
   3.117188,
   3.11875,
   3.120312,
   3.121875,
   3.123438,
   3.125,
   3.128125,
   3.13125,
   3.1375,
   3.140625,
   3.14375,
   3.146875,
   3.148438,
   3.15,
   3.151562,
   3.153125,
   3.154688,
   3.15625,
   3.157812,
   3.159375,
   3.160937,
   3.1625,
   3.164062,
   3.165625,
   3.167187,
   3.16875,
   3.170312,
   3.171875,
   3.173437,
   3.175,
   3.176562,
   3.178125,
   3.179688,
   3.18125,
   3.182812,
   3.184375,
   3.185938,
   3.1875,
   3.189062,
   3.190625,
   3.192188,
   3.19375,
   3.195312,
   3.196875,
   3.198438,
   3.2,
   3.201563,
   3.203125,
   3.204688,
   3.20625,
   3.209375,
   3.2125,
   3.215625,
   3.21875,
   3.221875,
   3.225,
   3.228125,
   3.23125,
   3.234375,
   3.2375,
   3.240625,
   3.24375,
   3.246875,
   3.25,
   3.25625,
   3.2625,
   3.26875,
   3.275,
   3.2875,
   3.3,
   3.325,
   3.35,
   3.375,
   3.4,
   3.425,
   3.45,
   3.475,
   3.5,
   3.5,
   3.525,
   3.525};
   Double_t JPSINOPR_fy8[218] = {
   0,
   390.2546,
   390.2546,
   388.336,
   386.425,
   384.5234,
   382.6337,
   380.7596,
   378.9061,
   377.0806,
   375.2945,
   373.5661,
   371.9261,
   371.1552,
   370.4309,
   369.7692,
   369.1945,
   368.951,
   368.745,
   368.5854,
   368.4839,
   368.4559,
   368.5233,
   368.6019,
   368.7169,
   368.8743,
   369.0817,
   369.3485,
   369.6864,
   370.1107,
   370.6407,
   370.9531,
   371.3021,
   371.6924,
   372.1293,
   372.6192,
   373.1694,
   373.789,
   374.4845,
   375.2277,
   376.0115,
   376.8369,
   377.7048,
   378.6162,
   379.5719,
   380.5725,
   381.6189,
   382.7117,
   383.8515,
   385.0388,
   386.2741,
   387.5577,
   388.89,
   390.2712,
   391.7015,
   393.1811,
   394.71,
   396.2885,
   397.9165,
   399.5942,
   401.322,
   403.1001,
   404.9293,
   406.8105,
   408.7455,
   410.7365,
   412.7871,
   414.9027,
   417.0733,
   419.2884,
   421.546,
   423.8441,
   426.1801,
   430.9552,
   435.8467,
   440.8256,
   450.9058,
   455.9259,
   460.8711,
   463.2998,
   465.6905,
   468.0363,
   470.3302,
   472.5651,
   474.7339,
   476.8296,
   478.8448,
   480.7727,
   482.6061,
   484.3384,
   485.9629,
   487.4733,
   488.8634,
   490.1278,
   491.2609,
   491.7768,
   492.258,
   492.7041,
   493.1147,
   493.4891,
   493.827,
   494.1279,
   494.3916,
   494.6176,
   494.8057,
   494.9556,
   495.0672,
   495.1402,
   495.1744,
   495.1699,
   495.1265,
   495.0443,
   494.9231,
   494.7632,
   494.5646,
   494.3274,
   494.0519,
   493.7382,
   493.3866,
   492.9974,
   492.5709,
   492.1075,
   491.6077,
   491.0717,
   490.5002,
   489.8936,
   489.2524,
   487.8688,
   486.3542,
   484.7139,
   482.9537,
   481.0795,
   479.0979,
   477.0153,
   474.8385,
   472.5746,
   470.2306,
   467.8136,
   465.3309,
   462.7894,
   460.1963,
   457.5586,
   454.8831,
   449.445,
   443.9329,
   432.8712,
   427.4026,
   422.022,
   416.7583,
   414.1781,
   411.636,
   409.1343,
   406.6753,
   404.2609,
   401.8927,
   399.5723,
   397.301,
   395.0798,
   392.9097,
   390.7914,
   388.7256,
   386.7127,
   384.7529,
   382.8465,
   380.9936,
   379.1941,
   377.4478,
   375.7545,
   374.1138,
   372.5254,
   370.9886,
   369.5029,
   368.0676,
   366.6821,
   365.3455,
   364.057,
   362.8157,
   361.6207,
   360.4709,
   359.3655,
   358.3034,
   357.2834,
   356.3045,
   355.3656,
   354.4655,
   353.603,
   352.7771,
   351.23,
   349.8147,
   348.522,
   347.3424,
   346.267,
   345.2868,
   344.3935,
   343.5788,
   342.8353,
   342.1555,
   341.5328,
   340.9609,
   340.4339,
   339.9467,
   339.0726,
   338.3051,
   337.6174,
   336.9887,
   335.8467,
   334.7925,
   332.7852,
   330.8146,
   328.8497,
   326.8855,
   324.9213,
   322.9572,
   320.993,
   319.0289,
   319.0289,
   319.0289,
   0};
   graph = new TGraph(218,JPSINOPR_fx8,JPSINOPR_fy8);
   graph->SetName("JPSINOPR");
   graph->SetTitle("Projection of pdfCTAUMASS_Tot_PbPb");
   graph->SetFillColor(1);

   ci = TColor::GetColor("#339933");
   graph->SetLineColor(ci);
   graph->SetLineStyle(4);
   graph->SetLineWidth(5);
   graph->SetMarkerStyle(20);
   
   TH1F *Graph_JPSINOPR8 = new TH1F("Graph_JPSINOPR8","Projection of pdfCTAUMASS_Tot_PbPb",218,2.48,3.62);
   Graph_JPSINOPR8->SetMinimum(0.5446919);
   Graph_JPSINOPR8->SetMaximum(544.6919);
   Graph_JPSINOPR8->SetDirectory(0);
   Graph_JPSINOPR8->SetStats(0);
   Graph_JPSINOPR8->SetLineStyle(0);
   Graph_JPSINOPR8->SetMarkerStyle(20);
   Graph_JPSINOPR8->GetXaxis()->SetLabelFont(42);
   Graph_JPSINOPR8->GetXaxis()->SetLabelOffset(0.007);
   Graph_JPSINOPR8->GetXaxis()->SetLabelSize(0.05);
   Graph_JPSINOPR8->GetXaxis()->SetTitleSize(0.06);
   Graph_JPSINOPR8->GetXaxis()->SetTitleOffset(0.9);
   Graph_JPSINOPR8->GetXaxis()->SetTitleFont(42);
   Graph_JPSINOPR8->GetYaxis()->SetLabelFont(42);
   Graph_JPSINOPR8->GetYaxis()->SetLabelOffset(0.007);
   Graph_JPSINOPR8->GetYaxis()->SetLabelSize(0.05);
   Graph_JPSINOPR8->GetYaxis()->SetTitleSize(0.06);
   Graph_JPSINOPR8->GetYaxis()->SetTitleOffset(1.25);
   Graph_JPSINOPR8->GetYaxis()->SetTitleFont(42);
   Graph_JPSINOPR8->GetZaxis()->SetLabelFont(42);
   Graph_JPSINOPR8->GetZaxis()->SetLabelOffset(0.007);
   Graph_JPSINOPR8->GetZaxis()->SetLabelSize(0.05);
   Graph_JPSINOPR8->GetZaxis()->SetTitleSize(0.06);
   Graph_JPSINOPR8->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_JPSINOPR8);
   
   graph->Draw("l");
   
   Double_t dOS_fx3005[36] = {
   2.6125,
   2.6375,
   2.6625,
   2.6875,
   2.7125,
   2.7375,
   2.7625,
   2.7875,
   2.8125,
   2.8375,
   2.8625,
   2.8875,
   2.9125,
   2.9375,
   2.9625,
   2.9875,
   3.0125,
   3.0375,
   3.0625,
   3.0875,
   3.1125,
   3.1375,
   3.1625,
   3.1875,
   3.2125,
   3.2375,
   3.2625,
   3.2875,
   3.3125,
   3.3375,
   3.3625,
   3.3875,
   3.4125,
   3.4375,
   3.4625,
   3.4875};
   Double_t dOS_fy3005[36] = {
   407,
   379,
   378,
   398,
   356,
   365,
   371,
   343,
   390,
   377,
   378,
   426,
   390,
   422,
   414,
   463,
   528,
   795,
   969,
   1008,
   933,
   736,
   600,
   453,
   383,
   376,
   330,
   371,
   326,
   342,
   324,
   344,
   297,
   345,
   306,
   288};
   Double_t dOS_felx3005[36] = {
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   Double_t dOS_fely3005[36] = {
   20.17424,
   19.46792,
   19.44222,
   19.94994,
   18.86796,
   19.10497,
   19.26136,
   18.52026,
   19.74842,
   19.41649,
   19.44222,
   20.63977,
   19.74842,
   20.54264,
   20.34699,
   21.51743,
   22.97825,
   28.19574,
   31.12876,
   31.74902,
   30.54505,
   27.12932,
   24.4949,
   21.2838,
   19.57039,
   19.39072,
   18.1659,
   19.26136,
   18.05547,
   18.49324,
   18,
   18.54724,
   17.23369,
   18.57418,
   17.49286,
   16.97056};
   Double_t dOS_fehx3005[36] = {
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   Double_t dOS_fehy3005[36] = {
   20.17424,
   19.46792,
   19.44222,
   19.94994,
   18.86796,
   19.10497,
   19.26136,
   18.52026,
   19.74842,
   19.41649,
   19.44222,
   20.63977,
   19.74842,
   20.54264,
   20.34699,
   21.51743,
   22.97825,
   28.19574,
   31.12876,
   31.74902,
   30.54505,
   27.12932,
   24.4949,
   21.2838,
   19.57039,
   19.39072,
   18.1659,
   19.26136,
   18.05547,
   18.49324,
   18,
   18.54724,
   17.23369,
   18.57418,
   17.49286,
   16.97056};
   grae = new TGraphAsymmErrors(36,dOS_fx3005,dOS_fy3005,dOS_felx3005,dOS_fehx3005,dOS_fely3005,dOS_fehy3005);
   grae->SetName("dOS");
   grae->SetTitle("Histogram of dOS_DATA_PbPb_plot__invMass");
   grae->SetFillColor(1);
   grae->SetMarkerStyle(8);
   grae->SetMarkerSize(1.2);
   
   TH1F *Graph_dOS3005 = new TH1F("Graph_dOS3005","Histogram of dOS_DATA_PbPb_plot__invMass",100,2.525,3.575);
   Graph_dOS3005->SetMinimum(194.1575);
   Graph_dOS3005->SetMaximum(1116.621);
   Graph_dOS3005->SetDirectory(0);
   Graph_dOS3005->SetStats(0);
   Graph_dOS3005->SetLineStyle(0);
   Graph_dOS3005->SetMarkerStyle(20);
   Graph_dOS3005->GetXaxis()->SetLabelFont(42);
   Graph_dOS3005->GetXaxis()->SetLabelOffset(0.007);
   Graph_dOS3005->GetXaxis()->SetLabelSize(0.05);
   Graph_dOS3005->GetXaxis()->SetTitleSize(0.06);
   Graph_dOS3005->GetXaxis()->SetTitleOffset(0.9);
   Graph_dOS3005->GetXaxis()->SetTitleFont(42);
   Graph_dOS3005->GetYaxis()->SetLabelFont(42);
   Graph_dOS3005->GetYaxis()->SetLabelOffset(0.007);
   Graph_dOS3005->GetYaxis()->SetLabelSize(0.05);
   Graph_dOS3005->GetYaxis()->SetTitleSize(0.06);
   Graph_dOS3005->GetYaxis()->SetTitleOffset(1.25);
   Graph_dOS3005->GetYaxis()->SetTitleFont(42);
   Graph_dOS3005->GetZaxis()->SetLabelFont(42);
   Graph_dOS3005->GetZaxis()->SetLabelOffset(0.007);
   Graph_dOS3005->GetZaxis()->SetLabelSize(0.05);
   Graph_dOS3005->GetZaxis()->SetTitleSize(0.06);
   Graph_dOS3005->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_dOS3005);
   
   grae->Draw("p");
   
   Double_t PDF_fx9[86] = {
   2.575,
   2.575,
   2.6,
   2.625,
   2.65,
   2.675,
   2.7,
   2.725,
   2.75,
   2.775,
   2.8,
   2.825,
   2.85,
   2.875,
   2.9,
   2.925,
   2.9375,
   2.95,
   2.95625,
   2.9625,
   2.96875,
   2.975,
   2.98125,
   2.9875,
   2.99375,
   3,
   3.00625,
   3.0125,
   3.01875,
   3.025,
   3.03125,
   3.0375,
   3.05,
   3.05625,
   3.059375,
   3.0625,
   3.065625,
   3.06875,
   3.071875,
   3.075,
   3.078125,
   3.08125,
   3.084375,
   3.0875,
   3.090625,
   3.09375,
   3.096875,
   3.1,
   3.103125,
   3.10625,
   3.109375,
   3.1125,
   3.115625,
   3.11875,
   3.125,
   3.1375,
   3.14375,
   3.15,
   3.15625,
   3.1625,
   3.16875,
   3.175,
   3.18125,
   3.1875,
   3.19375,
   3.2,
   3.20625,
   3.2125,
   3.21875,
   3.225,
   3.2375,
   3.25,
   3.2625,
   3.275,
   3.3,
   3.325,
   3.35,
   3.375,
   3.4,
   3.425,
   3.45,
   3.475,
   3.5,
   3.5,
   3.525,
   3.525};
   Double_t PDF_fy9[86] = {
   0,
   392.1903,
   392.1903,
   390.4427,
   388.7309,
   387.0637,
   385.4533,
   383.9168,
   382.4781,
   381.1725,
   380.0543,
   379.2096,
   378.785,
   379.0482,
   380.5409,
   384.5302,
   388.4,
   394.7343,
   399.448,
   405.8215,
   414.733,
   427.7605,
   444.8969,
   465.3285,
   489.2799,
   516.8991,
   548.2589,
   583.4081,
   622.571,
   666.0134,
   712.5446,
   761.266,
   860.1565,
   906.9244,
   928.7612,
   949.207,
   967.9922,
   984.856,
   999.5531,
   1011.861,
   1021.586,
   1028.571,
   1032.701,
   1033.903,
   1032.155,
   1027.48,
   1019.952,
   1009.686,
   996.8411,
   981.6106,
   964.2177,
   944.9081,
   923.943,
   901.5918,
   853.8101,
   752.9449,
   703.2569,
   655.7685,
   611.3332,
   570.5086,
   533.6083,
   500.7533,
   471.9156,
   446.9543,
   425.6435,
   407.6958,
   392.783,
   380.5543,
   370.653,
   362.7309,
   351.5408,
   344.7325,
   340.6176,
   338.0476,
   334.9814,
   332.8124,
   330.8178,
   328.85,
   326.8855,
   324.9213,
   322.9572,
   320.993,
   319.0289,
   319.0289,
   319.0289,
   0};
   graph = new TGraph(86,PDF_fx9,PDF_fy9);
   graph->SetName("PDF");
   graph->SetTitle("Projection of pdfCTAUMASS_Tot_PbPb");
   graph->SetFillColor(1);
   graph->SetLineWidth(3);
   graph->SetMarkerStyle(20);
   
   TH1F *Graph_PDF9 = new TH1F("Graph_PDF9","Projection of pdfCTAUMASS_Tot_PbPb",100,2.48,3.62);
   Graph_PDF9->SetMinimum(1.137294);
   Graph_PDF9->SetMaximum(1137.294);
   Graph_PDF9->SetDirectory(0);
   Graph_PDF9->SetStats(0);
   Graph_PDF9->SetLineStyle(0);
   Graph_PDF9->SetMarkerStyle(20);
   Graph_PDF9->GetXaxis()->SetLabelFont(42);
   Graph_PDF9->GetXaxis()->SetLabelOffset(0.007);
   Graph_PDF9->GetXaxis()->SetLabelSize(0.05);
   Graph_PDF9->GetXaxis()->SetTitleSize(0.06);
   Graph_PDF9->GetXaxis()->SetTitleOffset(0.9);
   Graph_PDF9->GetXaxis()->SetTitleFont(42);
   Graph_PDF9->GetYaxis()->SetLabelFont(42);
   Graph_PDF9->GetYaxis()->SetLabelOffset(0.007);
   Graph_PDF9->GetYaxis()->SetLabelSize(0.05);
   Graph_PDF9->GetYaxis()->SetTitleSize(0.06);
   Graph_PDF9->GetYaxis()->SetTitleOffset(1.25);
   Graph_PDF9->GetYaxis()->SetTitleFont(42);
   Graph_PDF9->GetZaxis()->SetLabelFont(42);
   Graph_PDF9->GetZaxis()->SetLabelOffset(0.007);
   Graph_PDF9->GetZaxis()->SetLabelSize(0.05);
   Graph_PDF9->GetZaxis()->SetTitleSize(0.06);
   Graph_PDF9->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_PDF9);
   
   graph->Draw("l");
   
   TH1D *frame_25fab720__6 = new TH1D("frame_25fab720__6","",36,2.6,3.5);
   frame_25fab720__6->SetBinContent(1,1091.736);
   frame_25fab720__6->SetMinimum(224.1707);
   frame_25fab720__6->SetMaximum(2746.091);
   frame_25fab720__6->SetEntries(1);
   frame_25fab720__6->SetDirectory(0);
   frame_25fab720__6->SetStats(0);
   frame_25fab720__6->SetLineStyle(0);
   frame_25fab720__6->SetMarkerStyle(20);
   frame_25fab720__6->GetXaxis()->CenterTitle(true);
   frame_25fab720__6->GetXaxis()->SetLabelFont(42);
   frame_25fab720__6->GetXaxis()->SetLabelOffset(3);
   frame_25fab720__6->GetXaxis()->SetLabelSize(0.05);
   frame_25fab720__6->GetXaxis()->SetTitleSize(0.045);
   frame_25fab720__6->GetXaxis()->SetTitleOffset(3);
   frame_25fab720__6->GetXaxis()->SetTitleFont(42);
   frame_25fab720__6->GetYaxis()->SetTitle("Events / ( 0.025 GeV/c^{2} )");
   frame_25fab720__6->GetYaxis()->SetLabelFont(42);
   frame_25fab720__6->GetYaxis()->SetLabelOffset(0.007);
   frame_25fab720__6->GetYaxis()->SetTitleOffset(1.7);
   frame_25fab720__6->GetYaxis()->SetTitleFont(42);
   frame_25fab720__6->GetZaxis()->SetLabelFont(42);
   frame_25fab720__6->GetZaxis()->SetLabelOffset(0.007);
   frame_25fab720__6->GetZaxis()->SetLabelSize(0.05);
   frame_25fab720__6->GetZaxis()->SetTitleSize(0.06);
   frame_25fab720__6->GetZaxis()->SetTitleFont(42);
   frame_25fab720__6->Draw("AXISSAME");
//   TLatex *   tex = new TLatex(0.2,0.735,"N_{bkg}^{PbPb} = 12688#pm89 ");
//tex->SetNDC();
//   tex->SetTextSize(0.026);
//   tex->SetLineWidth(2);
//   tex->Draw();
//      tex = new TLatex(0.2,0.69,"N_{J/#psi}^{PbPb} = 3622#pm59 ");
//tex->SetNDC();
//   tex->SetTextSize(0.026);
//   tex->SetLineWidth(2);
//   tex->Draw();
//      tex = new TLatex(0.2,0.645,"b_{J/#psi}^{PbPb} = 0.2105#pm0.0118");
//tex->SetNDC();
//   tex->SetTextSize(0.026);
//   tex->SetLineWidth(2);
//   tex->Draw();
//      tex = new TLatex(0.2,0.86,"2015 HI Soft Muon ID");
//tex->SetNDC();
//   tex->SetTextSize(0.03);
//   tex->SetLineWidth(2);
//   tex->Draw();
//      tex = new TLatex(0.2,0.815,"HLT_HIL1DoubleMu0_v1");
//tex->SetNDC();
//   tex->SetTextSize(0.03);
//   tex->SetLineWidth(2);
//   tex->Draw();
//      tex = new TLatex(0.5175,0.725,"1.8 < |y^{#mu#mu}| < 2.4");
//tex->SetNDC();
//   tex->SetTextSize(0.03);
//   tex->SetLineWidth(2);
//   tex->Draw();
//      tex = new TLatex(0.5175,0.68,"4.5 < p_{T}^{#mu#mu} < 5.5 GeV/c");
//tex->SetNDC();
//   tex->SetTextSize(0.03);
//   tex->SetLineWidth(2);
//   tex->Draw();
//      tex = new TLatex(0.5175,0.635,"Cent. 0-100%");
//tex->SetNDC();
//   tex->SetTextSize(0.03);
//   tex->SetLineWidth(2);
//   tex->Draw();
  
  TLatex * tex = new TLatex(0.21,0.78,"4.5 #leq p_{T}^{#scale[1.3]{#mu#mu}} < 5.5 GeV/c");
  tex->SetNDC();
  tex->SetTextSize(0.036);
  tex->SetLineWidth(2);
  tex->Draw();
  tex = new TLatex(0.21,0.84,"1.8 #leq |y^{#scale[1.3]{#mu#mu}}| < 2.4");
  tex->SetNDC();
  tex->SetTextSize(0.036);
  tex->SetLineWidth(2);
  tex->Draw();
  tex = new TLatex(0.21,0.72,"Cent. 0-100%");
  tex->SetNDC();
  tex->SetTextSize(0.036);
  tex->SetLineWidth(2);
  tex->Draw();
  
  //   TLegend *leg = new TLegend(0.49,0.70,0.718,0.8809,NULL,"brNDC");
//  TLegend *leg = new TLegend(0.6240602,0.498708,0.8508772,0.6795866,NULL,"brNDC");
//  TLegend *leg = new TLegend(0.2055138,0.4728682,0.4323308,0.6537468,NULL,"brNDC");
  TLegend *leg = new TLegend(0.1766917,0.4702842,0.4035088,0.6511628,NULL,"brNDC");
  leg->SetTextFont(62);
  leg->SetTextSize(0.036);
  leg->SetLineColor(0);
  leg->SetLineStyle(0);
  leg->SetLineWidth(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(1001);
  TLegendEntry *entry=leg->AddEntry("dOS","Data","pe");
  entry->SetLineColor(1);
  entry->SetLineStyle(1);
  entry->SetLineWidth(1);
  entry->SetMarkerColor(1);
  entry->SetMarkerStyle(8);
  entry->SetMarkerSize(1.2);
  entry->SetTextFont(62);
  entry=leg->AddEntry("PDF","Total fit","l");
  
  ci = TColor::GetColor("#0000ff");
  entry->SetLineColor(ci);
  entry->SetLineStyle(1);
  entry->SetLineWidth(3);
  entry->SetMarkerColor(1);
  entry->SetMarkerStyle(21);
  entry->SetMarkerSize(1);
  entry->SetTextFont(62);
  entry=leg->AddEntry("JPSIPR","Prompt J/#psi","l");
  
  ci = TColor::GetColor("#660000");
  entry->SetLineColor(ci);
  entry->SetLineStyle(9);
  entry->SetLineWidth(3);
  entry->SetMarkerColor(1);
  entry->SetMarkerStyle(21);
  entry->SetMarkerSize(1);
  entry->SetTextFont(62);
  entry=leg->AddEntry("JPSINOPR","J/#psi from b hadrons","l");
  
  ci = TColor::GetColor("#006600");
  entry->SetLineColor(ci);
  entry->SetLineStyle(4);
  entry->SetLineWidth(5);
  entry->SetMarkerColor(1);
  entry->SetMarkerStyle(21);
  entry->SetMarkerSize(1);
  entry->SetTextFont(62);
  
  ci = TColor::GetColor("#9966ff");
  entry->SetFillColor(ci);
  entry->SetFillStyle(1001);
  entry->SetLineColor(1);
  entry->SetLineStyle(1);
  entry->SetLineWidth(3);
  entry->SetMarkerColor(1);
  entry->SetMarkerStyle(21);
  entry->SetMarkerSize(1);
  entry->SetTextFont(62);
  entry=leg->AddEntry("BKG","Background","fl");
  
  ci = TColor::GetColor("#99ccff");
  entry->SetFillColor(ci);
  entry->SetFillStyle(1001);
  
  leg->Draw();
  //      tex = new TLatex(0.96,0.9424,"PbPb 351 #mub^{-1} (5.02 TeV)");
  tex = new TLatex(0.96,0.9424,"PbPb 368 #mub^{-1} (5.02 TeV)");
  tex->SetNDC();
  tex->SetTextAlign(31);
  tex->SetTextFont(42);
  tex->SetTextSize(0.04);
  tex->SetLineWidth(2);
  tex->Draw();
  tex = new TLatex(0.16,0.9424,"");
  tex->SetNDC();
  tex->SetTextFont(42);
  tex->SetTextSize(0.04);
  tex->SetLineWidth(2);
  tex->Draw();
  tex = new TLatex(0.914,0.869275,"CMS");
  tex->SetNDC();
  tex->SetTextAlign(33);
  tex->SetTextFont(61);
  tex->SetTextSize(0.06);
  tex->SetLineWidth(2);
  tex->Draw();
//  tex = new TLatex(0.914,0.797275,"Preliminary");
//  tex->SetNDC();
//  tex->SetTextAlign(33);
//  tex->SetTextFont(52);
//  tex->SetTextSize(0.0456);
//  tex->SetLineWidth(2);
//  tex->Draw();
  //   pad1_PbPb->Modified();
  cMassFig_PbPb->cd();
  
//  \\\\\//////
//  
//  TLatex * tex = new TLatex(0.21,0.77,"4.5 #leq p_{T}^{#mu#mu} < 5.5 GeV/c");
//  tex->SetNDC();
//  tex->SetTextSize(0.03);
//  tex->SetLineWidth(2);
//  tex->Draw();
//  tex = new TLatex(0.21,0.86,"1.8 #leq |y^{#mu#mu}| < 2.4");
//  tex->SetNDC();
//  tex->SetTextSize(0.03);
//  tex->SetLineWidth(2);
//  tex->Draw();
//  tex = new TLatex(0.21,0.725,"Cent. 0-100%");
//  tex->SetNDC();
//  tex->SetTextSize(0.03);
//  tex->SetLineWidth(2);
//  tex->Draw();
//  
//   TLegend *leg = new TLegend(0.5175,0.7602,0.718,0.8809,NULL,"brNDC");
//   leg->SetTextFont(62);
//   leg->SetTextSize(0.03);
//   leg->SetLineColor(1);
//   leg->SetLineStyle(1);
//   leg->SetLineWidth(1);
//   leg->SetFillColor(0);
//   leg->SetFillStyle(1001);
//   TLegendEntry *entry=leg->AddEntry("dOS","Data","pe");
//   entry->SetLineColor(1);
//   entry->SetLineStyle(1);
//   entry->SetLineWidth(1);
//   entry->SetMarkerColor(1);
//   entry->SetMarkerStyle(8);
//   entry->SetMarkerSize(1.2);
//   entry->SetTextFont(62);
//   entry=leg->AddEntry("PDF","Total fit","l");
//   entry->SetLineColor(1);
//   entry->SetLineStyle(1);
//   entry->SetLineWidth(3);
//   entry->SetMarkerColor(1);
//   entry->SetMarkerStyle(21);
//   entry->SetMarkerSize(1);
//   entry->SetTextFont(62);
//   entry=leg->AddEntry("JPSIPR","Prompt J/#psi","l");
//
//   ci = TColor::GetColor("#660000");
//   entry->SetLineColor(ci);
//   entry->SetLineStyle(1);
//   entry->SetLineWidth(3);
//   entry->SetMarkerColor(1);
//   entry->SetMarkerStyle(21);
//   entry->SetMarkerSize(1);
//   entry->SetTextFont(62);
//   entry=leg->AddEntry("JPSINOPR","Non-Prompt J/#psi","l");
//
//   ci = TColor::GetColor("#006600");
//   entry->SetLineColor(ci);
//   entry->SetLineStyle(1);
//   entry->SetLineWidth(3);
//   entry->SetMarkerColor(1);
//   entry->SetMarkerStyle(21);
//   entry->SetMarkerSize(1);
//   entry->SetTextFont(62);
//   entry=leg->AddEntry("BKG","Background","fl");
//
//   ci = TColor::GetColor("#99ccff");
//   entry->SetFillColor(ci);
//   entry->SetFillStyle(1001);
//
//   ci = TColor::GetColor("#0000ff");
//   entry->SetLineColor(ci);
//   entry->SetLineStyle(2);
//   entry->SetLineWidth(3);
//   entry->SetMarkerColor(1);
//   entry->SetMarkerStyle(21);
//   entry->SetMarkerSize(1);
//   entry->SetTextFont(62);
//   leg->Draw();

   cMassFig_PbPb->cd();
   cMassFig_PbPb->Modified();
   cMassFig_PbPb->cd();
   cMassFig_PbPb->SetSelected(cMassFig_PbPb);
}
