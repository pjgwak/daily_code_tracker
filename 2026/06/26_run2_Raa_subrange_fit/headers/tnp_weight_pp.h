#ifndef tnp_weight_h
#define tnp_weight_h

#include "TMath.h"
#include <tuple>

// IN THIS FILE YOU WILL FIND:
// ++++++++++++++

//BINNED SCALE FACTORS for pp 2017

//To extract the scale factor for a muon of PT=pt and ETA=eta, do:

//#include "tnp_weight.h"
//void testtnpweight(){
// auto SFwithError = tnp_weight_theEfficiencyYouWant_pp( pt, eta );
// double central_value = std::get<0>(SFwithError);
// double stat_error = std::get<1>(SFwithError);
// double syst_error = std::get<2>(SFwithError);
// double tot_error = std::get<3>(SFwithError);
//}

// ++++++++++++++

///////////////////////////////////////////////////
//                 G l o b a l    p p               //
///////////////////////////////////////////////////

//To be used for global muon efficiency, with tight muon kinematic acceptance cuts
//Central scale-factor value, statistical error, systematic error, total error

std::tuple<double,double,double,double> tnp_weight_GlobalMuon_TightAcceptance_pp(double pt, double eta){

  if( (fabs(eta) > 0) && (fabs(eta) < 1.2) ){
    if (pt<3.7) return std::make_tuple(1.05682, 0.0065723, 0.00108251, 0.00666085);
    else if (pt<3.9) return std::make_tuple(1.01678, 0.00586556, 0.000211988, 0.00586939);
    else if (pt<4.2) return std::make_tuple(1.00333, 0.00448914, 0.000743086, 0.00455023);
    else if (pt<4.6) return std::make_tuple(0.992136, 0.00379769, 0.0014736, 0.00407357);
    else if (pt<5.2) return std::make_tuple(0.99067, 0.00315441, 9.585e-05, 0.00315587);
    else if (pt<7) return std::make_tuple(0.996967, 0.00215854, 0.00149636, 0.00262647);
    else if (pt<10.5) return std::make_tuple(0.991201, 0.0024164, 0.000896526, 0.00257735);
    else if (pt<30) return std::make_tuple(0.980219, 0.0037124, 0.00096564, 0.00383593);
  }

  if( (fabs(eta) > 1.2) && (fabs(eta) < 1.8) ){
    if (pt<2.4) return std::make_tuple(1.00302, 0.0248043, 0.0104505, 0.0269159);
    else if (pt<2.8) return std::make_tuple(1.00085, 0.0107777, 0.003583, 0.0113577);
    else if (pt<3.2) return std::make_tuple(1.00227, 0.00740985, 0.00125896, 0.00751604);
    else if (pt<3.9) return std::make_tuple(0.994838, 0.00489704, 0.00233562, 0.00542551);
    else if (pt<5) return std::make_tuple(0.98849, 0.00412823, 0.00213652, 0.00464833);
    else if (pt<8) return std::make_tuple(0.996084, 0.00354905, 0.00058147, 0.00359637);
    else if (pt<30) return std::make_tuple(0.985733, 0.00521499, 0.000619272, 0.00525163);
  }

  if( (fabs(eta) > 1.8) && (fabs(eta) < 2.4) ){
    if (pt<1.8) return std::make_tuple(0.95216, 0.0192699, 0.0103348, 0.0218664);
    else if (pt<2.3) return std::make_tuple(0.976063, 0.0106322, 0.00791429, 0.0132544);
    else if (pt<2.9) return std::make_tuple(0.977207, 0.00820585, 0.00155826, 0.0083525);
    else if (pt<3.7) return std::make_tuple(0.978061, 0.00685112, 0.00197487, 0.00713008);
    else if (pt<6) return std::make_tuple(0.980203, 0.00476456, 0.000700885, 0.00481583);
    else if (pt<30) return std::make_tuple(0.967581, 0.00554748, 0.00538024, 0.00772797);
  }

  return std::make_tuple(1,0,0,0);
}

//To be used for HybridSoftID efficiency, with loose muon kinematic acceptance cuts
//Central scale-factor value, statistical error, systematic error, total error

std::tuple<double,double,double,double> tnp_weight_HybridSoftID_LooseAcceptance_pp(double pt, double eta){

  if( (fabs(eta) > 0) && (fabs(eta) < 1.2) ){
    if (pt<3.7) return std::make_tuple(0.991509, 0.00132469, 0.000169579, 0.0013355);
    else if (pt<4) return std::make_tuple(0.994812, 0.00098049, 9.3781e-05, 0.000984965);
    else if (pt<4.5) return std::make_tuple(0.996624, 0.000652936, 5.94331e-05, 0.000655635);
    else if (pt<5) return std::make_tuple(0.997641, 0.000608668, 0.000281891, 0.000670775);
    else if (pt<5.5) return std::make_tuple(0.998693, 0.000581153, 4.7311e-05, 0.000583075);
    else if (pt<6) return std::make_tuple(0.997692, 0.000592181, 4.28475e-05, 0.000593729);
    else if (pt<7) return std::make_tuple(0.99823, 0.000488995, 0.000133913, 0.000507);
    else if (pt<8) return std::make_tuple(0.998932, 0.000570471, 6.84542e-05, 0.000574564);
    else if (pt<10.5) return std::make_tuple(0.999942, 0.000461622, 0.000246741, 0.000523427);
    else if (pt<14) return std::make_tuple(1.00025, 0.000588468, 0.000412485, 0.000718636);
    else if (pt<18) return std::make_tuple(0.998795, 0.00103439, 0.000343453, 0.00108992);
    else if (pt<30) return std::make_tuple(1.00013, 0.00119241, 0.000258227, 0.00122005);
  }

  if( (fabs(eta) > 1.2) && (fabs(eta) < 1.8) ){
    if (pt<2.5) return std::make_tuple(0.999008, 0.00101804, 0.00022889, 0.00104346);
    else if (pt<3) return std::make_tuple(0.997385, 0.000670367, 0.000166048, 0.000690626);
    else if (pt<3.5) return std::make_tuple(0.99795, 0.000566551, 7.41528e-05, 0.000571383);
    else if (pt<4) return std::make_tuple(0.997328, 0.000565301, 0.000106841, 0.000575308);
    else if (pt<4.5) return std::make_tuple(0.998049, 0.000586026, 0.000166683, 0.00060927);
    else if (pt<5) return std::make_tuple(0.99893, 0.000624649, 0.000139086, 0.000639946);
    else if (pt<6) return std::make_tuple(0.99775, 0.000524755, 0.000129245, 0.000540437);
    else if (pt<7) return std::make_tuple(0.997207, 0.000631782, 0.000158721, 0.000651415);
    else if (pt<9) return std::make_tuple(0.997041, 0.000592354, 0.000262767, 0.00064802);
    else if (pt<14) return std::make_tuple(0.997027, 0.000643335, 0.000101755, 0.000651333);
    else if (pt<18) return std::make_tuple(0.998987, 0.00145434, 0.000105944, 0.00145819);
    else if (pt<30) return std::make_tuple(0.998832, 0.00176589, 0.000185321, 0.00177559);
  }

  if( (fabs(eta) > 1.8) && (fabs(eta) < 2.1) ){
    if (pt<1.75) return std::make_tuple(0.99501, 0.00190285, 0.000614269, 0.00199954);
    else if (pt<2.1) return std::make_tuple(0.995534, 0.000958298, 0.000938682, 0.00134144);
    else if (pt<2.5) return std::make_tuple(0.998697, 0.000606511, 0.000452902, 0.000756951);
    else if (pt<3) return std::make_tuple(0.998797, 0.000516932, 0.000241382, 0.000570512);
    else if (pt<3.5) return std::make_tuple(0.999229, 0.000521865, 1.82422e-05, 0.000522184);
    else if (pt<4) return std::make_tuple(0.998185, 0.000596435, 0.000199466, 0.000628905);
    else if (pt<4.5) return std::make_tuple(0.997687, 0.000654825, 0.000242263, 0.000698203);
    else if (pt<5.25) return std::make_tuple(0.998599, 0.000576324, 0.000257014, 0.000631036);
    else if (pt<6.5) return std::make_tuple(0.999244, 0.000485306, 0.000321732, 0.000582266);
    else if (pt<9) return std::make_tuple(0.998085, 0.000535479, 0.000314486, 0.000620998);
    else if (pt<13) return std::make_tuple(0.998346, 0.000716329, 0.000608368, 0.000939807);
    else if (pt<30) return std::make_tuple(0.998398, 0.00112506, 0.000231803, 0.00114869);
  }

  if( (fabs(eta) > 2.1) && (fabs(eta) < 2.4) ){
    if (pt<1.9) return std::make_tuple(0.99387, 0.00114168, 0.000306098, 0.001182);
    else if (pt<2.5) return std::make_tuple(0.996848, 0.00781191, 0.000367591, 0.00782055);
    else if (pt<3) return std::make_tuple(0.998103, 0.000651977, 0.000306688, 0.000720508);
    else if (pt<3.5) return std::make_tuple(1.00042, 0.000618082, 0.000608206, 0.000867145);
    else if (pt<4) return std::make_tuple(0.998227, 0.000779589, 0.000550001, 0.000954076);
    else if (pt<4.75) return std::make_tuple(0.999744, 0.000565002, 0.000220046, 0.000606339);
    else if (pt<5.5) return std::make_tuple(0.998818, 0.000747367, 0.000685803, 0.00101434);
    else if (pt<8) return std::make_tuple(0.997726, 0.000576727, 0.000259674, 0.000632491);
    else if (pt<12) return std::make_tuple(0.999967, 0.000812264, 0.000550096, 0.000981008);
    else if (pt<30) return std::make_tuple(1.00071, 1.49564e-12, 0.000545007, 0.000545007);
  }
  return std::make_tuple(1,0,0,0);

}


//To be used for HybridSoftID+(J/psi trigger) efficiency, with tight muon kinematic acceptance cuts (the loose cuts cannot be used with trigger efficiency)
//Central scale-factor value, statistical error, systematic error, total error

std::tuple<double,double,double,double> tnp_weight_HybridSoftIDTrigger_TightAcceptance_pp(double pt, double eta){
  if( (fabs(eta) > 0) && (fabs(eta) < 1.2) ){
    if (pt<3.75) return std::make_tuple(1.04685, 0.00307983, 0.000900168, 0.00320868);
    else if (pt<4) return std::make_tuple(1.03785, 0.00231993, 0.000865001, 0.00247594);
    else if (pt<4.25) return std::make_tuple(1.02122, 0.0019751, 0.000129247, 0.00197932);
    else if (pt<4.5) return std::make_tuple(1.01242, 0.00183717, 0.000452042, 0.00189196);
    else if (pt<4.75) return std::make_tuple(1.01333, 0.00178571, 0.000501881, 0.0018549);
    else if (pt<5) return std::make_tuple(1.01311, 0.00180214, 0.00038986, 0.00184382);
    else if (pt<5.25) return std::make_tuple(1.01071, 0.00179057, 0.000558043, 0.00187552);
    else if (pt<5.5) return std::make_tuple(1.01492, 0.00178375, 0.000621342, 0.00188887);
    else if (pt<5.75) return std::make_tuple(1.01162, 0.00190487, 0.000708795, 0.00203246);
    else if (pt<6) return std::make_tuple(1.00824, 0.00201297, 0.00144451, 0.00247763);
    else if (pt<6.5) return std::make_tuple(0.99932, 0.00150513, 0.000387211, 0.00155414);
    else if (pt<7) return std::make_tuple(0.990015, 0.00168554, 0.000334128, 0.00171833);
    else if (pt<7.5) return std::make_tuple(0.989662, 0.00184611, 0.000377601, 0.00188433);
    else if (pt<8) return std::make_tuple(0.991271, 0.00199633, 0.000397838, 0.00203559);
    else if (pt<9.25) return std::make_tuple(0.993863, 0.00145115, 0.000219424, 0.00146765);
    else if (pt<10.5) return std::make_tuple(0.98847, 0.00181365, 0.00046839, 0.00187316);
    else if (pt<12.25) return std::make_tuple(0.988261, 0.0020079, 0.000342797, 0.00203695);
    else if (pt<14) return std::make_tuple(0.986502, 0.00269587, 0.000312191, 0.00271388);
    else if (pt<16) return std::make_tuple(0.967524, 0.00360689, 0.000816031, 0.00369805);
    else if (pt<18) return std::make_tuple(0.974107, 0.0047847, 0.00173171, 0.00508843);
    else if (pt<24) return std::make_tuple(0.971851, 0.00434543, 0.000750224, 0.00440971);
    else if (pt<30) return std::make_tuple(0.95726, 0.00782987, 0.00260996, 0.0082534);
  }

  if( (fabs(eta) > 1.2) && (fabs(eta) < 1.8) ){
    if (pt<2.5) return std::make_tuple(1.05392, 0.00741907, 0.0086659, 0.0114079);
    else if (pt<2.75) return std::make_tuple(1.05418, 0.00510863, 0.000493877, 0.00513245);
    else if (pt<3) return std::make_tuple(1.03165, 0.00362626, 0.00051552, 0.00366272);
    else if (pt<3.25) return std::make_tuple(1.00137, 0.00324096, 0.000650532, 0.00330561);
    else if (pt<3.5) return std::make_tuple(0.985531, 0.00325619, 0.00101988, 0.00341217);
    else if (pt<3.75) return std::make_tuple(0.979026, 0.00333719, 0.00080565, 0.00343306);
    else if (pt<4) return std::make_tuple(0.978533, 0.00348675, 0.00128796, 0.00371702);
    else if (pt<4.5) return std::make_tuple(0.981693, 0.00262222, 0.00138638, 0.00296615);
    else if (pt<5) return std::make_tuple(0.988393, 0.00282795, 0.000574534, 0.00288572);
    else if (pt<5.5) return std::make_tuple(0.988535, 0.00296014, 0.0015193, 0.00332727);
    else if (pt<6) return std::make_tuple(0.995795, 0.003052, 0.000849741, 0.00316809);
    else if (pt<7) return std::make_tuple(0.986439, 0.00239136, 0.00132843, 0.00273557);
    else if (pt<8.5) return std::make_tuple(0.991989, 0.0023577, 0.00141202, 0.00274819);
    else if (pt<10) return std::make_tuple(0.993719, 0.00296632, 0.00125018, 0.003219);
    else if (pt<12) return std::make_tuple(0.990475, 0.0031633, 0.00082619, 0.00326941);
    else if (pt<14.5) return std::make_tuple(0.99984, 0.00391145, 0.00226697, 0.00452091);
    else if (pt<19) return std::make_tuple(0.999167, 0.00490761, 0.00146068, 0.00512038);
    else if (pt<30) return std::make_tuple(0.998484, 0.00721995, 0.000958602, 0.00728331);
  }

  if( (fabs(eta) > 1.8) && (fabs(eta) < 2.1) ){
    if (pt<1.9) return std::make_tuple(1.16258, 0.00663765, 0.00518218, 0.008421);
    else if (pt<2.2) return std::make_tuple(1.08194, 0.00434659, 0.00165702, 0.00465173);
    else if (pt<2.5) return std::make_tuple(1.04721, 0.00279163, 0.00198365, 0.00342462);
    else if (pt<2.8) return std::make_tuple(1.0117, 0.00250143, 0.00103727, 0.00270797);
    else if (pt<3.1) return std::make_tuple(1.0045, 0.00241179, 0.000301947, 0.00243062);
    else if (pt<3.4) return std::make_tuple(0.998153, 0.00255819, 0.00107157, 0.00277355);
    else if (pt<4) return std::make_tuple(0.991499, 0.00205979, 0.000162552, 0.00206619);
    else if (pt<5) return std::make_tuple(0.982154, 0.00203163, 0.00109106, 0.00230607);
    else if (pt<6) return std::make_tuple(0.979281, 0.00269371, 0.00139839, 0.00303506);
    else if (pt<8) return std::make_tuple(0.978722, 0.00284354, 0.00105435, 0.00303272);
    else if (pt<11) return std::make_tuple(0.982305, 0.00415391, 0.000812637, 0.00423265);
    else if (pt<15) return std::make_tuple(0.98089, 0.00737357, 0.000877749, 0.00742563);
    else if (pt<30) return std::make_tuple(0.985428, 0.0132041, 0.012101, 0.0179104);
  }

  if( (fabs(eta) > 2.1) && (fabs(eta) < 2.4) ){
    if (pt<1.9) return std::make_tuple(1.10953, 0.00267574, 0.000888127, 0.00281928);
    else if (pt<2.3) return std::make_tuple(1.04759, 0.00206322, 0.00049553, 0.0021219);
    else if (pt<2.7) return std::make_tuple(1.02649, 0.00198497, 0.00156344, 0.00252675);
    else if (pt<3) return std::make_tuple(1.01133, 0.00249575, 0.000958998, 0.00267366);
    else if (pt<3.5) return std::make_tuple(1.01379, 0.00208413, 0.00104966, 0.00233354);
    else if (pt<4.25) return std::make_tuple(1.00341, 0.00202201, 0.000761361, 0.0021606);
    else if (pt<5) return std::make_tuple(0.997733, 0.00254249, 0.00299566, 0.00392915);
    else if (pt<6) return std::make_tuple(0.991922, 0.00275362, 0.000262191, 0.00276608);
    else if (pt<7) return std::make_tuple(0.9968, 0.00362768, 0.00239638, 0.00434772);
    else if (pt<8.5) return std::make_tuple(1.00011, 0.00361481, 0.00168078, 0.00398646);
    else if (pt<13) return std::make_tuple(1.00663, 0.00411005, 0.00220707, 0.00466516);
    else if (pt<30) return std::make_tuple(0.992803, 0.010354, 0.0109342, 0.0150586);
  }
  return std::make_tuple(1,0,0,0);

}


#endif
