void read_tfile()
{
  // ------------------------------------------------------------------
  // read inputs
  // ------------------------------------------------------------------
  const char *DATA_INPUT_PATH = "/data/users/pjgwak/work/daily_code_tracker/2026/02/02_2d_fit_OO/onia_to_skim/roodataset_roots/RooDataSet_miniAOD_isMC0_Jpsi_cent0_200_Effw0_Accw0_PtW0_TnP0_OO25_HLT_OxyL1SingleMuOpen_v1.root";

  auto *data = TFile::Open(DATA_INPUT_PATH);

  data->Print("V");
}