void dumptree(const char* filename) {
   TFile *f = TFile::Open(filename);
   if (!f) return;
   TTree *tr = (TTree*) f->Get("fitresults");
   if (!tr) return;
   tr->SetScanField(0);
   tr->Scan("*","","colsize=30");
}
