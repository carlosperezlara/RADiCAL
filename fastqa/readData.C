
int readData(int run=0, Long64_t totalevents=100000, float UST_thr=1000, float threshold=270) {
  gStyle->SetOptStat(0);
  TFile *file = new TFile(  Form("/Volumes/uva/testbeam_2021_06_data/merged/Run_%d.root",run) );
  TTree *tree = (TTree*) file->Get("pulse");

  int event;
  ushort tc[2]; // trigger counter bin
  float times[2][1024]; // calibrated time
  float channel[18][1024]; // calibrated input (in V)
  tree->SetBranchAddress("event", &event);
  tree->SetBranchAddress("tc", tc);
  tree->SetBranchAddress("channel", channel);
  tree->SetBranchAddress("times", times);
  
  const int nchn = 8;
  int chnidx[ nchn ] =     {    8,    2,    1,    3,    4,    5,    6,    9};
  TString chnlab[ nchn ] = {"MCP","TRT","UST","TRE","BRE","TLE","BLE","USE"};

  TH2D *summary[nchn];
  for(int i=0; i!=nchn; ++i) {
    summary[i] = new TH2D( Form("STACK_%s", chnlab[i].Data() ),
			   Form("STACK_%s;ns;mV", chnlab[i].Data() ),
			   100,0.0,200.0, 100,-500.1,+500.1 );
  }
  
  Long64_t nentries = TMath::Min( totalevents, tree->GetEntries() );
  nentries = TMath::Max( nentries, Long64_t(0) );
  TGraph *trace[18];
  Long64_t nev=0;
  for(;nev!=nentries;++nev) {
    tree->GetEntry( nev );
    if(nev%500==0)
      cout << "Events read so far: " << nev << endl;
    for(int i=0; i!=18; ++i) {
      int group = i/9;
      trace[i] = new TGraph( 1024, times[group], channel[i] );
    }
    // ===================================
    // Access data through TRACE ( TGraph )

    for(int idx=0; idx!=nchn; ++idx) {
      for(int sample=0; sample!=1024; ++sample) {
	double x = trace[ chnidx[ idx ] ]->GetPointX( sample );
	double y = trace[ chnidx[ idx ] ]->GetPointY( sample );
	summary[ idx ]->Fill( x, y );
      }
    }
    
    // ===================================
    
    for(int i=0; i!=18; ++i) delete trace[i];
  }
  cout << "total number of events read: " << nev << endl;

  TCanvas *c1 = new TCanvas("c1",Form("RUN %d",run));
  c1->Divide(4,2);
  for(int idx=0; idx!=nchn; ++idx) {
    c1->cd(idx+1)->SetLogz(1);
    summary[idx]->Draw("colz");
  }
  
  return 0;
}
