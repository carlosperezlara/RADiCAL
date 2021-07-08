#include "WaveForm.cc"
#include "Reader.cc"
#include "DT5742DATReader.cc"

void FormatSummary(TH2D *obj) {
  obj->SetTitleSize(0.5);
  obj->GetYaxis()->SetLabelSize(0.06);
  obj->GetXaxis()->SetLabelSize(0.06);
}

void FormatSummary(TH1 *obj) {
  obj->SetTitleSize(0.5);
  obj->GetYaxis()->SetLabelSize(0.06);
  obj->GetXaxis()->SetLabelSize(0.06);
}

void FormatStats(TH1 *obj) {
  obj->SetTitleSize(0.5);
  obj->GetYaxis()->SetLabelSize(0.06);
  obj->GetXaxis()->SetLabelSize(0.13);
  obj->SetFillColor(kGreen-3);
}

int reader(int run, int totalevents=5000, float threshold=270) {
  gStyle->SetOptStat(0);
  DT5742DATReader file(  Form("/data/testbeam_2021_06_data/DRS_data/Run_%d.dat",run) );

  file.ReadHeader();

  WaveForm *trace[18];
  TH2D *summary[18];
  TH2D *corrAmp[18];
  TH1D *amplitudes[18];
  TH1D *times[18];
  TH1D *event7[18];
  TProfile *waveprof[18];
  TH1D *effs[18];
  
  for(int i=0; i!=18; ++i) {
    summary[i] = file.GetSummaryPlot(i);
    trace[i]  = file.GetTrace(i);
    trace[i]->CreateProfile();
    waveprof[i] = trace[i]->GetProfile();
    corrAmp[i] = new TH2D( Form("A%d",i), Form("A%d",i), 100, 0, -500, 100, 0, -500);
    amplitudes[i] = new TH1D( Form("AA%d",i), Form("AA%d",i), 1000, -500, +500 );
    times[i] = new TH1D( Form("TT%d",i), Form("TT%d",i), 1000, -0.5, 1023.5 );
    effs[i] = new TH1D( Form("EFFS%d",i), Form("EFFS%d",i), 4, -0.5, 3.5 );
    effs[i]->GetXaxis()->SetBinLabel(1,Form(">%.0f",threshold));
    effs[i]->GetXaxis()->SetBinLabel(2,Form("<%.0f",threshold));
    effs[i]->GetXaxis()->SetBinLabel(3,"<-300");
    FormatSummary(summary[i]);
    FormatSummary(waveprof[i]);
    FormatStats(effs[i]);
  }
  int nev=0;
  //for(;nev!=5000;++nev) {
  for(;nev!=totalevents;++nev) {
    //for(int nev=0;;++nev) {
    if(!file.ReadEvent()) break;
    if(nev%500==0)
      cout << "Events read so far: " << nev << endl;
    // i have access through trace
    double ampl[18];
    int time[18];
    for(int i=0; i!=18; ++i) {
      trace[i]->ComputeMin(1,1024,ampl[i],time[i]);
      amplitudes[i]->Fill( ampl[i] );
      if(ampl[i]<threshold) {
	trace[i]->FillProfile();
	times[i]->Fill( time[i] );
	effs[i]->Fill(1);
	if(ampl[i]<-300)
	  effs[i]->Fill(2);
      } else {
	effs[i]->Fill(0);
      }
    }
  }
  cout << "total number of events read: " << nev << endl;

  TCanvas *main_DS = new TCanvas("main_DS","Downstream",0,0,800,600);
  main_DS->Divide(3,3);
  summary[5]->SetTitle("TLE"); // TLE
  summary[3]->SetTitle("TRE"); // TRE
  summary[6]->SetTitle("BLE"); // BLE
  summary[4]->SetTitle("BRE"); // BRE
  summary[2]->SetTitle("TRT"); // TRT
  main_DS->cd(1)->SetLogz(1); summary[5]->DrawCopy("colz");
  main_DS->cd(3)->SetLogz(1); summary[3]->DrawCopy("colz");
  main_DS->cd(7)->SetLogz(1); summary[6]->DrawCopy("colz");
  main_DS->cd(9)->SetLogz(1); summary[4]->DrawCopy("colz");
  main_DS->cd(2)->SetLogz(1); summary[2]->DrawCopy("colz");
  TCanvas *main_US = new TCanvas("main_US","Upstream",800,0,800,600);
  main_US->Divide(3,3);
  summary[1]->SetTitle("UST"); // TLE
  summary[8]->SetTitle("USE"); // TLE
  summary[16]->SetTitle("MCP"); // TLE
  main_US->cd(2)->SetLogz(1); summary[1]->DrawCopy("colz"); // UST
  main_US->cd(3)->SetLogz(1); summary[8]->DrawCopy("colz"); // USE
  main_US->cd(5)->SetLogz(1); summary[16]->DrawCopy("colz"); // MCP  
  
  //================== 
  
  TCanvas *main2_DS = new TCanvas("main2_DS","Downstream",0,600,800,600);
  main2_DS->Divide(3,3);
  waveprof[5]->SetTitle("TLE"); // TLE
  waveprof[3]->SetTitle("TRE"); // TRE
  waveprof[6]->SetTitle("BLE"); // BLE
  waveprof[4]->SetTitle("BRE"); // BRE
  waveprof[2]->SetTitle("TRT"); // TRT
  main2_DS->cd(1); waveprof[5]->DrawCopy(); // TLE
  main2_DS->cd(3); waveprof[3]->DrawCopy(); // TRE
  main2_DS->cd(7); waveprof[6]->DrawCopy(); // BLE
  main2_DS->cd(9); waveprof[4]->DrawCopy(); // BRE
  main2_DS->cd(2); waveprof[2]->DrawCopy(); // TRT
  TCanvas *main2_US = new TCanvas("main2_US","Upstream",800,600,800,600);
  main2_US->Divide(3,3);
  waveprof[1]->SetTitle("UST"); // TLE
  waveprof[8]->SetTitle("USE"); // TLE
  waveprof[16]->SetTitle("MCP"); // TLE
  main2_US->cd(2); waveprof[1]->DrawCopy(); // UST
  main2_US->cd(3); waveprof[8]->DrawCopy(); // USE
  main2_US->cd(5); waveprof[16]->DrawCopy(); // MCP  

  //================== 
  
  TCanvas *main3_DS = new TCanvas("main3_DS","Downstream",0,2000,800,600);
  main3_DS->Divide(3,3);
  amplitudes[5]->SetTitle("TLE"); // TLE
  amplitudes[3]->SetTitle("TRE"); // TRE
  amplitudes[6]->SetTitle("BLE"); // BLE
  amplitudes[4]->SetTitle("BRE"); // BRE
  amplitudes[2]->SetTitle("TRT"); // TRT
  main3_DS->cd(1)->SetLogy(1); amplitudes[5]->DrawCopy(); // TLE
  main3_DS->cd(3)->SetLogy(1); amplitudes[3]->DrawCopy(); // TRE
  main3_DS->cd(7)->SetLogy(1); amplitudes[6]->DrawCopy(); // BLE
  main3_DS->cd(9)->SetLogy(1); amplitudes[4]->DrawCopy(); // BRE
  main3_DS->cd(2)->SetLogy(1); amplitudes[2]->DrawCopy(); // TRT
  TCanvas *main3_US = new TCanvas("main3_US","Upstream",800,2000,800,600);
  main3_US->Divide(3,3);
  amplitudes[1]->SetTitle("UST"); // TLE
  amplitudes[8]->SetTitle("USE"); // TLE
  amplitudes[16]->SetTitle("MCP"); // TLE
  main3_US->cd(2)->SetLogy(1); amplitudes[1]->DrawCopy(); // UST
  main3_US->cd(3)->SetLogy(1); amplitudes[8]->DrawCopy(); // USE
  main3_US->cd(5)->SetLogy(1); amplitudes[16]->DrawCopy(); // MCP  

  //================== 

  TCanvas *main4_DS = new TCanvas("main4_DS","Downstream",0,2000,800,600);
  main4_DS->Divide(3,3);
  for(int i=0; i!=18; ++i) {
    double integ = effs[i]->Integral(1,2);
    effs[i]->Scale(1/integ);
  }
  effs[5]->SetTitle("TLE"); // TLE
  effs[3]->SetTitle("TRE"); // TRE
  effs[6]->SetTitle("BLE"); // BLE
  effs[4]->SetTitle("BRE"); // BRE
  effs[2]->SetTitle("TRT"); // TRT
  main4_DS->cd(1); effs[5]->DrawCopy("hist"); // TLE
  main4_DS->cd(3); effs[3]->DrawCopy("hist"); // TRE
  main4_DS->cd(7); effs[6]->DrawCopy("hist"); // BLE
  main4_DS->cd(9); effs[4]->DrawCopy("hist"); // BRE
  main4_DS->cd(2); effs[2]->DrawCopy("hist"); // TRT
  TCanvas *main4_US = new TCanvas("main4_US","Upstream",800,2000,800,600);
  main4_US->Divide(3,3);
  effs[1]->SetTitle("UST"); // TLE
  effs[8]->SetTitle("USE"); // TLE
  effs[16]->SetTitle("MCP"); // TLE
  main4_US->cd(2); effs[1]->DrawCopy("hist"); // UST
  main4_US->cd(3); effs[8]->DrawCopy("hist"); // USE
  main4_US->cd(5); effs[16]->DrawCopy("hist"); // MCP  
  

  return 0;
}
