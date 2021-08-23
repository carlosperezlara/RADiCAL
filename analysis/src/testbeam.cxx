#include "testbeam.h"
#include <TFile.h>
#include <TList.h>
#include <TString.h>
#include <TTree.h>
#include <TCanvas.h>
#include <iostream>
#include <TMath.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TProfile.h>

#include "mcp.h"
#include "energy.h"
#include "timing.h"

//=================================
testbeam::testbeam() {
  fList = NULL;
  fTree = NULL;
  fMCP0 = NULL;
  fEventDisplayMCP=NULL;
  fEvents = NULL;
  for(int i=0; i!=2; ++i) {
    fEventDisplayEnergy[i]=NULL;
  }
}
//=================================
testbeam::~testbeam() {
  if(fList) delete fList;
}
//=================================
void testbeam::SetInputFileName(TString val) {
  TFile *inputfile = new TFile( val.Data() );
  fTree = (TTree*) inputfile->Get("pulse");
}
//=================================
void testbeam::SetOutputFileName(TString val) {
  fOutputFileName = val;
}
//=================================
TGraph* testbeam::DisplayWave(DRSWAVE *A, DRSWAVE *B, DRSWAVE *C, int showfit) {
  TGraph *a, *b, *c;
  double tmp;
  double aplot = 999;
  if(A!=NULL) {
    a = new TGraph(1024,A->GetX(),A->GetY());
    a->SetLineColor(kBlue-3);
    a->SetMarkerColor(kBlue-3);
    a->SetMarkerStyle(24);
    A->GetMinimum(aplot, tmp, 0, 1024, 0);
  }
  double bplot = 999;
  if(B!=NULL) {
    b = new TGraph(1024,B->GetX(),B->GetY());
    b->SetLineColor(kGreen-3);
    b->SetMarkerColor(kGreen-3);
    b->SetMarkerStyle(24);
    B->GetMinimum(bplot, tmp, 0, 1024, 0);
  }
  double cplot = 999;
  if(C!=NULL) {
    c = new TGraph(1024,C->GetX(),C->GetY());
    c->SetLineColor(kRed-3);
    c->SetMarkerColor(kRed-3);
    c->SetMarkerStyle(24);
    C->GetMinimum(cplot, tmp, 0, 1024, 0);
  }

  TGraph *ret = NULL;
  if(aplot>0 && bplot>0 && cplot>0) return NULL;
  if(aplot<bplot&&aplot<cplot)
    ret = a;
  if(bplot<aplot&&bplot<cplot)
    ret = b;
  if(cplot<aplot&&cplot<bplot)
    ret = c;
  if(ret==NULL) return NULL;
  ret->Draw("AP");
  if(A!=NULL) a->Draw("PSAME");
  if(B!=NULL) b->Draw("PSAME");
  if(C!=NULL) c->Draw("PSAME");
  if(showfit==2) {
    if(A!=NULL) {
      fQuadratic->SetParameter(0,A->GetQFa());
      fQuadratic->SetParameter(1,A->GetQFb());
      fQuadratic->SetParameter(2,A->GetQFc());
      fQuadratic->SetRange( A->GetTime()-4, A->GetTime()+4 );
      fQuadratic->SetLineColor(kBlack);
      fQuadratic->DrawCopy("same");
    }
    if(B!=NULL) {
      fQuadratic->SetParameter(0,B->GetQFa());
      fQuadratic->SetParameter(1,B->GetQFb());
      fQuadratic->SetParameter(2,B->GetQFc());
      fQuadratic->SetRange( B->GetTime()-4, B->GetTime()+4 );
      fQuadratic->SetLineColor(kBlack);
      fQuadratic->DrawCopy("same");
    }
    if(C!=NULL) {
      fQuadratic->SetParameter(0,C->GetQFa());
      fQuadratic->SetParameter(1,C->GetQFb());
      fQuadratic->SetParameter(2,C->GetQFc());
      fQuadratic->SetRange( C->GetTime()-4, C->GetTime()+4 );
      fQuadratic->SetLineColor(kBlack);
      fQuadratic->DrawCopy("same");
    }
  }
  if(showfit==1) {
    if(A!=NULL) {
      fQuadratic->SetParameter(0,0);
      fQuadratic->SetParameter(1,A->GetLFa());
      fQuadratic->SetParameter(2,A->GetLFb());
      fQuadratic->SetRange( A->GetTime()-1, A->GetTime()+1 );
      fQuadratic->SetLineColor(kBlack);
      fQuadratic->DrawCopy("same");
    }
    if(B!=NULL) {
      fQuadratic->SetParameter(0,0);
      fQuadratic->SetParameter(1,B->GetLFa());
      fQuadratic->SetParameter(2,B->GetLFb());
      fQuadratic->SetRange( B->GetTime()-1, B->GetTime()+1 );
      fQuadratic->SetLineColor(kBlack);
      fQuadratic->DrawCopy("same");
    }
    if(C!=NULL) {
      fQuadratic->SetParameter(0,0);
      fQuadratic->SetParameter(1,C->GetLFa());
      fQuadratic->SetParameter(2,C->GetLFb());
      fQuadratic->SetRange( C->GetTime()-1, C->GetTime()+1 );
      fQuadratic->SetLineColor(kBlack);
      fQuadratic->DrawCopy("same");
    }
  }
  return ret;
}
//=================================
void testbeam::Init() {
  fList = new TList();
  fList->SetName( "output" );
  fList->SetOwner();
  fEvents = new TH1D( "fEvents","fEvents",10,-0.5,9.5 );
  fList->Add( fEvents );
  //--------------
  TList *listMCP = new TList();
  listMCP->SetName("mcp");
  listMCP->SetOwner();
  fMCP0 = new mcp("MCP0");
  listMCP->Add( fMCP0->GetList() );
  fEventDisplayMCP = new TCanvas("ED_MCP");
  fEventDisplayMCP->Divide(5,5);
  listMCP->Add( fEventDisplayMCP );
  fList->Add(listMCP);
  //--------------
  TList *listENE = new TList();
  listENE->SetName("energy");
  listENE->SetOwner();
  for(int i=0; i!=5; ++i) {
    fEnergy[i] = new energy( Form("SiPM%dE",i) );
    listENE->Add( fEnergy[i]->GetList() );
  }
  fEventDisplayEnergy[0] = new TCanvas( "ED_SiPME_3D" );
  fEventDisplayEnergy[0]->Divide(5,5);
  listENE->Add( fEventDisplayEnergy[0] );
  fEventDisplayEnergy[1] = new TCanvas( "ED_SiPME_2D" );
  fEventDisplayEnergy[1]->Divide(5,5);
  listENE->Add( fEventDisplayEnergy[1] );
  fEnergyCorrelation3D = new TH3D("EnergyCorrelation_3D",
                                  "EnergyCorrelation_3D;Bottom Right  [mV];Top Left  [mV];Bottom Left  [mV]",
                                  100,0,800,100,0,800,100,0,800);
  listENE->Add( fEnergyCorrelation3D );
  fEnergyCorrelation2D = new TH2D("EnergyCorrelation_2D",
                                  "EnergyCorrelation_2D;TEnergy Downstream  [mV];TEnergy Upstream  [mV]",
                                  100,0,800,100,0,800);
  listENE->Add( fEnergyCorrelation2D );
  fList->Add(listENE);
  //--------------
  TList *listTIM = new TList();
  listTIM->SetName("timing");
  listTIM->SetOwner();
  for(int i=0; i!=2; ++i) {
    fTiming[i] = new timing( Form("SiPM%dT",i) );
    listTIM->Add( fTiming[i]->GetList() );
  }
  fEventDisplayTiming = new TCanvas( "ED_SiPMT" );
  fEventDisplayTiming->Divide(5,5);
  listTIM->Add( fEventDisplayTiming );
  fTimingCorrelation2D = new TH2D("TimingCorrelation_2D",
                                  "TimingCorrelation_2D;Downstream Time - MCP  [nS];Upstream Time - MCP  [ns]",
                                  100,-3,+10,100,-3,+10);
  listTIM->Add( fTimingCorrelation2D );
  fATUSCorrelation2D = new TH2D("fATUSCorrelation_2D",
                                  "fATUSCorrelation_2D;Amplitude  [mV];Time - MCP  [ns]",
                                  100,0,400,100,-3,+10);
  listTIM->Add( fATUSCorrelation2D );
  fATDSCorrelation2D = new TH2D("fATDSCorrelation_2D",
                                  "fATDSCorrelation_2D;Amplitude  [mV];Time - MCP  [ns]",
                                  100,0,800,100,-3,+10);
  listTIM->Add( fATDSCorrelation2D );
  fList->Add(listTIM);
  fQuadratic = new TF1("fQuadratic","[0]*x*x+[1]*x+[2]");
}
//=================================
void testbeam::Terminate() {
  //std::cout << "Terminate() was called" << std::endl;
  TFile *outputfile = new TFile( fOutputFileName.Data(), "RECREATE" );
  fList->Write( fList->GetName(), TObject::kSingleKey );
  outputfile->Close();
  std::cout << "File " << fOutputFileName.Data() << " was saved." << std::endl;

}
//=================================
void testbeam::Process() {
  if(!fTree) return;
  Init();
  fTree->Print();
  int event;
  unsigned short tc[2]; // trigger counter bin
  float times[2][1024]; // calibrated time
  float channel[18][1024]; // calibrated input (in V)
  fTree->SetBranchAddress("event", &event);
  fTree->SetBranchAddress("tc", tc);
  fTree->SetBranchAddress("channel", channel);
  fTree->SetBranchAddress("times", times);

  Long64_t iEvent = 0;
  Long64_t iEvent1 = 0;
  Long64_t iEvent2 = 0;
  Long64_t iEvent3 = 0;
  Long64_t nEvents = fTree->GetEntries();
  fEvents->GetXaxis()->SetBinLabel(1,"RAW");
  fEvents->GetXaxis()->SetBinLabel(2,"MCP passed");
  fEvents->GetXaxis()->SetBinLabel(3,"ENERGY passed");
  fEvents->GetXaxis()->SetBinLabel(4,"TIMING passed");
  fEvents->GetXaxis()->SetBinLabel(5,"USE passed");
  DRSWAVE *ref = new DRSWAVE();
  for(;iEvent!=nEvents; ++iEvent) {
    fTree->GetEntry( iEvent );
    fEvents->Fill( 0 );
    if(iEvent%1000==0) std::cout << "Events read so far: " << iEvent << std::endl;



    fMCP0->Fill( times[0], channel[8] ); // group 0 channel 8
    fMCP0->ShiftWave( channel[7] ); // CMN
    bool passMCP0 =  fMCP0->Process();
    if(iEvent1<25) {
      fEventDisplayMCP->cd(iEvent1+1);
      TGraph *ret = DisplayWave(fMCP0, NULL, NULL, 1);
      ret->SetTitle( Form("Event  %lld",iEvent) );
    }
    double T0 = fMCP0->GetTime();
    double xtime[6];
    xtime[0] = 70 + 180 - T0;
    xtime[1] = 70 + 160 - T0;
    xtime[2] = 70 + 160 - T0;
    xtime[3] = 70 + 160 - T0;
    xtime[4] = 90 + 160 - T0;

    fEnergy[0]->Fill( times[0], channel[3] ); // group 1 channel 3 TRE
    fEnergy[1]->Fill( times[0], channel[4] ); // group 1 channel 4 BRE
    fEnergy[2]->Fill( times[0], channel[5] ); // group 1 channel 5 TLE
    fEnergy[3]->Fill( times[0], channel[6] ); // group 1 channel 6 BLE
    fEnergy[4]->Fill( times[1], channel[8+1] ); // group 1 channel 8+1 USE
    bool passEnergy[5] = {false,false,false,false,false};
    for(int i=0; i!=5; ++i) {
      fEnergy[i]->ShiftWave( channel[7] ); // CMN
      fEnergy[i]->SetRefTime( xtime[i] );
      passEnergy[i] =  fEnergy[i]->Process();
    }
    double amplitude[5];
    for(int i=0; i!=5; ++i) {
      //amplitude[i] = -fEnergy[i]->GetAmplitude();
      amplitude[i] = fEnergy[i]->GetAvgEne();
    }
    if(iEvent2<25) {
      fEventDisplayEnergy[0]->cd(iEvent2+1);
      TGraph *ret = DisplayWave(fEnergy[1], fEnergy[2], fEnergy[3], 2);
      ret->SetTitle( Form("Event  %lld",iEvent) );
      fEventDisplayEnergy[1]->cd(iEvent2+1);
      ret = DisplayWave(fEnergy[0], fEnergy[4], NULL, 2);
      ret->SetTitle( Form("Event  %lld",iEvent) );
    }

    fTiming[0]->Fill( times[0], channel[1] ); // group 0 channel 1 TRT
    fTiming[1]->Fill( times[0], channel[2] ); // group 0 channel 2 UST

    double fintime[2];
    bool passTiming[2] = {false,false};
    for(int i=0; i!=2; ++i) {
      fTiming[i]->SetThreshold( -300 ); // 20mV threshold
      passTiming[i] =  fTiming[i]->Process();
      fintime[i] = 108.5 + fTiming[i]->GetTime() - T0;
    }


    iEvent1++;
    iEvent2++;


    if( !passMCP0 ) continue;
    fEvents->Fill( 1 );
    if(!passEnergy[1]||!passEnergy[2]||!passEnergy[3]) continue;
    fEvents->Fill( 2 );
    if( !passTiming[0] || !passTiming[1] ) continue;
    fEvents->Fill( 3 );

    if(amplitude[4]<10) continue;
    fEvents->Fill( 4 );

    fEnergyCorrelation3D->Fill( amplitude[1], amplitude[2], amplitude[3] );
    fEnergyCorrelation2D->Fill( amplitude[0], amplitude[4] );
    fTimingCorrelation2D->Fill( fintime[0], fintime[1] );
    fATDSCorrelation2D->Fill( amplitude[0], fintime[0] );
    fATUSCorrelation2D->Fill( amplitude[4], fintime[1] );

    if(iEvent3<25) {
      fEventDisplayTiming->cd(iEvent3+1);
      TGraph *ret = DisplayWave(fTiming[0], fTiming[1], NULL ,1);
      ret->SetTitle( Form("Event  %lld",iEvent) );
    }
    iEvent3++;

  }
  Terminate();
  fTree = NULL; //don't call me again
}
