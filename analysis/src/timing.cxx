#include "timing.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TString.h>

//=================================
timing::timing(TString name) {
  fList = new TList();
  fList->SetName( name.Data() );
  fList->SetOwner();
  fBaseline_mean = new TH1D(Form("%s_BL_mean",name.Data()),
                            Form("%s_BL_mean;[mV]",name.Data()),
                            100,-500,+500);
  fList->Add( fBaseline_mean );
  fBaseline_sdev = new TH1D(Form("%s_BL_sdev",name.Data()),
                            Form("%s_BL_sdev;[mV]",name.Data()),
                            100,0,0.1);
  fList->Add( fBaseline_sdev );
  fAmplitude_rawmaxx = new TH1D(Form("%s_A_rmaxx",name.Data()),
                                Form("%s_A_rmaxX;[ns]",name.Data()),
                                100,0,200);
  fList->Add( fAmplitude_rawmaxx );
  fAmplitude_rawmaxy = new TH1D(Form("%s_A_rmaxy",name.Data()),
                                Form("%s_A_rmaxY;[mV]",name.Data()),
                                100,0,1000);
  fList->Add( fAmplitude_rawmaxy );
  fAmplitude_Thr = new TH1D(Form("%s_A_thr",name.Data()),
                                Form("%s_A_thr;[mV]",name.Data()),
                                100,0,1000);
  fList->Add( fAmplitude_Thr );
  fTime_Thr = new TH1D(Form("%s_T_thr",name.Data()),
                            Form("%s_T_thr;[ns]",name.Data()),
                            100,0,200);
  fList->Add( fTime_Thr );
  fHisto_LFa = new TH1D(Form("%s_LFa",name.Data()),
                        Form("%s_LFa;[mV/ns]",name.Data()),
                        100,-1000,+1000);
  fList->Add( fHisto_LFa );
  fHisto_LFb = new TH1D(Form("%s_LFb",name.Data()),
                        Form("%s_LFb;[mV]",name.Data()),
                        100,-1000,+100000);
  fList->Add( fHisto_LFb );
  fHisto_LFrab = new TH1D(Form("%s_LFrab",name.Data()),
                          Form("%s_fHisto_LFrab",name.Data()),
                          100,-2,+2);
  fList->Add( fHisto_LFrab );
}
//=================================
timing::~timing() {
  delete fList;
}
//=================================
bool timing::Process() {
  SetAmplitude(-9999);
  SetTime(-9999);

  double bl_mean, bl_sdev;
  GetBaseline( bl_mean, bl_sdev, 0, 200 );
  fBaseline_mean->Fill( bl_mean );
  fBaseline_sdev->Fill( bl_sdev );

  ShiftWave( -bl_mean );

  double rawmaxy, rawmaxx;
  int imiddle = GetMinimum( rawmaxy, rawmaxx, 0, 1024, 0);
  
  fAmplitude_rawmaxx->Fill( rawmaxx );
  fAmplitude_rawmaxy->Fill( -rawmaxy );

  SetAmplitude( rawmaxy );
  SetTime(rawmaxx);
  
  fAmplitude_Thr->Fill( -fThreshold );
  double timeX;
  FindThreshold_P1( timeX, fThreshold, 1, 1);
  fTime_Thr->Fill( timeX );
  fHisto_LFa->Fill( GetLFa() );
  fHisto_LFb->Fill( GetLFb() );
  fHisto_LFrab->Fill( GetLFrab() );
  if(GetLFrab()>-0.8) return false;
  if(GetLFb()<0) return false;
  SetTime(timeX); // first order
  
  // improving timing

  return true;
}
