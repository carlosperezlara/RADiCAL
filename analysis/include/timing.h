#ifndef TIMING_H
#define TIMING_H

#include "DRSWAVE.h"
#include <TGraph.h>
#include <TList.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TString.h>

class timing : public DRSWAVE {
public:
  timing(TString name);
  ~timing();
  void SetThreshold(double val) {fThreshold = val;}
  bool Process();
  TList* GetList() {return fList;}

private:
  double fThreshold;
  
  TList *fList;
  TH1D *fBaseline_mean;
  TH1D *fBaseline_sdev;
  TH1D *fAmplitude_rawmaxx;
  TH1D *fAmplitude_rawmaxy;
  TH1D *fAmplitude_Thr;
  TH1D *fTime_Thr;
  TH1D *fHisto_LFa;
  TH1D *fHisto_LFb;
  TH2D *fHisto_LFab;
  TH1D *fHisto_LFrab;

};

#endif
