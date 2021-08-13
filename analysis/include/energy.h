#ifndef ENERGY_H
#define ENERGY_H

#include "DRSWAVE.h"
#include <TGraph.h>
#include <TList.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TString.h>

class energy : public DRSWAVE {
public:
  energy(TString name);
  ~energy();
  bool Process();
  TList* GetList() {return fList;}
  double GetAvgEne() {return fAvgEne;}
  void SetRefTime(double val) {fRefTime=val;}

private:
  double fAvgEne;
  double fRefTime;
  
  TList *fList;
  TH1D *fBaseline_mean;
  TH1D *fBaseline_sdev;
  TH1D *fAmplitude_rawmaxx;
  TH1D *fAmplitude_rawmaxy;
  TH1D *fSampleSum;
  TH1D *fHisto_QFa;
  TH1D *fHisto_QFb;
  TH1D *fHisto_QFc;
  TH1D *fHisto_QFx0;
  TH1D *fHisto_QFy0;
  TH2D *fHisto_QFy0rawmaxy;
};

#endif
