#ifndef TESTBEAM_H
#define TESTBEAM_H
#include <TString.h>
#include <TTree.h>
#include <TList.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TProfile.h>
#include <TGraph.h>

#include "mcp.h"
#include "energy.h"
#include "timing.h"

class testbeam {
public:
  testbeam();
  ~testbeam();
  void SetInputFileName(TString);
  void SetOutputFileName(TString);
  void Process();

private:
  void Init();
  void Terminate();
  void DisplayWave(int idx, int iEvent2, DRSWAVE *A, DRSWAVE *B);
  TGraph* DisplayWave(DRSWAVE *A, DRSWAVE *B, DRSWAVE *C, int showfit=0);

  TTree *fTree;
  TList *fList;
  TString fOutputFileName;

  mcp *fMCP0;
  energy *fEnergy[5];
  timing *fTiming[2];

  TH1D *fEvents;

  TCanvas *fEventDisplayMCP;
  TCanvas *fEventDisplayEnergy[2];
  TCanvas *fEventDisplayTiming;

  TH3D *fEnergyCorrelation3D;
  TH2D *fEnergyCorrelation2D;

  TH2D *fTimingCorrelation2D;
  TH2D *fATDSCorrelation2D;
  TH2D *fATUSCorrelation2D;

  TF1 *fQuadratic;
};

#endif
