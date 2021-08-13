#ifndef MCP_H
#define MCP_H

#include "DRSWAVE.h"
#include <TGraph.h>
#include <TList.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TString.h>

class mcp : public DRSWAVE {
public:
  mcp(TString name);
  ~mcp();
  bool Process();
  TList* GetList() {return fList;}

private:
  TList *fList;
  TH1D *fBaseline_mean;
  TH1D *fBaseline_sdev;
  TH1D *fAmplitude_rawmaxx;
  TH1D *fAmplitude_rawmaxy;
  TH1D *fHisto_LFa;
  TH1D *fHisto_LFb;
  TH1D *fHisto_LFrab;
  TH2D *fTime_t500;
  TH2D *fTime_t501;
};

#endif
