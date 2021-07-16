const int kUSE = 9;
const int kTRT = 2;
const int kUST = 1;
const int kMCP = 8;
const int kTRE = 3;
const int kBRE = 4;
const int kTLE = 5;
const int kBLE = 6;

const int kCOMM0 = 7;
const int kCOMM1 = 10;

//==================
void shiftY(TGraph *gr,double shift) {
  double ret=999;
  int n=gr->GetN();
  for(int i=0; i!=n; ++i) {
    double val = gr->GetPointY(i);
    gr->SetPointY(i,val+shift);
  }
}
//==================
double getmin(TGraph *gr, double &atime, int ini=0, int fin=-1) {
  double ret=999;
  int n=gr->GetN();
  if(fin>ini) n = fin;
  for(int i=ini; i!=n; ++i) {
    double val = gr->GetPointY(i);
    if(ret>val) {
      ret=val;
      atime = gr->GetPointX(i);
    }
  }
  return ret;
}
//==================
double getminAVG(TGraph *gr, double &atime, int ini=0, int fin=-1, int hwidth=2) {
  int n=gr->GetN();
  if(fin>ini) n = fin;
  int imiddle = -1;
  double ret=999;
  for(int i=ini; i!=n; ++i) {
    double val = gr->GetPointY(i);
    if(ret>val) {
      imiddle = i;
      ret = val;
    }
  }
  int ist = imiddle - hwidth;
  int ind = imiddle + hwidth;
  if(ist<0) ist = 0;
  if(ind>n) ind = n;
  atime = gr->GetPointX(imiddle);
  ret=0;
  for(int i=ist; i<=ind; ++i) {
    double wy = gr->GetPointY(i);
    double wx = gr->GetPointY(i);
    ret+=wy;
  }
  ret /= hwidth*2;
  return ret;
}
//==================
double getsdev(TGraph *gr, int ini=0, int fin=-1) {
  double ret=999;
  int n=gr->GetN();
  if(fin>ini) n = fin;
  double sx=0;
  double sxx=0;
  int vals=n-ini;
  if(vals<2) return -1;
  for(int i=ini; i!=n; ++i) {
    double val = gr->GetPointY(i);
    sx += val;
    sxx += val*val;
  }
  double sdev=sqrt(vals*sxx-sx*sx)/(vals*(vals-1));
  return sdev;
}
//==================
double getmean(TGraph *gr, int ini=0, int fin=-1) {
  double ret=999;
  int n=gr->GetN();
  if(fin>ini) n = fin;
  double sx=0;
  int vals=n-ini;
  if(vals<2) return -1;
  for(int i=ini; i!=n; ++i) {
    double val = gr->GetPointY(i);
    sx += val;
  }
  double mean = sx/vals;
  return mean;
}
//==================
//==================
//==================
int analysis(int run, Long64_t totalevents=100000, float UST_thr=400) {
  //gStyle->SetOptStat(1111111);
  gStyle->SetOptStat(0);
  TFile *file = new TFile(  Form  ("/Volumes/uva/testbeam_2021_06_data/merged/Run_%d.root",run) );
  //TFile *file = new TFile(  Form("Run_%d.root",run) );
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
  int chnidx[ nchn ] =     { kUSE, kTRT, kUST, kMCP, kTRE, kBRE, kTLE, kBLE};
  TString chnlab[ nchn ] = {"USE","TRT","UST","MCP","TRE","BRE","TLE","BLE"};

  TH2D *summary[nchn];
  TH1D *rms[nchn];
  TH1D *minimum[nchn];
  TH1D *minimumtime[nchn];
  TH1D *baseline_rms[nchn];
  TLatex *tex = new TLatex();
  TH1D *sumamp = new TH1D("sumamp","Three Capilary Sum;[mV]",200,0,2500);
  TProfile *proamp[nchn][4];
  int color[4] = {kCyan-3,kOrange-3,kGreen-3,kMagenta-3};
  for(int i=0; i!=nchn; ++i) {
    summary[i] = new TH2D(Form("STACK_%s", chnlab[i].Data() ),
                          Form("STACK_%s;ns;mV", chnlab[i].Data() ),
                          100,0.0,200.0, 100,-800.1,+200.1 );
    baseline_rms[i] = new TH1D(Form("BASELINE_RMS_%s", chnlab[i].Data() ),
                               Form("BASELINE_RMS_%s;rms{ mV }", chnlab[i].Data() ),
                               100,0.0,0.05 );
    rms[i] = new TH1D(Form("RMS_%s", chnlab[i].Data() ),
                      Form("RMS_%s;rms{ mV }", chnlab[i].Data() ),
                      100,0.0,25.0 );
    minimum[i] = new TH1D(Form("MIN_%s", chnlab[i].Data() ),
                          Form("MIN_%s;mV", chnlab[i].Data() ),
                          100,-1020,+20 );
    minimumtime[i] = new TH1D( Form("MINTIME_%s", chnlab[i].Data() ),
                              Form("MINTIME_%s;mV", chnlab[i].Data() ),
                              100,0,200 );
    for(int ii=0; ii!=4; ++ii) {
      proamp[i][ii] = new TProfile(Form("proamp_%d_%d",i,ii),"profile;[ns];[mV]",200,0,200);
      proamp[i][ii]->SetLineColor( color[ii] );
      proamp[i][ii]->SetMarkerColor( color[ii] );
    }
  }
  
  TString outputpdf = Form("analysis_Run%d.pdf",run);
  TCanvas *c0 = new TCanvas("c0",Form("RUN %d",run));
  TCanvas *c1 = new TCanvas("c1",Form("RUN %d",run));
  c1->Divide(4,2);
  c1->SaveAs( Form("%s[",outputpdf.Data()), "PDF" );

  Long64_t nentries = TMath::Min( totalevents, tree->GetEntries() );
  nentries = TMath::Max( nentries, Long64_t(0) );
  TGraph *trace[18];
  Long64_t nev=0;
  Long64_t ngoodev=0;
  for(;nev!=nentries;++nev) {
    tree->GetEntry( nev );
    if(nev%500==0)
      cout << "Events read so far: " << nev << endl;
    //get data from tree
    for(int i=0; i!=18; ++i) {
      int group = i/9;
      trace[i] = new TGraph( 1024, times[group], channel[i] );
    }
    //common mode noise
    for(int i=0; i!=nchn; ++i) {
      for(int is=0; is!=trace[kCOMM0]->GetN(); ++ is) {
        double newval = trace[ chnidx[i] ]->GetPointY(is);
        newval -= trace[kCOMM0]->GetPointY(is);
        trace[ chnidx[i] ]->SetPointY(is, newval);
      }
    }
    //baseline removal
    double bas_rms[nchn];
    double max_time[nchn];
    double max_amp[nchn];
    for(int i=0; i!=nchn; ++i) {
      double mean = getmean(trace[ chnidx[i] ],0,150);
      double sdev = getsdev(trace[ chnidx[i] ],0,150);
      baseline_rms[i]->Fill( sdev );
      bas_rms[i] = sdev;
      shiftY( trace[ chnidx[i] ], -mean );
    }
    double maxtime[nchn];
    double maxamp[nchn];
    for(int idx=0; idx!=nchn; ++idx) {
      maxamp[ idx ] = getminAVG( trace[ chnidx[ idx ] ], maxtime[idx] );
    }
    

    // ===================================
    // Access data through TRACE ( TGraph )

    //CUTS
    if(maxtime[1]<50 || maxtime[2]>90) continue;
    if(maxtime[1]<50 || maxtime[2]>90) continue;

    if(maxtime[5]<50 || maxtime[5]>90) continue;
    if(maxtime[6]<50 || maxtime[6]>90) continue;
    if(maxtime[7]<50 || maxtime[7]>90) continue;
    
    if(maxamp[2]>-UST_thr) continue;
      
    //double rmsUST = trace[ kUST ]->GetRMS(2);

    //double rmsBRE = trace[ kBRE ]->GetRMS(2);
    //double rmsTLE = trace[ kTLE ]->GetRMS(2);
    //double rmsBLE = trace[ kBLE ]->GetRMS(2);
    
    //if(rmsUST<10) continue;

    //if(rmsBRE<3) continue;
    //if(rmsTLE<3) continue;
    //if(rmsBLE<3) continue;

    ngoodev++;
    
    // FILLING HISTOGRAMS
    double energy = maxamp[5] + maxamp[6] + maxamp[7];
    sumamp->Fill( -energy );

    for(int idx=0; idx!=nchn; ++idx) {
      double tmp;
      double peak = -getmin( trace[ chnidx[idx] ], tmp );
      for(int sample=0; sample!=1024; ++sample) {
        double x = trace[ chnidx[ idx ] ]->GetPointX( sample );
        double y = trace[ chnidx[ idx ] ]->GetPointY( sample );
        summary[ idx ]->Fill( x, y );
        if(peak>50 && peak<200) proamp[idx][0]->Fill(x,y);
        if(peak>200 && peak<400) proamp[idx][1]->Fill(x,y);
        if(peak>400 && peak<600) proamp[idx][2]->Fill(x,y);
        if(peak>600) proamp[idx][3]->Fill(x,y);
      }
      rms[ idx ]->Fill( trace[ chnidx[ idx ] ]->GetRMS(2) );
      minimum[ idx ]->Fill( maxamp[idx] );
      minimumtime[ idx ]->Fill( maxtime[idx] );
    }

    if(ngoodev<50) {
      for(int idx=0; idx!=nchn; ++idx) {
        c1->cd(idx+1);
        trace[ chnidx[ idx ] ]->Draw("A*L");
        //cout << idx << " "  << trace[ chnidx[ idx ] ]->GetRMS(2) << endl;
        trace[ chnidx[ idx ] ]->SetTitle( chnlab[ idx ].Data() );
        trace[ chnidx[ idx ] ]->GetXaxis()->SetTitle( "ns" );
        trace[ chnidx[ idx ] ]->GetYaxis()->SetTitle( "mV" );
        tex->DrawLatexNDC(0.65,0.2,Form("min %.1f",maxamp[idx]));
      }
      c1->SaveAs( Form("%s[",outputpdf.Data()), "PDF" );
    }
    // ===================================
    
    for(int i=0; i!=18; ++i) delete trace[i];
  }
  cout << "total number of events read: " << nev << ". Events selected " << ngoodev << endl;

  for(int idx=0; idx!=nchn; ++idx) {
    c1->cd(idx+1);
    baseline_rms[idx]->Draw();
  }
  c1->SaveAs( Form("%s[",outputpdf.Data()), "PDF" );

  for(int idx=0; idx!=nchn; ++idx) {
    c1->cd(idx+1)->SetLogz(1);
    summary[idx]->Draw("colz");
  }
  c1->SaveAs( Form("%s[",outputpdf.Data()), "PDF" );

  TH2D *axis = new TH2D("axis","Pulse Shape;[ns];[mV]",100,0,200,100,-900,+300);
  TLegend *leg = new TLegend(0.7,0.1,0.9,0.4);
  leg->AddEntry(proamp[0][0],"[50-200]");
  leg->AddEntry(proamp[0][1],"[200-400]");
  leg->AddEntry(proamp[0][2],"[400-600]");
  leg->AddEntry(proamp[0][3],"[600-]");
  for(int idx=0; idx!=nchn; ++idx) {
    c1->cd(idx+1);
    axis->SetTitle( Form("PulseShape %s",chnlab[idx].Data()) );
    axis->DrawCopy();
    proamp[idx][3]->Draw("lsame");
    proamp[idx][2]->Draw("lsame");
    proamp[idx][1]->Draw("lsame");
    proamp[idx][0]->Draw("lsame");
    leg->Draw();
  }
  c1->SaveAs( Form("%s[",outputpdf.Data()), "PDF" );
    
  for(int idx=0; idx!=nchn; ++idx) {
    c1->cd(idx+1);
    rms[idx]->Draw();
  }
  c1->SaveAs( Form("%s[",outputpdf.Data()), "PDF" );

  for(int idx=0; idx!=nchn; ++idx) {
    c1->cd(idx+1);
    minimum[idx]->Draw();
  }
  c1->SaveAs( Form("%s[",outputpdf.Data()), "PDF" );

  for(int idx=0; idx!=nchn; ++idx) {
    c1->cd(idx+1);
    minimumtime[idx]->Draw();
  }
  c1->SaveAs( Form("%s[",outputpdf.Data()), "PDF" );

  c0->cd();
  sumamp->Draw();
  c0->SaveAs( Form("%s[",outputpdf.Data()), "PDF" );

  c1->SaveAs( Form("%s]",outputpdf.Data()), "PDF" );
  
  return 0;
}
