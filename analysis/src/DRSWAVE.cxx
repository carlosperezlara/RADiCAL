#include "DRSWAVE.h"
#include <TList.h>
#include <TF1.h>

//=================================
DRSWAVE::DRSWAVE() : TGraph(1024) {
}
//=================================
DRSWAVE::~DRSWAVE() {
}
//=================================
void DRSWAVE::Fill(float x[1024], float y[1024]) {
  for(int i=0; i!=1024; ++i) {
    SetPointX(i,x[i]);
    SetPointY(i,y[i]);
  }
}
//=================================
void DRSWAVE::GetBaseline(double &mean, double &sdev, int irange, int frange) {
  double sx=0;
  double sxx=0;
  int vals=frange-irange;
  for(int i=irange; i!=frange; ++i) {
    double val = GetPointY(i);
    sx += val;
    sxx += val*val;
  }
  mean = sx/vals;
  sdev=sqrt(vals*sxx-sx*sx)/(vals*(vals-1));
}
//=================================
void DRSWAVE::GetSampleSum(double &sum, int irange, int frange) {
  sum=0;
  for(int i=irange; i!=frange; ++i) {
    double val = GetPointY(i);
    sum += val;
  }
}
//=================================
int DRSWAVE::GetMinimum(double &ymin, double &xmin, int irange, int frange, int idelta) {
  double temp = 2000;
  int imiddle = -1;
  for(int i=irange; i!=frange; ++i) {
    double val = GetPointY(i);
    if(temp>val) {
      imiddle = i;
      temp = val;
    }
  }

  int ist = imiddle - idelta;
  int ind = imiddle + idelta +1;
  if(ist<0) ist = 0;
  if(ind>1024) ind = 1024;
  int vals = ind - ist;
  double sxy = 0;
  double sy = 0;
  for(int i=ist; i!=ind; ++i) {
    double y = GetPointY(i);
    double x = GetPointX(i);
    sxy += x * y;
    sy += y;
  }
  ymin = sy / vals;
  xmin = sxy / sy;
  return imiddle;
}
//=================================
void DRSWAVE::ShiftWave(double cte) {
  for(int i=0; i!=1024; ++i) {
    double val = GetPointY(i);
    SetPointY(i,val+cte);
  }
}
//=================================
void DRSWAVE::ShiftWave(double ref[1024]) {
  for(int i=0; i!=1024; ++i) {
    double val = GetPointY(i);
    SetPointY(i,val-ref[i]);
  }
}
//=================================
void DRSWAVE::ShiftWave(float ref[1024]) {
  for(int i=0; i!=1024; ++i) {
    double val = GetPointY(i);
    SetPointY(i,val-ref[i]);
  }
}
//=================================
void DRSWAVE::ShiftWave(TGraph *ref) {
  for(int i=0; i!=1024; ++i) {
    double val = GetPointY(i);
    SetPointY(i,val+ref->GetPointY(i));
  }
}
//=================================
void DRSWAVE::FindThreshold_P1(double &x, double thr, int leftpoints, int rightpoints) {
  int i=0;
  for(; i!=1024; ++i) {
    double val = GetPointY(i);
    if(val<thr) break;
  }
  x=-999;
  if(i==1023) return;
  int imin = i-leftpoints;
  int imax = i+rightpoints-1;
  if(0) {
    Fit( "pol1", "rq", "", GetPointX(imin), GetPointX(imax) );
    x = ((TF1*)GetListOfFunctions()->At(0))->GetX( thr, GetPointX(imin), GetPointX(imax) );
  } else {
    LinearFit(imin,imax+1);
    x = (thr - fLFb) / fLFa;
  }
}
//=================================
void DRSWAVE::FindMinimum_P2(double &ymin, double &xmin, int ini, int fin) {
  QuadraticFit(ini,fin);
  xmin = -fQFb/fQFa/2.0;
  //ymin = fQFb*(fQFc/fQFb+xmin/2.0);
  ymin = fQFa*( xmin*xmin + fQFb/fQFa*xmin + fQFc/fQFa);
}
//=================================
void DRSWAVE::LinearFit(int ini, int fin) {
  double Sx=0;
  double Sy=0;
  double Sxx=0;
  double Sxy=0;
  double S=0;
  for(int i=ini; i!=fin; ++i) {
    double x = GetPointX( i );
    double y = GetPointY( i );
    //double dx = gr->GetErrorX( i );
    double dy = 1;//gr->GetErrorY( i );
    S += 1/dy;
    Sx += x/dy;
    Sy += y/dy;
    Sxx += x*x/dy;
    Sxy += x*y/dy;
  }
  double den = Sx*Sx - S*Sxx;
  fLFa = (Sy*Sx - Sxy*S) / den;
  fLFb = (Sxy*Sx - Sy*Sxx) / den;
  double var_a = -S / den;
  double var_b = -Sxx / den;
  double covar = Sx / den;
  fLFrab = covar / TMath::Sqrt( var_a*var_b );

  // correction for sigma_i = 1
  if(0) {
    double chi2 = 0;
    for(int i=ini; i!=fin; ++i) {
      double x = GetPointX( i );
      double y = GetPointY( i );
      chi2 += TMath::Power( (y-fLFb-fLFa*x)/1, 2 );
    }
    chi2 /= (fin-ini-2);
    var_a *= chi2;
    var_b *= chi2;
  }
}
//=================================
void DRSWAVE::QuadraticFit(int ini, int fin) {
  double Sx=0;
  double Sy=0;
  double Sxx=0;
  double Sxy=0;
  double Sxxx=0;
  double Sxxxx=0;
  double Sxxy=0;
  double S=0;
  for(int i=ini; i!=fin; ++i) {
    double x = GetPointX( i );
    double y = GetPointY( i );
    S += 1;
    Sx += x;
    Sy += y;
    Sxx += x*x;
    Sxxx += x*x*x;
    Sxxxx += x*x*x*x;
    Sxy += x*y;
    Sxxy += x*x*y;
  }
  /// Sxxxx  Sxxx  Sxx     |  a  |    Sxxy
  /// Sxxx    Sxx    Sx      |  b  |    Sxy
  /// Sxx      Sx      S       |  c  |     Sy
  double det = Sxxxx*(Sxx*S-Sx*Sx) - Sxxx*(Sxxx*S-Sx*Sxx) + Sxx*(Sxxx*Sx-Sxx*Sxx);
  /// S Sxx - Sx Sx         | Sx Sxx - S Sxxx        | Sx Sxxx - Sxx Sxx
  /// Sx Sxx - S Sxxx     | S Sxxxx - Sxx Sxx     | Sxx Sxxx - Sx Sxxxx
  /// Sx Sxxx - Sxx Sxx  | Sxx Sxxx - Sx Sxxxx | Sxx Sxxxx - Sxxx Sxxx
  fQFa = ( Sxxy*(Sxx*S-Sx*Sx)     + Sxy*(Sxx*Sx-Sxxx*S)     + Sy*(Sxxx*Sx-Sxx*Sxx) )     / det;
  fQFb = ( Sxxy*(Sx*Sxx-Sxxx*S)   + Sxy*(Sxxxx*S-Sxx*Sxx)   + Sy*(Sxx*Sxxx-Sxxxx*Sx) )   / det;
  fQFc = ( Sxxy*(Sxxx*Sx-Sxx*Sxx) + Sxy*(Sxxx*Sxx-Sxxxx*Sx) + Sy*(Sxxxx*Sxx-Sxxx*Sxxx) ) / det;
}

