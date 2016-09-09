#include "common.h"
#include "debug.h"
#include "gundersen.h"
#include "matrix.h"
#include "tensor.h"
#include "vector.h"


double valueH[] = {1,1,2,3,1,4,2,5,1,1,1,6};
int indexH[] = {0,0,1,2,1,3,2,4,0,2,4,5};
int pointerH[] = {0,1,3,4,6,8,12};

double valueT[] = {1,2,2,2,3,4,4,4,5,5,5,6,6,6,6,6,6,6,6};
int indexT[] = {0,0,0,1,2,1,1,3,2,2,4,0,2,2,4,0,2,4,5};
int pointerT[] = {0,1,2,4,5,6,8,9,11,12,13,15,19};


void gundersen() {
  int N = 6;
  double p[] = {1,1,1,1,1,1};
  double g[] = {0,0,0,0,0,0};
  int start = 0, stop = 0, i = 0, j = 0, k = 0;
  int ind = 0, tp = 0;
  int pi = 0, pj = 0, pk = 0, pipj = 0, kpipj = 0, kpi = 0, kpj = 0, pipi = 0, pkTijk=0;
  int startTubek=0, stopTubek=0;
  double Tijk = 0, Tijj = 0, Tiii = 0, Tiik = 0;
  double ga = 0;
  for(i = 0;i<N;i++,ind++,tp++){
    start = pointerH[i];
    stop = pointerH[i+1]-1;
    j = indexH[start];
    pi = p[i];
    kpi = 2*pi;
    pipi = pi*pi;
    for(;j<i;start++,ind++,tp++){
      j = indexH[start];
      pj = p[j];
      pipj = pi*pj;
      kpipj = 2*pipj;
      kpj = 2*pj;
      startTubek = pointerT[tp];
      stopTubek = pointerT[tp+1]-1;
      for(;startTubek<stopTubek;startTubek++,ind++){
        //Handle the case when no indices are equal: i!=j!=k!=i
        k = indexT[ind];
        Tijk = valueT[ind];
        pk = p[k];
        pkTijk = pk*Tijk;
        ga += pkTijk;
        g[k] += kpipj*Tijk;
      }
      //Handle the case when two indices are equal: k=j
      Tijj = valueT[ind];
      g[i] += (kpj*ga+pj*pj*Tijj);
      g[j] += (kpi*ga+kpipj*Tijj);
      ga = 0.0;
      j = indexH[start+1];
    }
    startTubek = pointerT[tp];
    stopTubek = pointerT[tp+1]-1;
    for(;startTubek<stopTubek;startTubek++,ind++){
      //Handle the case when two indices are equal: j=i
      k = indexT[ind];
      Tiik = valueT[ind];
      pk = p[k];
      ga += pk*Tiik;
      g[k] += pipi*Tiik;
    }
    //Handle the case when all three indices are equal
    g[i] += (kpi*ga + pipi*valueT[ind]);
    ga = 0.0;
  }
}
