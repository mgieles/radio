#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Compute potential on host with double loop

void getphi_func(float *r, float *m, float *phi, int N)
{
  for (int i=0; i<N; ++i){
    phi[i] = 0.0;
    for (int j=0; j<i; ++j)
      phi[i] -= m[j]/r[i];
    for (int j=i+1; j<N; ++j)
      phi[i] -= m[j]/r[j];
  }
}
