#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <cuda.h>

#define BLOCKSIZE 512

__global__ void gpu_phi(float *r, float *m, float *phi, int N)
{
  int i; 

  i = threadIdx.x + blockIdx.x*blockDim.x;
  
  if (i < N)
    {
      phi[i] = 0.0;
      for (int j=0; j<i; ++j)
	phi[i] -= m[j]/r[i];
      
      for (int j=i+1; j<N; ++j)
	phi[i] -= m[j]/r[j];
    }
}

extern "C" void getphi_func(float *r, float *m, float *phi, int N)
{
  float *r_d, *m_d, *phi_d;
  
  cudaMalloc(&r_d   , sizeof(float)*N); 
  cudaMalloc(&m_d   , sizeof(float)*N); 
  cudaMalloc(&phi_d , sizeof(float)*N); 
  
  cudaMemcpy(r_d  , r  , sizeof(float)*N, cudaMemcpyHostToDevice); // Host -> Device
  cudaMemcpy(m_d  , m  , sizeof(float)*N, cudaMemcpyHostToDevice); // Host -> Device
  cudaMemcpy(phi_d, phi, sizeof(float)*N, cudaMemcpyHostToDevice); // Host -> Device

  //  gpu_phi <<< 256, BLOCKSIZE >>>(r_d, m_d, phi_d, N);
  gpu_phi <<< ((N+BLOCKSIZE-1))/BLOCKSIZE,BLOCKSIZE >>>(r_d, m_d, phi_d, N);
  
  cudaMemcpy(phi, phi_d, sizeof(float)*N, cudaMemcpyDeviceToHost); // Device -> Host

  cudaFree(r_d);
  cudaFree(m_d);
  cudaFree(phi_d);

  return;
}

