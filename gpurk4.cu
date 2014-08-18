#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <cuda.h>

#define BLOCKSIZE 512

__device__ float gpu_compute_acc(float r, float J2, float *rp, float *cm, int N)
{
  float acc;
  float M = 0.0;
  int i = 0;
  float r2 = r*r;
  
  while ( (r > rp[i])&&(i<N))
    ++i;	
  
  if (i>0)
    M = cm[i-1];

  acc = -M/r2 + J2/(r*r2);
  return acc;
}

__global__ void gpu_rk4(float *r, float *vr, float *J2, float *rp, float *cm, int N, float dt)
{
  float kr[4], kv[4];
  float r0, v0, r1;
  int i; 

  i = threadIdx.x + blockIdx.x*blockDim.x;
  
  if (i < N)
    {
      r0 = r[i];
      v0 = vr[i];
      
      // Step 1
      r1 = r0;
      kr[0] = v0 ;
      kv[0] = gpu_compute_acc(r1, J2[i], rp, cm, N); 
      
      // Step 2        
      r1 = r0 + 0.5*dt*kr[0];
      kr[1] = v0 + 0.5*dt*kv[0];
      kv[1] = gpu_compute_acc(r1, J2[i], rp, cm, N); 
      
      
      // Step 3
      r1 = r0 + 0.5*dt*kr[1];
      kr[2] = v0 + 0.5*dt*kv[1];
      kv[2] = gpu_compute_acc(r1, J2[i], rp, cm, N); 
      
      
      // Step 4
      r1 = r0 + dt*kr[2];
      kr[3] = v0 + dt*kv[2];
      kv[3] = gpu_compute_acc(r1, J2[i], rp, cm, N);
      
      
      r[i] = r0  + dt*(kr[0] + 2.0*kr[1] + 2.0*kr[2] + kr[3])/6.0;
      vr[i] = v0 + dt*(kv[0] + 2.0*kv[1] + 2.0*kv[2] + kv[3])/6.0;
  }
  
}

extern "C" void rk4(float *r, float *vr, float *J2,  float *rp, float *cm, int N, float dt)
{
  float *r_d, *vr_d, *J2_d, *rp_d, *cm_d;
  
  cudaMalloc(&r_d  , sizeof(float)*N); 
  cudaMalloc(&vr_d , sizeof(float)*N); 
  cudaMalloc(&J2_d , sizeof(float)*N); 
  cudaMalloc(&rp_d , sizeof(float)*N); 
  cudaMalloc(&cm_d , sizeof(float)*N); 
  
  cudaMemcpy(r_d,  r , sizeof(float)*N, cudaMemcpyHostToDevice); // Host -> Device
  cudaMemcpy(vr_d, vr, sizeof(float)*N, cudaMemcpyHostToDevice); // Host -> Device
  cudaMemcpy(J2_d, J2, sizeof(float)*N, cudaMemcpyHostToDevice); // Host -> Device
  cudaMemcpy(rp_d, rp, sizeof(float)*N, cudaMemcpyHostToDevice); // Host -> Device
  cudaMemcpy(cm_d, cm, sizeof(float)*N, cudaMemcpyHostToDevice); // Host -> Device

  gpu_rk4 <<<256,BLOCKSIZE >>>(r_d, vr_d, J2_d, rp_d, cm_d, N, dt);

  cudaMemcpy(r , r_d, sizeof(float)*N, cudaMemcpyDeviceToHost); // Device -> Host
  cudaMemcpy(vr, vr_d, sizeof(float)*N, cudaMemcpyDeviceToHost); // Device -> Host

  cudaFree(r_d);
  cudaFree(vr_d);
  cudaFree(J2_d);
  cudaFree(rp_d);
  cudaFree(cm_d);

  return;
}

