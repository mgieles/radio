#include <thrust/sort.h>
#include <thrust/device_vector.h>
#include <thrust/device_ptr.h>

extern "C" void sort_func(float *r, int *jlist, int N)
{
  float *r_p;
  int *jlist_p;
  
  cudaMalloc(&r_p    , N*sizeof(float));
  cudaMalloc(&jlist_p, N*sizeof(int));
  
  cudaMemcpy(r_p    ,  r    , sizeof(float)*N, cudaMemcpyHostToDevice); // Host -> Device
  cudaMemcpy(jlist_p,  jlist, sizeof(int)*N, cudaMemcpyHostToDevice); // Host -> Device
  
  thrust::device_ptr<float> r_d(r_p);
  thrust::device_ptr<int> jlist_d(jlist_p);
  
  thrust::sort_by_key(r_d, r_d + N, jlist_d);
  
  r_p     = thrust::raw_pointer_cast(r_d);
  jlist_p = thrust::raw_pointer_cast(jlist_d);
  
  cudaMemcpy(r    , r_p, sizeof(float)*N, cudaMemcpyDeviceToHost); // Device -> Host
  cudaMemcpy(jlist, jlist_p, sizeof(int)*N, cudaMemcpyDeviceToHost); // Device -> Host
  
  cudaFree(r_p);
  cudaFree(jlist_p);
}
