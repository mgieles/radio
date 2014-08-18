#include <stdio.h>
#include <stdlib.h>
#include <math.h>

float compute_acc(float r, float J2, float *rp, float *cm, int N);

void rk4(float *r, float *vr, float *J2, float *rp, float *cm,  int N, float dt)
{
  float kr[4], kv[4];
  float r0, v0, r1;
  for (int i = 0; i<N; ++i)
    {
      float t = 0;
      while (t < dt)
	{
	  r0 = r[i];
	  v0 = vr[i];

	  /* Step 1 */
	  r1 = r0;
	  kr[0] = v0 ;
	  kv[0] = compute_acc(r1, J2[i], rp, cm, N);

	  /* Step 2 */
	  r1 = r0 + 0.5*dt*kr[0];
	  kr[1] = v0 + 0.5*dt*kv[0];
	  kv[1] = compute_acc(r1, J2[i], rp, cm, N);
	  
	  /* Step 3 */
	  r1 = r0 + 0.5*dt*kr[1];
	  kr[2] = v0 + 0.5*dt*kv[1];
	  kv[2] = compute_acc(r1, J2[i], rp, cm, N);
	  
	  /* Step 4 */
	  r1 = r0 + dt*kr[2];
	  kr[3] = v0 + dt*kv[2];
	  kv[3] = compute_acc(r1, J2[i], rp, cm, N);

	  r[i] = r0  + dt*(kr[0] + 2.0*kr[1] + 2.0*kr[2] + kr[3])/6.0;
	  vr[i] = v0 + dt*(kv[0] + 2.0*kv[1] + 2.0*kv[2] + kv[3])/6.0;
	  t+=dt;
	}
    }
}


float compute_acc(float r, float J2, float *rp, float *cm, int N)
{
  /***
      The equation of motion in a spherically symmetric potential
      dr/dt = v_r
      dv_r/dt = -dPsi/dr, with Psi = phi + J^2/2r^2 = effective potential
              = -dphi/dr + J^2/r^3
              = -GM(<r)/r^2 + J^2/r^3
  ***/

  float M = 0.0;
  float r2 = pow(r,2.0);
  int i = 0;

  while ((r > rp[i])&&(i<N))
    ++i;

  if (i>0)
    M = cm[i-1];
  float acc = -M/r2 + J2/(r*r2);
  return acc;
}


