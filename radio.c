#include <sys/types.h>
#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "radio.h"
#define BUFFERSIZE 300
#define LARGE_N 10000000
#ifndef M_PI
#define M_PI  3.14159265358979323846
#endif


/*************************************************/
void get_args(int argc, char** argv, INPUT *params)
{
  // User defined parameters
  params->dt    = pow(2.0,-9.0);    
  params->dtadj = pow(2.0,-9.0);
  params->dtout = pow(2.0,-3.0);
  params->tend  = 1;    

  for(int i = 1; i < argc; i++){
    if (argv[i][0] == '-'){
      switch (argv[i][1]){
      case 't': params->tend = atof(argv[++i]);
	break;
      case 'd': params->dt = pow(2.0,-atof(argv[++i]));
	break;
      case 'a': params->dtadj = pow(2.0,-atof(argv[++i]));
	break;
      }
    }  
  }
  params->tout     = params->dtout ; 
  params->tadjout  = 0.0; 
  params->tadj     = params->dtadj ; 
  params->wtime0   = wtime();
}

/*************************************************/
void readdata(CLUSTER **cluster)
{
  STAR *stars = malloc(sizeof(STAR));
  double** data = malloc(7 * sizeof(double*));      // allocate the rows
  for (int i = 0; i < 7; ++i)
    data[i] = malloc(LARGE_N * sizeof(double));   // allocate the columns
  char *tmp;
  int i=0,k;
  char buffer[BUFFERSIZE];
  float r2, v2;

  // Read data in data array
  while(fgets(buffer, BUFFERSIZE, stdin)) {
    tmp = buffer;
    for (int j=0; j<7;j++)
      data[j][i] = strtod(tmp, &tmp);
    i++;
  } 

  // Allocate memory for cluster and stars
  *cluster = malloc(sizeof(CLUSTER));
  (*cluster)->N = i;
  (*cluster)->t = 0.0;
  (*cluster)->stars = calloc((*cluster)->N, sizeof(*stars));

  // Copy data to cluster and convert to spherical coordinates
  for (i=0; i<(*cluster)->N; i++)
    {
      (*cluster)->stars[i].mass = data[0][i];

      r2 = 0;
      v2 = 0;
      (*cluster)->stars[i].vr  = 0.0; 
      for (k=0; k<3; k++){
	r2 += pow(data[1+k][i],2.0);
	v2 += pow(data[4+k][i],2.0);
	(*cluster)->stars[i].vr += data[1+k][i]*data[4+k][i];
      }
      (*cluster)->stars[i].r = sqrt(r2);
      (*cluster)->stars[i].vr /= (*cluster)->stars[i].r;
      (*cluster)->stars[i].vt = sqrt(v2 - pow((*cluster)->stars[i].vr,2.0));
      (*cluster)->stars[i].J2 = pow((*cluster)->stars[i].r * (*cluster)->stars[i].vt, 2.0);    
      (*cluster)->stars[i].E0 = 0.5*v2; // Add potential later
    }

  for (i=0; i<(*cluster)->N; ++i){
    (*cluster)->stars[i].id = i+1;
  }

  sort(*cluster);
  getphi(*cluster);

  for (i=0; i<(*cluster)->N; i++){
    (*cluster)->stars[i].E0 += (*cluster)->stars[i].phi;
    (*cluster)->stars[i].E = (*cluster)->stars[i].E0;

    /*
      Individual time step based on radial period: NOTE USED
      float dt = (2.0*M_PI/pow(-2.0*(*cluster)->stars[i].E0,1.5))/1024.0;
      float dt_discrete = pow(2.0,-3.0);      
      while(dt_discrete > dt)
      dt_discrete /= 2.0;
      (*cluster)->stars[i].dt = pow(2.0, -7.0); //dt_discrete;
    */

  }

  for (i=0; i<(*cluster)->N; ++i){
    (*cluster)->stars[i].id = i+1;
  }

  free(data);  
}

/*************************************************/
void adjust(CLUSTER *cluster, INPUT *params)
{
  float mvr2 = 0.0, mvt2 = 0.0, vr2, vt2;

  if (cluster->t > 0.0)
    {
      sort(cluster);
      getphi(cluster);
    }

  // Energy
  cluster->M = 0.0;
  cluster->E = 0.0;
  cluster->K = 0.0;
  cluster->W = 0.0;
  cluster->J = 0.0;

  // Cumulative mass 
  cluster->stars[0].cmass = cluster->stars[0].mass;

  for (int i=0; i<cluster->N; ++i){
    vr2 = pow(cluster->stars[i].vr,2.0);
    vt2 = pow(cluster->stars[i].vt,2.0);

    mvr2 += vr2;
    mvt2 += vt2;

    cluster->M += cluster->stars[i].mass;
    cluster->W += 0.5*cluster->stars[i].mass*cluster->stars[i].phi;
    cluster->K += 0.5*cluster->stars[i].mass*(vr2 + vt2);
    cluster->J += cluster->stars[i].mass*cluster->stars[i].vt*cluster->stars[i].r;
    cluster->stars[i].E = 0.5*(vr2 + vt2) + cluster->stars[i].phi;
    cluster->stars[i].rp = cluster->stars[i].r;
    
    if (i>0)
      cluster->stars[i].cmass = cluster->stars[i-1].cmass + cluster->stars[i].mass;
  }

  // Half-mass radius
  cluster->rh=0.0;

  for (int i=0; i<cluster->N; ++i){
    if ((cluster->rh == 0.0) && (cluster->stars[i].cmass >= 0.5*cluster->M))
      cluster->rh = cluster->stars[i].r;
  }

  // Some global diagnostics to error stream
  if (cluster->t >= params->tadjout){
    fprintf(stderr, " t = %9.4f  M = %5.3f  K = %8.5f  W = %8.5f  E = %8.5f  J = %8.5f  <vr2> = %7.5f  <vt2> = %7.5f  rh = %6.3f  wtime = %9.3f sec \n",
	    cluster->t, cluster->M, cluster->K, cluster->W,  cluster->W+ cluster->K, cluster->J, mvr2/(float)cluster->N, mvt2/(float)cluster->N, cluster->rh, wtime()-params->wtime0);
    params->tadjout += params->dtout;
  }
}

/*************************************************/
void sort(CLUSTER *cluster)
{
  // Sort stars in order of increasing r on CPU or GPU

  CLUSTER *cl;
  cl = malloc(sizeof(CLUSTER));
  cl->stars = calloc(cluster->N, sizeof(STAR));

  // Copy data to arrays for (possible) use of GPU
  int N = (int)cluster->N;
  float *r = (float *)malloc(N*sizeof(float));
  int *jlist = (int *)malloc(N*sizeof(int));

  for (int i=0; i<cluster->N; ++i){
    cl->stars[i] = cluster->stars[i];
    r[i] = (float)cluster->stars[i].r;
    jlist[i] = i;
  }

  sort_func(r, jlist, N);

  // Put the stars in correct order
  for (int i=0; i<cluster->N; ++i){
    cluster->stars[i] = cl->stars[jlist[i]];
  }
}


/*************************************************/
void getphi(CLUSTER *cluster)
{
  // Compute potential on CPU or GPU
  int N = (int)cluster->N;

  float *r = (float *)malloc(N*sizeof(float));
  float *m = (float *)malloc(N*sizeof(float));
  float *phi = (float *)malloc(N*sizeof(float));

  for (int i=0; i<cluster->N; ++i){
    r[i] = (float)cluster->stars[i].r;
    m[i] = (float)cluster->stars[i].mass;
    phi[i] = 0.0;
  }

  getphi_func(r, m, phi, N);

  for (int i=0; i<cluster->N; ++i){
    cluster->stars[i].phi = phi[i];
  }
  free(r);
  free(m);
  free(phi);
}


/*************************************************/
void integrate(CLUSTER *cluster, INPUT *params)
{
  // Output at t = 0;
  adjust(cluster, params);
  output(cluster);

  // Integrate cluster
  while (cluster->t < params->tend)
    {
      // Copy data to arrays for (possible) use of GPU
      int N = (int)cluster->N;
      float *r = (float *)malloc(N*sizeof(float));
      float *vr = (float *)malloc(N*sizeof(float));
      float *rp = (float *)malloc(N*sizeof(float));
      float *cm = (float *)malloc(N*sizeof(float));
      float *J2 = (float *)malloc(N*sizeof(float));
      
      for (int i=0; i<cluster->N; ++i){
	r[i] = (float)cluster->stars[i].r;
	vr[i] = (float)cluster->stars[i].vr;
	rp[i] = (float)cluster->stars[i].r;
	cm[i] = (float)cluster->stars[i].cmass;
	J2[i] = (float)cluster->stars[i].J2;
      }
      
      // Take RK4 step, can be on GPU
      rk4(r, vr, J2, rp, cm, cluster->N, params->dt); 

      // Copy data back to cluster
      for (int i=0; i<cluster->N; ++i){
	cluster->stars[i].r = r[i];
	cluster->stars[i].vr = vr[i];
	cluster->stars[i].vt = sqrt(J2[i])/r[i];
      }
      
      free(r);
      free(vr);
      free(rp);
      free(cm);
      free(J2);
      free(dt);
      
      cluster->t += params->dt;
      
      if (cluster->t >= params->tadj){
	params->tadj += params->dtadj;
	adjust(cluster, params);
      }
      
      if (cluster->t >= params->tout){
	params->tout += params->dtout;
	output(cluster);
      }
    }
}


/*************************************************/
void output(CLUSTER *cluster)
{
  for (int i=0;i<cluster->N;i++){
    float dE = (cluster->stars[i].E - cluster->stars[i].E0)/cluster->stars[i].E0; 
    printf("%10.3e %7d %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e \n",
	   cluster->t,                    // (1)
	   cluster->stars[i].id,          // (2)
	   cluster->stars[i].mass,        // (3)
	   cluster->stars[i].r,           // (4)
	   cluster->stars[i].vr,          // (5)
	   cluster->stars[i].vt,          // (6)
	   cluster->stars[i].phi,         // (7)
	   cluster->stars[i].E,           // (8)
	   cluster->stars[i].J2,          // (9)
	   cluster->stars[i].cmass,       // (10)
	   cluster->stars[i].E0,          // (11)
	   dE,                            // (12)
	   cluster->stars[i].dt);         // (13)
  }
}

/*************************************************/
void free_memory(CLUSTER *cluster)
{
  free(cluster->stars);
  free(cluster);
}

/*************************************************/
double wtime()
{
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return (double)((tv.tv_usec + (tv.tv_sec*1.0e6))/1.e6);
}

/*************************************************/
int main(int argc, char** argv)
{
  CLUSTER *cluster;
  INPUT params;
  
  get_args(argc, argv, &params);   
  readdata(&cluster);
  integrate(cluster, &params);
  free_memory(cluster);
}

