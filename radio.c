#include <sys/types.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "radio.h"
#define BUFFERSIZE 300
#define LARGE_N 10000000
#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif


/*************************************************/
void get_args(int argc, char** argv, INPUT *params)
{
  // User defined parameters
  params->dt    = pow(2.0,-8.0);    
  params->dtadj = pow(2.0,-8.0);
  params->tend  = 1;    
  params->dtout = pow(2.0,-1.0);


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

  params->tout  = params->dtout ; 
  params->tadj  = params->dtadj ; 

}

/*************************************************/
void readdata(CLUSTER **cluster)
{
  STAR *stars = malloc(sizeof(STAR));
  float** data = malloc(7 * sizeof(float*));      // allocate the rows
  for (int i = 0; i < 7; ++i)
    data[i] = malloc(LARGE_N * sizeof(float));   // allocate the columns
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
  //  fprintf(stderr, " N = %10d \n",i);
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

  sort(*cluster);
  get_phi(*cluster);
  for (i=0; i<(*cluster)->N; i++){
    (*cluster)->stars[i].E0 += (*cluster)->stars[i].phi;
    (*cluster)->stars[i].E = (*cluster)->stars[i].E0;
    
    float dt = (2.0*M_PI/pow(-2.0*(*cluster)->stars[i].E0,1.5))/256.;
    float power = 8.0;      
    while(power < 1.0/dt)
      power*=2.0;
    dt = 1.0/power;
    (*cluster)->stars[i].dt = dt;
  }

  adjust(*cluster);

  for (i=0; i<(*cluster)->N; ++i)
    (*cluster)->stars[i].id = i+1;

  output(*cluster);
  free(data);
  
}

/*************************************************/
void adjust(CLUSTER *cluster)
{
  float vr2 = 0.0, vt2 = 0.0, vr2i, vt2i;
  sort(cluster);
  get_phi(cluster);

  // Energy
  cluster->M = 0.0;
  cluster->E = 0.0;
  cluster->K = 0.0;
  cluster->W = 0.0;
  cluster->J = 0.0;

  // Cumulative mass 
  cluster->stars[0].cmass = cluster->stars[0].mass;

  for (int i=0; i<cluster->N; ++i){
    vr2i = pow(cluster->stars[i].vr,2.0);
    vt2i = pow(cluster->stars[i].vt,2.0);

    vr2 += vr2i;
    vt2 += vt2i;

    cluster->M += cluster->stars[i].mass;
    cluster->W += 0.5*cluster->stars[i].mass*cluster->stars[i].phi;
    cluster->K += 0.5*cluster->stars[i].mass*(vr2i + vt2i);
    cluster->J += cluster->stars[i].mass*cluster->stars[i].vt*cluster->stars[i].r;
    cluster->stars[i].E = 0.5*(vr2i + vt2i) + cluster->stars[i].phi;
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

  fprintf(stderr, " t = %7.3f  M = %10.8f  K = %10.7f  W = %10.7f  E = %10.7f  J = %10.7f  <vr2> = %10.7f  <vt2> = %10.7f  rh = %10.7f  \n",
	  cluster->t, cluster->M, cluster->K, cluster->W,  cluster->W+ cluster->K, cluster->J, vr2/(float)cluster->N, vt2/(float)cluster->N, cluster->rh);
}

/*************************************************/
void sort(CLUSTER *cluster)
{
  // (Shell) sort stars in order of increasing r
  int n = cluster->N;
  int i,j,inc; 
  float v;
  STAR vstar;
  inc = 1;

  do {
    inc *= 3;
    inc++;
  } while (inc <= n); 
  do {
    inc /= 3;
    for (i=inc; i<n; i++) {
      v=cluster->stars[i].r;
      vstar=cluster->stars[i];
      j=i;
      while (cluster->stars[j-inc].r > v) {
	cluster->stars[j]=cluster->stars[j-inc];
	j -= inc;
	if (j < inc) break;
      }
      cluster->stars[j]=vstar; 
    }
  } while (inc >= 1);  
}

/*************************************************/
void get_phi(CLUSTER *cluster)
{
  for (int i=0; i<cluster->N; ++i){
    cluster->stars[i].phi = 0.0;
    for (int j=0; j<i; ++j)
      cluster->stars[i].phi -= cluster->stars[j].mass/cluster->stars[i].r;
    for (int j=i+1; j<cluster->N; ++j)
      cluster->stars[i].phi -= cluster->stars[j].mass/cluster->stars[j].r;
  }
}

/*************************************************/
void integrate(CLUSTER *cluster, INPUT params)
{
  while (cluster->t < params.tend)
    {
      //Action
      int N = (int)cluster->N;
      float *r = (float *)malloc(N*sizeof(float));
      float *vr = (float *)malloc(N*sizeof(float));
      float *rp = (float *)malloc(N*sizeof(float));
      float *cm = (float *)malloc(N*sizeof(float));
      float *J2 = (float *)malloc(N*sizeof(float));
      
      for (int i=0; i<cluster->N; ++i){
	r[i] = (float)cluster->stars[i].r;
	vr[i] = (float)cluster->stars[i].vr;
	rp[i] = (float)cluster->stars[i].rp;
	cm[i] = (float)cluster->stars[i].cmass;
	J2[i] = (float)cluster->stars[i].J2;
      }
      
      rk4(r, vr, J2, rp, cm, N, params.dt); // On GPU

      for (int i=0; i<cluster->N; ++i){
	cluster->stars[i].r = r[i];
	cluster->stars[i].vr = vr[i];
	cluster->stars[i].vt = sqrt(J2[i])/r[i];
	cluster->stars[i].rp =  rp[i];
	cluster->stars[i].cmass = cm[i];
	cluster->stars[i].J2 = J2[i];
      }

      free(r);
      free(vr);
      free(rp);
      free(cm);
      free(J2);
      
      cluster->t += params.dt;
      
      if (cluster->t >= params.tadj){
	params.tadj += params.dtadj;
	adjust(cluster);
      }
      
      if (cluster->t >= params.tout){
	params.tout += params.dtout;
	output(cluster);
      }
    }

  
}


/*************************************************/
void output(CLUSTER *cluster)
{
  for (int i=0;i<cluster->N;i++){
    printf("%10.3e %7d %10.3e   %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e \n",
	   cluster->t,                    // (1)
	   cluster->stars[i].id,          // (2)
	   cluster->stars[i].mass,        // (3)
	   cluster->stars[i].r,           // (4)
	   cluster->stars[i].vr,          // (5)
	   cluster->stars[i].phi,         // (6)
	   cluster->stars[i].E,           // (7)
	   cluster->stars[i].J2,          // (8)
	   cluster->stars[i].E,           // (9)
	   cluster->stars[i].vt,          // (10)
	   cluster->stars[i].dt,          // (11)
	   cluster->stars[i].cmass,       // (12)
	   cluster->stars[i].E0,          // (13)
	   //	   cluster->stars[i].E);          // (13)
	   (cluster->stars[i].E - cluster->stars[i].E0)/cluster->stars[i].E0);    //  (14)

  }
}


/*************************************************/
void free_memory(CLUSTER *cluster)
{
  free(cluster->stars);
  free(cluster);
}
/*************************************************/
int main(int argc, char** argv)
{
  CLUSTER *cluster;
  INPUT params;

  time_t t0, t1;
  t0 = time(NULL);
  get_args(argc, argv, &params);   
  readdata(&cluster);
  integrate(cluster, params);
  free_memory(cluster);
  t1 = time(NULL);
  fprintf(stderr, " WALL TIME = %ld\n",(long)(t1-t0));
}

