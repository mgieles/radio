typedef struct input{
  float tend;
  float dt;
  float dtout;
  float dtadj;
  float tout;
  float tadj;
}INPUT;

typedef struct star{
  int id;
  float mass;
  float cmass;
  float phi, dt;
  float r, vr, vt, E0, E, J2;
  float rp;
  float rnew;
} STAR;

typedef struct cluster{
  int N;
  float M, E, J, K, W;
  float t ;
  float rh;
  float *rp;
  STAR *stars;
} CLUSTER;

float myrand();
void get_args(int argc, char* argv[], INPUT *params);
void readdata(CLUSTER **cluster);
void sort(CLUSTER *cluster);
void adjust(CLUSTER *cluster);
void get_phi(CLUSTER *cluster);
void integrate(CLUSTER *cluster, INPUT params);
void sort(CLUSTER *cluster);
void output(CLUSTER *cluster);
void free_memory(CLUSTER *cluster);
void rk4(float *r, float *vr, float *J2, float *rp, float *cm, int N, float dt); 



