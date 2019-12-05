#include "gray.h"
#include <cstdlib>
#include <cmath>
#include <stdio.h>

//relevent constants, probably are some more that we need
#define DT_DUMP (-64)
#define N_R     264 // large enough to hold our grids
#define N_RS    8   // number of cells to skip near the outer boundary
#define SINGLE 1  //i've turned everything into floats so this is unimportant
#define EPSILON 1e-32
#define CONST_c 2.9979e10
#define CONST_mp_me 1836.15
#define CONST_me 9.109e-28
#define R_SCHW 2 

//DEFINING THE STRUCTURES

//defining the Coord and Field Structs:
typedef struct {
  float t, r, theta, phi;
  float kr, ktheta;
  float bimpact;   // impact parameter defined as L / E, constant
  float padding;
} State; 


typedef struct {
  // size_t i, j, k;
  // double x1, x2, x3;
  // double v1, v2, v3;
  // double gcon[4][4];
  // double gcov[4][4];
  // double gdet;
  // double ck[4];
  // double dxdxp[4][4];
  float theta;
  float gcov [4][4];
  float dxdxp[4][4];
} Coord;

//output parameters
typedef struct {
  float r, theta, phi;
  float b, ne, te, br, btheta, bphi;
  float ur;
  float uphi;
  float beta;


} outparams;

typedef struct {
  float rho;
  float u;
  float jet1; // = -u_t
  float jet2; // = -T^r_t / (rho u^r)
  float ucont;
  float v1, v2, v3;
  float B1, B2, B3;
} Field;



//defining the Const struct

typedef struct {
  // Parameters for geodesic integration
  
  float a_spin;    // dimensionless spin j/mc
  float dt_scale;  // typical step size
  float epsilon;   // stop photon
  float tolerance; // if xi+1 > tolerance, fall back to forward Euler

  // Parameters for radiative transfer
  float   m_BH, ne_rho, threshold, tgas_max, Ti_Te_d, Ti_Te_f;

  Coord *coord;
  Field *field;
  size_t n_rx, n_r, n_theta, n_phi;
  float   Gamma;
  float   r[N_R];
} Const;



static Const c; //making c a global variable so that all functions can see it



//DEFINING THE FUNCTIONS

//load_coord, reads in a Const structure (c) and the usgdump2d file, these 
//parameters are the same across all sims, not dependent on time
//returns a Coord struct

//declaring the functions

//Coord *harm::load_coord(Const &c, const char *name);
Coord*  load_coord(char name[]);

Field* load_field(char name[]);

outparams transform(State s);

//Coord *harm::load_coord(Const &c, const char *name)
Coord*  load_coord(char name[])
{

  printf("%s", "testing print in load_coord");
 
 
  FILE *file = fopen(name, "r");
  if(!file)
    printf("%s", "error, fail to open file");
  //  error("ERROR: fail to open file \"%s\".\n", name);

  size_t count, n_r, n_theta,n_phi = 0;
  //size_t count, n_r, n_theta, n_phi;
  fseek(file, 12, SEEK_CUR);
  if(fread(&n_r,     sizeof(size_t), 1, file) != 1 ||
     fread(&n_theta, sizeof(size_t), 1, file) != 1 ||
     fread(&n_phi,   sizeof(size_t), 1, file) != 1)
    printf("%s", "Error: fail to read grid dimensions\n");

//fill global structure w/ grid info
  c.n_r = n_r;
  c.n_theta = n_theta;
  c.n_phi = n_phi;

  printf("%s", "loaded grid dimensions");

  c.n_rx = c.n_r - N_RS - 1; // use power-law extrapolation for r > c.r[c.n_rx]
  count  = c.n_r * c.n_theta;

  double Gamma, a_spin = 0; // a_spin is uninitialized if Gamma cannot be read
  fseek(file, 56, SEEK_CUR);
  if(fread(&Gamma,    sizeof(double), 1, file) != 1 ||
     fread(&a_spin,   sizeof(double), 1, file) != 1)
    printf("%s", "ERROR: fail to read Gamma or spin\n");
  //error("ERROR: fail to read Gamma or spin\n");
  fseek(file, 52, SEEK_CUR);
  c.Gamma  = Gamma;
  c.a_spin = a_spin;

  Coord *host;
  if(!(host = (Coord *)malloc(sizeof(Coord) * count)))
    printf("%s", "ERROR: fail to allocate host memory\n");
    //error("ERROR: fail to allocate host memory\n");
  for(size_t i = 0; i < count; ++i) {
    double in[16];

    fseek(file, 4 + 3 * sizeof(size_t) + 3 * sizeof(double), SEEK_CUR);

    if(fread(in, sizeof(double), 2, file) != 2)
      printf("%s", "ERROR: fail to read grid coordinates\n");
    //printf("%s", "ERROR: fail to read grid coordinates\n");
      //error("ERROR: fail to read grid coordinates\n");
      //error("ERROR: fail to read grid coordinates\n");
    if(i < c.n_r && i < N_R) c.r[i] = in[0];
    host[i].theta = in[1];
    //printf("%f\n", host[i].theta);
    fseek(file, 17 * sizeof(double), SEEK_CUR);

    if(fread(in, sizeof(double), 16, file) != 16)
      printf("%s", "ERROR: fail to read metric\n");
    //error("ERROR: fail to read metric\n");
    for(size_t j = 0; j < 16; ++j)
      (&(host[i].gcov[0][0]))[j] = in[j];

    fseek(file, 5 * sizeof(double), SEEK_CUR);

    if(fread(in, sizeof(double), 16, file) != 16)
      printf("%s", "ERROR: fail to read dxdxp\n");
    //error("ERROR: fail to read dxdxp\n");
    for(size_t j = 0; j < 16; ++j)
      (&(host[i].dxdxp[0][0]))[j] = in[j];

    fseek(file, 52, SEEK_CUR);
  }
  fclose(file);

  printf("Data size = %zu x %zu x %zu\n"
        "Gamma = %g, spin parameter a = %g, rmin = %g, rmax = %g\n",
        c.n_r, c.n_theta, c.n_phi, c.Gamma, c.a_spin, c.r[0], c.r[c.n_r-1]);
  
  return host;
}

//load field takes the Const struct, same accross timesteps in a sim, and
//the fieldline*.bin, with time-varying quantites (densities, rho, B)
//returns a Field structure



//loading field values, this is where density, temp, etc. are found after we've filled c coordinate struct w/ grid info
Field*  load_field(char name[])
{
  FILE *file = fopen(name, "r");
  if(!file)
    
    printf("ERROR: fail to open file \"%s\".\n", name);
    //error("ERROR: fail to open file \"%s\".\n", name);

  double time;
  size_t n1, n2, n3, count;
  if(fscanf(file, "%lf %zu %zu %zu", &time, &n1, &n2, &n3) != 4)
    printf("ERROR: fail to read time or grid size\n");
    //error("ERROR: fail to read time or grid size\n");
  while('\n' != fgetc(file));
  count = c.n_r * c.n_theta * c.n_phi; // num of grid points
  printf("%s", "n_r:\n");
  printf("%zu\n", c.n_r);
  printf("%s", "n_theta:\n");
  printf("%zu\n", c.n_theta);
  printf("%s", "n_phi:\n");
  printf("%zu\n", c.n_phi);
  if(n1 != c.n_r || n2 != c.n_theta || n3 != c.n_phi)
    printf("ERROR: inconsistent grid size\n");
    //error("ERROR: inconsistent grid size\n");

  Field *host;
  if(!(host = (Field *)malloc(sizeof(Field) * count)))
    printf("ERROR: fail to allocate host memory\n");
    //error("ERROR: fail to allocate host memory\n");
  if(fread(host, sizeof(Field), count, file) != count)
    printf("Error: fail to read data\n");
    //error("ERROR: fail to read data\n");
  fclose(file);

  float dmax = 0, Emax = 0, vmax = 0, Bmax = 0, Tmax = 0;
  for(size_t i = 0; i < count; ++i) {
    const float v = sqrt(host[i].v1 * host[i].v1 +
                        host[i].v2 * host[i].v2 +
                        host[i].v3 * host[i].v3);
    const float B = sqrt(host[i].B1 * host[i].B1 +
                        host[i].B2 * host[i].B2 +
                        host[i].B3 * host[i].B3);
    const float T = host[i].u / host[i].rho;

    if(dmax < host[i].rho) dmax = host[i].rho;
    if(Emax < host[i].u  ) Emax = host[i].u;
    if(vmax < v          ) vmax = v;
    if(Bmax < B          ) Bmax = B;
    if(Tmax < T          ) Tmax = T;
  }

  printf("Maximum density        = %g\n"
        "Maximum energy         = %g\n"
        "Maximum temperature    = %g\n"
        "Maximum speed          = %g\n"
        "Maximum magnetic field = %g\n",
        dmax, Emax, Tmax, vmax, Bmax);

  return host;
}



//Defining our transformation function, takes in state, returns state:

outparams transform(State s)

{
//defining some indices to access HARM data

outparams output = {0}; 
  output.theta = s.theta;
  output.r = s.r;
  output.phi = s.phi;

  //defining the quantities needed for metric calculations, new as of 3/8
  float a2 = c.a_spin * c.a_spin; 
  float r2 = s.r * s.r;           

  float sin_theta, cos_theta, c2, cs, s2;
  {
    //SINCOS(s.theta, &sin_theta, &cos_theta);//not totally sure if this line will work
    sin_theta=sinf(s.theta);
    cos_theta=cosf(s.theta);
    c2 = cos_theta * cos_theta;
    cs = cos_theta * sin_theta;
    s2 = sin_theta * sin_theta;
  }

  float g00, g11, g22, g33, g30, g33_s2, Dlt;//0 is t, 1 is r, 2 is theta, 3 is phi
  {
    //these are all the metric elements in BL coordinates
        float sum, tmp;
    sum    = r2 + a2;
    Dlt    = sum - R_SCHW * s.r; //equation (3) in Psaltis & Johannsen R_SCHW * s.r = 2Mr

    g22    = sum - a2 * s2; // = r2 + a2 * [cos(theta)]^2 = Sigma corresponds to gthetatheta
    g11    = g22 / Dlt; //grr

    tmp    = R_SCHW * s.r / g22; //   2Mr / (r^2 + a^2(1-sin^2(theta))) or just 2Mr/(r^2+a^2cos^2(theta)
    g00    = tmp - 1; //tt component of Kerr metic in BL coords
    g30    = -c.a_spin * tmp * s2;
    g33_s2 = sum - c.a_spin * g30;
    g33    = g33_s2 * s2;

    tmp    = 1 / (g33 * g00 - g30 * g30);
    
   } 

  float d001, d111, d212, d313, d020, d121, d222, d323, sum1, sum2;
    
  float a4 = a2*a2;
  float cot2 = c2/s2;
  float csc2 = 1/s2;
  float sum = a2+r2;



 float fg, Fg, fG, FG, fgh, Fgh, fGh, FGh, fgH, FgH, fGH, FGH;
 int  ij, Ij, iJ, IJ, ijk, Ijk, iJk, IJk, ijK, IjK, iJK, IJK;

{
    float F = 0.5;
    int  I = c.n_r - N_RS - 1, i = I; // assume c.n_r > N_RS + 1
    if(c.r[i] > s.r) {
      do I = i--; while(i && c.r[i] > s.r); // assume s.r >= c.r[0]
      F = (s.r - c.r[i]) / (c.r[I] - c.r[i]);
    } // else, constant extrapolate

    float G = (float)0.5;
    int  J = c.n_theta-1, j = J;
    if(s.theta < c.coord[i].theta * (1-F) +
                 c.coord[I].theta *    F)
      j = J = 0; // constant extrapolation
    else {
      float theta = c.coord[j*c.n_r + i].theta * (1-F)  +
                   c.coord[j*c.n_r + I].theta *    F;
      if(theta > s.theta) {
	float Theta;
        do {
          J = j--;
          Theta = theta;
          theta = c.coord[j*c.n_r + i].theta * (1-F) +
                  c.coord[j*c.n_r + I].theta *    F;
        } while(j && theta > s.theta);
        G = (s.theta - theta) / (Theta - theta);
      } // else, constant extrapolation
    }


    //float H = (s.phi - (M_PI/180) * c.j_obs) / (2*M_PI);
    float H = s.phi/(2*M_PI); // effectively setting j_obs to 0, don't see why it should matter here...
    //so H is basically what fraction of the way around phi is 

    H -= floor(H); //then we subtract off any bit that goes over 1, effectively making it wrap around 2pi

    H *=c.n_phi; //now multiply by how many grid cells in phi we have, which gives our grid value in phi
    int k = H, K = k == c.n_phi-1 ? 0 : k+1;
    H -= k;

    fg = (1-F) * (1-G);
    Fg =    F  * (1-G);
    fG = (1-F) *    G ;
    FG =    F  *    G ;

    ij = j * c.n_r + i;
    Ij = j * c.n_r + I;
    iJ = J * c.n_r + i;
    IJ = J * c.n_r + I;

    if(i > c.n_rx) {
      I = i = c.n_rx;
      F = 0.5;
    }

    //constants for interpolation
    fgh = (1-F) * (1-G) * (1-H);
    Fgh =    F  * (1-G) * (1-H);
    fGh = (1-F) *    G  * (1-H);
    FGh =    F  *    G  * (1-H);
    fgH = (1-F) * (1-G) *    H ;
    FgH =    F  * (1-G) *    H ;
    fGH = (1-F) *    G  *    H ;
    FGH =    F  *    G  *    H ;

    //grid indices given r, theta, phi
    ijk = (k * c.n_theta + j) * c.n_r + i;
    Ijk = (k * c.n_theta + j) * c.n_r + I;
    iJk = (k * c.n_theta + J) * c.n_r + i;
    IJk = (k * c.n_theta + J) * c.n_r + I;
    ijK = (K * c.n_theta + j) * c.n_r + i;
    IjK = (K * c.n_theta + j) * c.n_r + I;
    iJK = (K * c.n_theta + J) * c.n_r + i;
    IJK = (K * c.n_theta + J) * c.n_r + I;
}

  // Construct the four vectors u^\mu and b^\mu in modified KS coordinates
  float ut, ur, utheta, uphi, u; // u is energy
  float bt, br, btheta, bphi, b, ti_te;
  {
    const float gKSP00 = (fg*c.coord[ij].gcov[0][0]+Fg*c.coord[Ij].gcov[0][0]+
                         fG*c.coord[iJ].gcov[0][0]+FG*c.coord[IJ].gcov[0][0]);
    const float gKSP01 = (fg*c.coord[ij].gcov[0][1]+Fg*c.coord[Ij].gcov[0][1]+
                         fG*c.coord[iJ].gcov[0][1]+FG*c.coord[IJ].gcov[0][1]);
    const float gKSP02 = (fg*c.coord[ij].gcov[0][2]+Fg*c.coord[Ij].gcov[0][2]+
                         fG*c.coord[iJ].gcov[0][2]+FG*c.coord[IJ].gcov[0][2]);
    const float gKSP03 = (fg*c.coord[ij].gcov[0][3]+Fg*c.coord[Ij].gcov[0][3]+
                         fG*c.coord[iJ].gcov[0][3]+FG*c.coord[IJ].gcov[0][3]);
    const float gKSP11 = (fg*c.coord[ij].gcov[1][1]+Fg*c.coord[Ij].gcov[1][1]+
                         fG*c.coord[iJ].gcov[1][1]+FG*c.coord[IJ].gcov[1][1]);
    const float gKSP12 = (fg*c.coord[ij].gcov[1][2]+Fg*c.coord[Ij].gcov[1][2]+
                         fG*c.coord[iJ].gcov[1][2]+FG*c.coord[IJ].gcov[1][2]);
    const float gKSP13 = (fg*c.coord[ij].gcov[1][3]+Fg*c.coord[Ij].gcov[1][3]+
                         fG*c.coord[iJ].gcov[1][3]+FG*c.coord[IJ].gcov[1][3]);
    const float gKSP22 = (fg*c.coord[ij].gcov[2][2]+Fg*c.coord[Ij].gcov[2][2]+
                         fG*c.coord[iJ].gcov[2][2]+FG*c.coord[IJ].gcov[2][2]);
    const float gKSP23 = (fg*c.coord[ij].gcov[2][3]+Fg*c.coord[Ij].gcov[2][3]+
                         fG*c.coord[iJ].gcov[2][3]+FG*c.coord[IJ].gcov[2][3]);
    const float gKSP33 = (fg*c.coord[ij].gcov[3][3]+Fg*c.coord[Ij].gcov[3][3]+
                         fG*c.coord[iJ].gcov[3][3]+FG*c.coord[IJ].gcov[3][3]);

    ur     = (fgh * c.field[ijk].v1 + Fgh * c.field[Ijk].v1 +
              fGh * c.field[iJk].v1 + FGh * c.field[IJk].v1 +
              fgH * c.field[ijK].v1 + FgH * c.field[IjK].v1 +
              fGH * c.field[iJK].v1 + FGH * c.field[IJK].v1);
    
    

    utheta = (fgh * c.field[ijk].v2 + Fgh * c.field[Ijk].v2 +
              fGh * c.field[iJk].v2 + FGh * c.field[IJk].v2 +
              fgH * c.field[ijK].v2 + FgH * c.field[IjK].v2 +
              fGH * c.field[iJK].v2 + FGH * c.field[IJK].v2);
    uphi   = (fgh * c.field[ijk].v3 + Fgh * c.field[Ijk].v3 +
              fGh * c.field[iJk].v3 + FGh * c.field[IJk].v3 +
              fgH * c.field[ijK].v3 + FgH * c.field[IjK].v3 +
              fGH * c.field[iJK].v3 + FGH * c.field[IJK].v3);
    u      = (fgh * c.field[ijk].u  + Fgh * c.field[Ijk].u  +
              fGh * c.field[iJk].u  + FGh * c.field[IJk].u  +
              fgH * c.field[ijK].u  + FgH * c.field[IjK].u  +
              fGH * c.field[iJK].u  + FGH * c.field[IJK].u );
    br     = (fgh * c.field[ijk].B1 + Fgh * c.field[Ijk].B1 +
              fGh * c.field[iJk].B1 + FGh * c.field[IJk].B1 +
              fgH * c.field[ijK].B1 + FgH * c.field[IjK].B1 +
              fGH * c.field[iJK].B1 + FGH * c.field[IJK].B1);
    btheta = (fgh * c.field[ijk].B2 + Fgh * c.field[Ijk].B2 +
              fGh * c.field[iJk].B2 + FGh * c.field[IJk].B2 +
              fgH * c.field[ijK].B2 + FgH * c.field[IjK].B2 +
              fGH * c.field[iJK].B2 + FGH * c.field[IJK].B2);
    bphi   = (fgh * c.field[ijk].B3 + Fgh * c.field[Ijk].B3 +
              fGh * c.field[iJk].B3 + FGh * c.field[IJk].B3 +
              fgH * c.field[ijK].B3 + FgH * c.field[IjK].B3 +
              fGH * c.field[iJK].B3 + FGH * c.field[IJK].B3);
  
  
 // Vector u
    ut     = 1 / sqrt(EPSILON
                      -(gKSP00                   +
                        gKSP11 * ur     * ur     +
                        gKSP22 * utheta * utheta +
                        gKSP33 * uphi   * uphi   +
                        2 * (gKSP01 * ur              +
                             gKSP02 * utheta          +
                             gKSP03 * uphi            +
                             gKSP12 * ur     * utheta +
                             gKSP13 * ur     * uphi   +
                             gKSP23 * utheta * uphi)));
    ur     *= ut;
    utheta *= ut;
    uphi   *= ut;

    // Vector B
    bt     = (br     * (gKSP01 * ut     + gKSP11 * ur    +
                        gKSP12 * utheta + gKSP13 * uphi) +
              btheta * (gKSP02 * ut     + gKSP12 * ur    +
                        gKSP22 * utheta + gKSP23 * uphi) +
              bphi   * (gKSP03 * ut     + gKSP13 * ur    +
                        gKSP23 * uphi   + gKSP33 * uphi));
    br     = (br     + bt * ur    ) / ut;
    btheta = (btheta + bt * utheta) / ut;
    bphi   = (bphi   + bt * uphi  ) / ut;

    const float bb = (bt     * (gKSP00 * bt     + gKSP01 * br    +
                               gKSP02 * btheta + gKSP03 * bphi) +
                     br     * (gKSP01 * bt     + gKSP11 * br    +
                               gKSP12 * btheta + gKSP13 * bphi) +
                     btheta * (gKSP02 * bt     + gKSP12 * br    +
                               gKSP22 * btheta + gKSP23 * bphi) +
                     bphi   * (gKSP03 * bt     + gKSP13 * br    +
                               gKSP23 * btheta + gKSP33 * bphi));
    const float ibeta = bb / (2 * (c.Gamma-1) * u + EPSILON);
    ti_te = (ibeta > c.threshold) ? c.Ti_Te_f : (1 + c.Ti_Te_d * R_SCHW / s.r);
    b = sqrt(bb);
    output.beta = 1.0/ibeta;
    
  } 

 float rho, tgas;
  {
    const float Gamma = (1 + (ti_te+1) / (ti_te+2) / 1.5 + c.Gamma) / 2;

    rho  = (fgh * c.field[ijk].rho + Fgh * c.field[Ijk].rho +
            fGh * c.field[iJk].rho + FGh * c.field[IJk].rho +
            fgH * c.field[ijK].rho + FgH * c.field[IjK].rho +
            fGH * c.field[iJK].rho + FGH * c.field[IJK].rho);
    tgas = (Gamma - 1) * u / (rho + EPSILON);

  } 

 // Skip cell if tgas is above the threshold
outparams d = {0}; //define this null state to return for a cell if tgas is too high
//if(tgas > c.tgas_max) return d; // let's try commenting this out, see what happens 

 // Transform vector u and b from KSP to KS coordinates
  {
    const float dxdxp00=(fg*c.coord[ij].dxdxp[0][0]+Fg*c.coord[Ij].dxdxp[0][0]+
                        fG*c.coord[iJ].dxdxp[0][0]+FG*c.coord[IJ].dxdxp[0][0]);
    const float dxdxp11=(fg*c.coord[ij].dxdxp[1][1]+Fg*c.coord[Ij].dxdxp[1][1]+
                        fG*c.coord[iJ].dxdxp[1][1]+FG*c.coord[IJ].dxdxp[1][1]);
    const float dxdxp12=(fg*c.coord[ij].dxdxp[1][2]+Fg*c.coord[Ij].dxdxp[1][2]+
                        fG*c.coord[iJ].dxdxp[1][2]+FG*c.coord[IJ].dxdxp[1][2]);
    const float dxdxp21=(fg*c.coord[ij].dxdxp[2][1]+Fg*c.coord[Ij].dxdxp[2][1]+
                        fG*c.coord[iJ].dxdxp[2][1]+FG*c.coord[IJ].dxdxp[2][1]);
    const float dxdxp22=(fg*c.coord[ij].dxdxp[2][2]+Fg*c.coord[Ij].dxdxp[2][2]+
                        fG*c.coord[iJ].dxdxp[2][2]+FG*c.coord[IJ].dxdxp[2][2]);
    const float dxdxp33=(fg*c.coord[ij].dxdxp[3][3]+Fg*c.coord[Ij].dxdxp[3][3]+
                        fG*c.coord[iJ].dxdxp[3][3]+FG*c.coord[IJ].dxdxp[3][3]);
    

    float temp1, temp2;

    temp1  = ur;
    temp2  = utheta;
    ut    *= dxdxp00;
    ur     = (dxdxp11 * temp1 + dxdxp12 * temp2);
    utheta = (dxdxp21 * temp1 + dxdxp22 * temp2);
    uphi  *= dxdxp33;

    temp1  = br;
    temp2  = btheta;
    bt    *= dxdxp00;
    br     = (dxdxp11 * temp1 + dxdxp12 * temp2);
    btheta = (dxdxp21 * temp1 + dxdxp22 * temp2);
    bphi  *= dxdxp33;
  } 

 
  // Transform vector u and b from KS to BL coordinates
  {
    float temp0 = -R_SCHW * s.r / Dlt; // Note that s.r and Dlt are
    float temp3 = -c.a_spin     / Dlt; // evaluated at photon position

    ut   += ur * temp0;
    uphi += ur * temp3;

    bt   += br * temp0;
    bphi += br * temp3;
  } // 13 FLOP



//rewriting the above block with only the necessary quantities
//essentially just computing the cgs values of b, ne, and te. We are now in BL coordinates 
  float ne, te;
  {
  
    b *= sqrt(c.ne_rho) *
         (CONST_c * sqrt(4 * M_PI * (CONST_mp_me + 1) * CONST_me));

    //here i'm adding the transformation to cgs of btheta, br, and bphi
    //
    br *= sqrt(c.ne_rho) *
         (CONST_c * sqrt(4 * M_PI * (CONST_mp_me + 1) * CONST_me));
    btheta *= sqrt(c.ne_rho) *
         (CONST_c * sqrt(4 * M_PI * (CONST_mp_me + 1) * CONST_me));
    bphi *= sqrt(c.ne_rho) *
         (CONST_c * sqrt(4 * M_PI * (CONST_mp_me + 1) * CONST_me));
    //end of new stuff that I added

    ne = c.ne_rho * rho;
    te = ti_te < 0 ? -ti_te : tgas * CONST_mp_me / (ti_te+1);
  } 




//now that we've calculated these values the same way we did in rhs.h, let's output them

    output.b = b;
    output.ne = ne;
    output.te = te;
    output.ur = ur;
    output.uphi = uphi;
    output.br = br;
    output.btheta = btheta;
    output.bphi = bphi;
    return output;
    }

//maybe write a new function that does this and takes state as an input rather than doing it within the transform function

int main(void)
{
  //let's define our constants necessary for using img and load_field:
  //note that these are part of a GLOBAL structure
{
  //setting normalization quantities to 1
  c.a_spin    = 0.999;
  c.dt_scale  = 1.0 / 32;
  c.epsilon   = 1e-3;
  c.tolerance = 1e-1; 
  c.m_BH      = 1; // in unit of solar mass
  c.ne_rho    = 1;
  c.threshold = 5; // = 1 / beta_threshold
  c.tgas_max  = 1;
  c.Ti_Te_d   = 3;
  c.Ti_Te_f   = 3;
  c.coord = NULL;
  c.field = NULL;
}
  
  char fieldname[] = "testdata/fieldline3590.bin";
    
  char coordname[] = "testdata/usgdump2d";
  

 
  
    printf("%s", "output from load_coord:\n");
    c.coord = load_coord(coordname);
    
    printf("%s", "output from load_field:\n");
    c.field = load_field(fieldname);
 

    // so we want to modify this so that we specify what x,y,z we want, convert this to r,theta, phi, and then save into array given by x-y-z

    
    char outfile[] = "CSV_out/3D_test_out.txt";
    FILE *fp;
    fp = fopen(outfile, "w");
    


    // just make it small to begin with since we'll be saving a 3d array
    int array_size = 100.0;
    float max = 100.0; //in 2GM/c^2
    float grid_spacing = max / array_size;


//initializing instate, the values of which we will modify at every grid point to interpolate the spherical polar coords to rectangular 3d grid
State instate = {0};
    {
      instate.r = 0; 
      instate.theta = 0;
      instate.phi = 0; 
    } 
   
   float outarray[array_size][array_size][array_size];
   
   //first line has info on physical length of grid
   fprintf(fp, "%s", "physical length of grid (which is cubic): ");
   fprintf(fp, "%f", max);
   fprintf(fp, "%s", "\n");

   //second line has info on size of array

   fprintf(fp, "%s", "size of one side of array: ");
   fprintf(fp, "%d", array_size);
   fprintf(fp, "%s", "\n");




   fprintf(fp, "%s", "i , j , k,  b , ne , te \n");
    for(int i=0; i< array_size; i++)
    {
      float x = i*grid_spacing - max/2.;
      for(int j=0; j<array_size; j++)
      {
        float y = j*grid_spacing - max/2.;
        for (int k=0; k<array_size; k++)
        {
          float z = k*grid_spacing - max/2.;
          // so now we have the x,y,z that we want to find our values at.  Now we convert this to r theta phi for the instate
          float radius = sqrt(x*x+y*y+z*z);
          float phi = 2*atan(y/x);
          // make phi go from 0 to 2pi?
          //phi += M_PI;

          float theta = acos(z/radius);

          instate.r = radius;
          instate.theta=theta;
          instate.phi=phi;

          outparams out = transform(instate);
          if(radius<1)//zeroing out grid spaces within event horizon
	      {
		        out.b = 0;
		        out.ne = 0;
		        out.te = 0;
	      }

         fprintf(fp, "%d , %d , %d , %f , %f , %f \n", i, j, k, out.b, out.ne,out.te);
        }
      }
    }
    fclose(fp);

  return 0;
}
