#ifndef ENGINE_H
#define ENGINE_H

// Definitions
#include "myLib.h"
#include "MLayers.h"
#include "Bessel.h"
#include "QuadL.h"

// Struct GfuncArgs
typedef struct GfuncArgs GfuncArgs;
struct GfuncArgs{
	Config *myConfig;
	double z, z_, rho;
	double a, b;
};

// Functions
void spectralParameters(Config *myConfig, int n, 
	char pol, complex double *Z, complex double *kz, 
	complex double *Theta, complex double k_rho);
complex double Refl(Config *myConfig, int n, char pol, 
	char dir, complex double k_rho);
int findLayer(Config *myConfig, double z);
void TLGF(Config *myConfig, char pol, char sor, 
	double z, double z_, complex double k_rho, 
	complex double *V, complex double *I);
void TLGFr(Config *myConfig, char pol, char sor, 
	double z, double z_, complex double k_rho, 
	complex double *V, complex double *I);
// GEJ
complex double GEJxx(Config *myConfig, double x, double x_, 
	double y,  double y_, double z, double z_, double tol);
complex double GEJxy(Config *myConfig, double x, double x_, 
	double y,  double y_, double z, double z_, double tol);
complex double GEJxz(Config *myConfig, double x, double x_, 
	double y,  double y_, double z, double z_, double tol);
complex double GEJyx(Config *myConfig, double x, double x_, 
	double y,  double y_, double z, double z_, double tol);
complex double GEJyy(Config *myConfig, double x, double x_, 
	double y,  double y_, double z, double z_, double tol);
complex double GEJyz(Config *myConfig, double x, double x_, 
	double y,  double y_, double z, double z_, double tol);
complex double GEJzx(Config *myConfig, double x, double x_, 
	double y,  double y_, double z, double z_, double tol);
complex double GEJzy(Config *myConfig, double x, double x_, 
	double y,  double y_, double z, double z_, double tol);
complex double GEJzz(Config *myConfig, double x, double x_, 
	double y,  double y_, double z, double z_, double tol);
// GEM
complex double GEMxx(Config *myConfig, double x, double x_, 
	double y,  double y_, double z, double z_, double tol);
complex double GEMxy(Config *myConfig, double x, double x_, 
	double y,  double y_, double z, double z_, double tol);
complex double GEMxz(Config *myConfig, double x, double x_, 
	double y,  double y_, double z, double z_, double tol);
complex double GEMyx(Config *myConfig, double x, double x_, 
	double y,  double y_, double z, double z_, double tol);
complex double GEMyy(Config *myConfig, double x, double x_, 
	double y,  double y_, double z, double z_, double tol);
complex double GEMyz(Config *myConfig, double x, double x_, 
	double y,  double y_, double z, double z_, double tol);
complex double GEMzx(Config *myConfig, double x, double x_, 
	double y,  double y_, double z, double z_, double tol);
complex double GEMzy(Config *myConfig, double x, double x_, 
	double y,  double y_, double z, double z_, double tol);
complex double GEMzz(Config *myConfig, double x, double x_, 
	double y,  double y_, double z, double z_, double tol);
// GHJ
complex double GHJxx(Config *myConfig, double x, double x_, 
	double y,  double y_, double z, double z_, double tol);
complex double GHJxy(Config *myConfig, double x, double x_, 
	double y,  double y_, double z, double z_, double tol);
complex double GHJxz(Config *myConfig, double x, double x_, 
	double y,  double y_, double z, double z_, double tol);
complex double GHJyx(Config *myConfig, double x, double x_, 
	double y,  double y_, double z, double z_, double tol);
complex double GHJyy(Config *myConfig, double x, double x_, 
	double y,  double y_, double z, double z_, double tol);
complex double GHJyz(Config *myConfig, double x, double x_, 
	double y,  double y_, double z, double z_, double tol);
complex double GHJzx(Config *myConfig, double x, double x_, 
	double y,  double y_, double z, double z_, double tol);
complex double GHJzy(Config *myConfig, double x, double x_, 
	double y,  double y_, double z, double z_, double tol);
complex double GHJzz(Config *myConfig, double x, double x_, 
	double y,  double y_, double z, double z_, double tol);
// GHM
complex double GHMxx(Config *myConfig, double x, double x_, 
	double y,  double y_, double z, double z_, double tol);
complex double GHMxy(Config *myConfig, double x, double x_, 
	double y,  double y_, double z, double z_, double tol);
complex double GHMxz(Config *myConfig, double x, double x_, 
	double y,  double y_, double z, double z_, double tol);
complex double GHMyx(Config *myConfig, double x, double x_, 
	double y,  double y_, double z, double z_, double tol);
complex double GHMyy(Config *myConfig, double x, double x_, 
	double y,  double y_, double z, double z_, double tol);
complex double GHMyz(Config *myConfig, double x, double x_, 
	double y,  double y_, double z, double z_, double tol);
complex double GHMzx(Config *myConfig, double x, double x_, 
	double y,  double y_, double z, double z_, double tol);
complex double GHMzy(Config *myConfig, double x, double x_, 
	double y,  double y_, double z, double z_, double tol);
complex double GHMzz(Config *myConfig, double x, double x_, 
	double y,  double y_, double z, double z_, double tol);
	
#endif