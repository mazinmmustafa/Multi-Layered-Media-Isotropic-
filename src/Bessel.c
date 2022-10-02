//
#include "../include/Bessel.h"

// Fortran Functions
extern void zbesj_f77_(double *fnu, double *zr, 
	double *zi, double *cyr, double *cyi);
	
extern void zbesy_f77_(double *fnu, double *zr, 
	double *zi, double *cyr, double *cyi);
	
extern void zbesi_f77_(double *fnu, double *zr, 
	double *zi, double *cyr, double *cyi);
	
extern void zbesk_f77_(double *fnu, double *zr, 
	double *zi, double *cyr, double *cyi);
	
extern void zbesh_f77_(int *k, double *fnu, double *zr, 
	double *zi, double *cyr, double *cyi);
	
complex double besselj(double n, complex double z){
	double zr=creal(z);
	double zi=cimag(z);
	double cyr, cyi;
	complex double j=csqrt(-1.0);
	zbesj_f77_(&n, &zr, &zi, &cyr, &cyi);
	return cyr+j*cyi;
}

complex double bessely(double n, complex double z){
	double zr=creal(z);
	double zi=cimag(z);
	double cyr, cyi;
	complex double j=csqrt(-1.0);
	zbesy_f77_(&n, &zr, &zi, &cyr, &cyi);
	return cyr+j*cyi;
}

complex double besseli(double n, complex double z){
	double zr=creal(z);
	double zi=cimag(z);
	double cyr, cyi;
	complex double j=csqrt(-1.0);
	zbesi_f77_(&n, &zr, &zi, &cyr, &cyi);
	return cyr+j*cyi;
}

complex double besselk(double n, complex double z){
	double zr=creal(z);
	double zi=cimag(z);
	double cyr, cyi;
	complex double j=csqrt(-1.0);
	zbesk_f77_(&n, &zr, &zi, &cyr, &cyi);
	return cyr+j*cyi;
}

complex double besselh(int k, double n, complex double z){
	assert(k==1||k==2);
	double zr=creal(z);
	double zi=cimag(z);
	double cyr, cyi;
	complex double j=csqrt(-1.0);
	zbesh_f77_(&k, &n, &zr, &zi, &cyr, &cyi);
	return cyr+j*cyi;
}