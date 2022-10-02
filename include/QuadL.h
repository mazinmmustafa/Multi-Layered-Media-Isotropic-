#ifndef QUADL_H
#define QUADL_H

// Definitions 
#include "myLib.h"

// Functions
void QuadL(int N, double *x, double *w);
void printQuadLRule(int Nmax);
complex double QuadL_1D(complex double func(complex double, void*), 
	void *args, double a, double b, double tol);
	
#endif