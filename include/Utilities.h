#ifndef UTILITIES_H
#define UTILITIES_H

// Definitions
#include "myLib.h"

// Timer
typedef struct Timer Timer;
struct Timer{
	int flag;
	time_t start, stop;
	double elapsed;
};

// Functions
void printComplex(complex double z);
double roundn(double x, int n);
void ticTimer(Timer *T);
void tocTimer(Timer *T);
void progressBar(int n, int N);
double deg2rad(double theta);
double rad2deg(double theta);

#endif