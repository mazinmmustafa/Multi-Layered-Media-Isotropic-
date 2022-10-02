#ifndef BESSEL_H
#define BESSEL_H

// Definitions
#include "myLib.h"

// Functions
complex double besselj(double n, complex double z);
complex double bessely(double n, complex double z);
complex double besseli(double n, complex double z);
complex double besselk(double n, complex double z);
complex double besselh(int k, double n, complex double z);

#endif