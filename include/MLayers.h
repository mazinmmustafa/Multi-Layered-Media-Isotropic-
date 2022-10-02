#ifndef MLAYERS_H
#define MLAYERS_H

// Definitions
#include "myLib.h"

// Layer
typedef struct Layer Layer;
struct Layer{
	double z_u, z_d, d;
	complex double mu, eps, k, eta, sigma_s;
};

// Configuration
typedef struct Config Config;
struct Config{
	int N;
	double lambda0, k0, freq, omega;
	complex double Gamma_u, Gamma_d;
	Layer *Layers;
};

// Functions
Config* createPaulus();
Config* createGoldKretschmann();
Config* createChew();
Config* createFEKO();
void setConfig(Config *myConfig);
void unsetConfig(Config *myConfig);
void saveConfig(Config *myConfig);

#endif