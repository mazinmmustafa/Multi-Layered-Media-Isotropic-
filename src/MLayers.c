//
#include "../include/MLayers.h"

void setConfig(Config *myConfig){
	int N=myConfig->N;
	assert(N>0);
	for (int i=0; i<N-1; i++){
		assert(myConfig->Layers[i].z_d==myConfig->Layers[i+1].z_u);
	}
	double z_u, z_d;
	complex double mu, eps; 
	double k0=myConfig->k0;
	for (int i=0; i<N; i++){
		z_u = myConfig->Layers[i].z_u;
		z_d = myConfig->Layers[i].z_d;
		mu = myConfig->Layers[i].mu;
		eps = myConfig->Layers[i].eps;
		myConfig->Layers[i].d = z_u-z_d;
		assert(myConfig->Layers[i].d>=0.0);
		myConfig->Layers[i].k = k0*csqrt(mu*eps);
		myConfig->Layers[i].eta = eta0*csqrt(mu/eps);
	}
}

void unsetConfig(Config *myConfig){
	assert(myConfig->Layers!=NULL);
	free(myConfig->Layers);
	myConfig->Layers = NULL;
	free(myConfig);
	myConfig = NULL;
}

void saveConfig(Config *myConfig){
	int N=myConfig->N;
	FILE *file=fopen("log.txt", "w");
	assert(file!=NULL);
	fprintf(file, "N = %d\n", myConfig->N);
	fprintf(file, "k0 = %21.14E\n", myConfig->k0);
	fprintf(file, "lambda0 = %21.14E\n", myConfig->lambda0);
	fprintf(file, "freq = %21.14E\n", myConfig->freq);
	fprintf(file, "omega = %21.14E\n", myConfig->omega);
	fprintf(file, "Gamma+ = (%21.14E, %21.14E)\n",
		creal(myConfig->Gamma_u), cimag(myConfig->Gamma_u));
	fprintf(file, "Gamma- = (%21.14E, %21.14E)\n",
		creal(myConfig->Gamma_d), cimag(myConfig->Gamma_d));
	fprintf(file, "\n");
	for (int i=0; i<N; i++){
		fprintf(file, "Layer %d:\n", i+1);
		fprintf(file, "z+ = %21.14E:\n", myConfig->Layers[i].z_u);
		fprintf(file, "z- = %21.14E:\n", myConfig->Layers[i].z_d);
		fprintf(file, "mu = (%21.14E, %21.14E)\n",
		creal(myConfig->Layers[i].mu), cimag(myConfig->Layers[i].mu));
		fprintf(file, "eps = (%21.14E, %21.14E)\n",
		creal(myConfig->Layers[i].eps), cimag(myConfig->Layers[i].eps));
		fprintf(file, "k = (%21.14E, %21.14E)\n",
		creal(myConfig->Layers[i].k), cimag(myConfig->Layers[i].k));
		fprintf(file, "eta = (%21.14E, %21.14E)\n",
		creal(myConfig->Layers[i].eta), cimag(myConfig->Layers[i].eta));
		fprintf(file, "\n");
	}
	fclose(file);
}

Config* createPaulus(){
	// Paulus configuration
	double nm=1.0E-9;
	complex double j=csqrt(-1.0);
	int N=4;
	double lambda0=633.0*nm;
	double k0=2.0*pi/lambda0;
	double freq=c0/lambda0;
	double omega=2.0*pi*freq;
	complex double Gamma_u=0.0;
	complex double Gamma_d=0.0;
	// Set general parameters
	Config *myConfig=(Config*)malloc(sizeof(Config));
	myConfig->Layers = (Layer*)malloc(N*sizeof(Layer));
	myConfig->N = N;
	myConfig->lambda0 = lambda0;
	myConfig->k0 = k0;
	myConfig->freq = freq;
	myConfig->omega = omega;
	myConfig->Gamma_u = Gamma_u;
	myConfig->Gamma_d = Gamma_d;
	int n=0;
	// Layer 1
	myConfig->Layers[n].z_u = +1000.0*nm;
	myConfig->Layers[n].z_d = +500.0*nm;
	myConfig->Layers[n].mu = 1.0-j*0.0;
	myConfig->Layers[n].eps = 1.0-j*0.0;
	myConfig->Layers[n].sigma_s = 0.0-j*0.0;
	n++;
	// Layer 2
	myConfig->Layers[n].z_u = +500.0*nm;
	myConfig->Layers[n].z_d = +0.0*nm;
	myConfig->Layers[n].mu = 1.0-j*0.0;
	myConfig->Layers[n].eps = 2.0-j*0.0;
	myConfig->Layers[n].sigma_s = 0.0-j*0.0;
	n++;
	// Layer 3
	myConfig->Layers[n].z_u = +0.0*nm;
	myConfig->Layers[n].z_d = -500.0*nm;
	myConfig->Layers[n].mu = 1.0-j*0.0;
	myConfig->Layers[n].eps = 10.0-j*0.0;
	myConfig->Layers[n].sigma_s = 0.0-j*0.0;
	n++;
	// Layer 4
	myConfig->Layers[n].z_u = -500.0*nm;
	myConfig->Layers[n].z_d = -1000.0*nm;
	myConfig->Layers[n].mu = 1.0-j*0.0;
	myConfig->Layers[n].eps = 1.0-j*0.0;
	myConfig->Layers[n].sigma_s = 0.0-j*0.0;
	n++;
	//
	return myConfig;
}

Config* createGoldKretschmann(){
	// Gold Kretschmann configuration
	double nm=1.0E-9;
	complex double j=csqrt(-1.0);
	int N=3;
	double lambda0=633.0*nm;
	double k0=2.0*pi/lambda0;
	double freq=c0/lambda0;
	double omega=2.0*pi*freq;
	complex double Gamma_u=0.0;
	complex double Gamma_d=0.0;
	// Set general parameters
	Config *myConfig=(Config*)malloc(sizeof(Config));
	myConfig->Layers = (Layer*)malloc(N*sizeof(Layer));
	myConfig->N = N;
	myConfig->lambda0 = lambda0;
	myConfig->k0 = k0;
	myConfig->freq = freq;
	myConfig->omega = omega;
	myConfig->Gamma_u = Gamma_u;
	myConfig->Gamma_d = Gamma_d;
	int n=0;
	// Layer 1
	myConfig->Layers[n].z_u = +2000.0*nm;
	myConfig->Layers[n].z_d = +0.0*nm;
	myConfig->Layers[n].mu = 1.0-j*0.0;
	myConfig->Layers[n].eps = 2.3013-j*0.0;
	myConfig->Layers[n].sigma_s = 0.0-j*0.0;
	n++;
	// Layer 2
	myConfig->Layers[n].z_u = +0.0*nm;
	myConfig->Layers[n].z_d = -50.0*nm;
	myConfig->Layers[n].mu = 1.0-j*0.0;
	myConfig->Layers[n].eps = -11.753-j*1.2596;
	myConfig->Layers[n].sigma_s = 0.0-j*0.0;
	n++;
	// Layer 3
	myConfig->Layers[n].z_u = -50.0*nm;
	myConfig->Layers[n].z_d = -2000.0*nm;
	myConfig->Layers[n].mu = 1.0-j*0.0;
	myConfig->Layers[n].eps = 1.0-j*0.0;
	myConfig->Layers[n].sigma_s = 0.0-j*0.0;
	n++;
	//
	return myConfig;
}

Config* createChew(){
	// Chew configuration
	complex double j=csqrt(-1.0);
	int N=7;
	double lambda0=1.0;
	double k0=2.0*pi/lambda0;
	double freq=c0/lambda0;
	double omega=2.0*pi*freq;
	complex double Gamma_u=0.0;
	complex double Gamma_d=0.0;
	// Set general parameters
	Config *myConfig=(Config*)malloc(sizeof(Config));
	myConfig->Layers = (Layer*)malloc(N*sizeof(Layer));
	myConfig->N = N;
	myConfig->lambda0 = lambda0;
	myConfig->k0 = k0;
	myConfig->freq = freq;
	myConfig->omega = omega;
	myConfig->Gamma_u = Gamma_u;
	myConfig->Gamma_d = Gamma_d;
	int n=0;
	// Layer 1
	myConfig->Layers[n].z_u = +10.0;
	myConfig->Layers[n].z_d = +0.0;
	myConfig->Layers[n].mu = 1.0-j*0.0;
	myConfig->Layers[n].eps = 1.0-j*0.0;
	myConfig->Layers[n].sigma_s = 0.0-j*0.0;
	n++;
	// Layer 2
	myConfig->Layers[n].z_u = +0.0;
	myConfig->Layers[n].z_d = -0.2;
	myConfig->Layers[n].mu = 1.0-j*0.0;
	myConfig->Layers[n].eps = 2.6-j*0.0;
	myConfig->Layers[n].sigma_s = 0.0-j*0.0;
	n++;
	// Layer 3
	myConfig->Layers[n].z_u = -0.2;
	myConfig->Layers[n].z_d = -0.5;
	myConfig->Layers[n].mu = 3.2-j*0.0;
	myConfig->Layers[n].eps = 6.5-j*0.0;
	myConfig->Layers[n].sigma_s = 0.0-j*0.0;
	n++;
	// Layer 4
	myConfig->Layers[n].z_u = -0.5;
	myConfig->Layers[n].z_d = -1.0;
	myConfig->Layers[n].mu = 6.0-j*0.0;
	myConfig->Layers[n].eps = 4.2-j*0.0;
	myConfig->Layers[n].sigma_s = 0.0-j*0.0;
	n++;
	// Layer 5
	myConfig->Layers[n].z_u = -1.0;
	myConfig->Layers[n].z_d = -1.3;
	myConfig->Layers[n].mu = 3.2-j*0.0;
	myConfig->Layers[n].eps = 6.5-j*0.0;
	myConfig->Layers[n].sigma_s = 0.0-j*0.0;
	n++;
	// Layer 6
	myConfig->Layers[n].z_u = -1.3;
	myConfig->Layers[n].z_d = -1.5;
	myConfig->Layers[n].mu = 1.0-j*0.0;
	myConfig->Layers[n].eps = 2.6-j*0.0;
	myConfig->Layers[n].sigma_s = 0.0-j*0.0;
	n++;
	// Layer 7
	myConfig->Layers[n].z_u = -1.5;
	myConfig->Layers[n].z_d = -10.0;
	myConfig->Layers[n].mu = 1.0-j*0.0;
	myConfig->Layers[n].eps = 1.0-j*0.0;
	myConfig->Layers[n].sigma_s = 0.0-j*0.0;
	n++;
	//
	return myConfig;
}

Config* createFEKO(){
	// FEKO configuration
	complex double j=csqrt(-1.0);
	int N=4;
	double lambda0=1.0;
	double k0=2.0*pi/lambda0;
	double freq=c0/lambda0;
	double omega=2.0*pi*freq;
	complex double Gamma_u=0.0;
	complex double Gamma_d=0.0;
	// Set general parameters
	Config *myConfig=(Config*)malloc(sizeof(Config));
	myConfig->Layers = (Layer*)malloc(N*sizeof(Layer));
	myConfig->N = N;
	myConfig->lambda0 = lambda0;
	myConfig->k0 = k0;
	myConfig->freq = freq;
	myConfig->omega = omega;
	myConfig->Gamma_u = Gamma_u;
	myConfig->Gamma_d = Gamma_d;
	int n=0;
	// Layer 1
	myConfig->Layers[n].z_u = +10.0;
	myConfig->Layers[n].z_d = +0.0;
	myConfig->Layers[n].mu = 1.0-j*0.0;
	myConfig->Layers[n].eps = 1.0-j*0.0;
	myConfig->Layers[n].sigma_s = 0.0-j*0.0;
	n++;
	// Layer 2
	myConfig->Layers[n].z_u = +0.0;
	myConfig->Layers[n].z_d = -0.2;
	myConfig->Layers[n].mu = 3.2-j*0.0;
	myConfig->Layers[n].eps = 6.5-j*0.0;
	myConfig->Layers[n].sigma_s = 0.0-j*0.0;
	n++;
	// Layer 3
	myConfig->Layers[n].z_u = -0.2;
	myConfig->Layers[n].z_d = -0.5;
	myConfig->Layers[n].mu = 6.0-j*0.0;
	myConfig->Layers[n].eps = 4.2-j*0.0;
	myConfig->Layers[n].sigma_s = 0.0-j*0.0;
	n++;
	// Layer 4
	myConfig->Layers[n].z_u = -0.5;
	myConfig->Layers[n].z_d = -10.0;
	myConfig->Layers[n].mu = 1.0-j*0.0;
	myConfig->Layers[n].eps = 1.0-j*0.0;
	myConfig->Layers[n].sigma_s = 0.0-j*0.0;
	n++;
	//
	return myConfig;
}