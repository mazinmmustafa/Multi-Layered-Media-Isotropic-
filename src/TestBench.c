//
#include "../include/TestBench.h"

void TestPaulus(){
	Config *myConfig=createPaulus();
	setConfig(myConfig);
	saveConfig(myConfig);
	//
	double nm=1.0E-9;
	int Ns=1000;
	double z_min=-1000.0*nm;
	double z_max=+1000.0*nm;
	double phi=deg2rad(45.0);
	double rho=633.0*nm;
	double x_=0.0*nm;
	double y_=0.0*nm;
	double z_=+750.0*nm;
	double tol=1.0E-3;
	//
	double dz=(z_max-z_min)/(Ns-1.0);
	double x, y, z;
	x = rho*cos(phi);
	y = rho*sin(phi);
	complex double Gxx, Gxy, Gxz;
	complex double Gyx, Gyy, Gyz;
	complex double Gzx, Gzy, Gzz;
	Timer T;
	ticTimer(&T);
	FILE *file=fopen("Data/DGF/Data.dat", "w");
	assert(file!=NULL);
	for (int i=0; i<Ns; i++){
		z = z_min+i*dz;
		Gxx = GEJxx(myConfig, x, x_, y, y_, z, z_, tol);
		Gxy = GEJxy(myConfig, x, x_, y, y_, z, z_, tol);
		Gxz = GEJxz(myConfig, x, x_, y, y_, z, z_, tol);
		Gyx = GEJyx(myConfig, x, x_, y, y_, z, z_, tol);
		Gyy = GEJyy(myConfig, x, x_, y, y_, z, z_, tol);
		Gyz = GEJyz(myConfig, x, x_, y, y_, z, z_, tol);
		Gzx = GEJzx(myConfig, x, x_, y, y_, z, z_, tol);
		Gzy = GEJzy(myConfig, x, x_, y, y_, z, z_, tol);
		Gzz = GEJzz(myConfig, x, x_, y, y_, z, z_, tol);
		fprintf(file, "%21.14E ", z);
		fprintf(file, "%21.14E %21.14E %21.14E ", cabs(Gxx), cabs(Gxy), cabs(Gxz));
		fprintf(file, "%21.14E %21.14E %21.14E ", cabs(Gyx), cabs(Gyy), cabs(Gyz));
		fprintf(file, "%21.14E %21.14E %21.14E ", cabs(Gzx), cabs(Gzy), cabs(Gzz));
		fprintf(file, "\n");
		progressBar(i, Ns);
	}
	fclose(file);
	tocTimer(&T);
	//
	unsetConfig(myConfig);
}

void TestChew1(){
	Config *myConfig=createChew();
	setConfig(myConfig);
	saveConfig(myConfig);
	//
	int Ns=1000;
	double x_min=-3.0;
	double x_max=+3.0;
	double y=+1.0;
	double z=-0.3;
	double x_=+0.0;
	double y_=+0.0;
	double z_=-1.4;
	double theta0=deg2rad(20.0);
	double phi0=deg2rad(30.0);
	double tol=1.0E-3;
	//
	double Jx=sin(theta0)*cos(phi0);
	double Jy=sin(theta0)*sin(phi0);
	double Jz=cos(theta0);
	double dx=(x_max-x_min)/(Ns-1.0);
	double x;
	complex double Ex, Ey, Ez;
	Timer T;
	ticTimer(&T);
	FILE *file=fopen("Data/Field/Data1.dat", "w");
	assert(file!=NULL);
	for (int i=0; i<Ns; i++){
		
		x = x_min+i*dx;
		Ex = GEJxx(myConfig, x, x_, y, y_, z, z_, tol)*Jx+
		     GEJxy(myConfig, x, x_, y, y_, z, z_, tol)*Jy+
			 GEJxz(myConfig, x, x_, y, y_, z, z_, tol)*Jz;
		Ey = GEJyx(myConfig, x, x_, y, y_, z, z_, tol)*Jx+
		     GEJyy(myConfig, x, x_, y, y_, z, z_, tol)*Jy+
			 GEJyz(myConfig, x, x_, y, y_, z, z_, tol)*Jz;
		Ez = GEJzx(myConfig, x, x_, y, y_, z, z_, tol)*Jx+
		     GEJzy(myConfig, x, x_, y, y_, z, z_, tol)*Jy+
			 GEJzz(myConfig, x, x_, y, y_, z, z_, tol)*Jz;
		fprintf(file, "%21.14E ", x);
		fprintf(file, "%21.14E %21.14E %21.14E ", cabs(Ex), cabs(Ey), cabs(Ez));
		fprintf(file, "\n");
		progressBar(i, Ns);
	}
	fclose(file);
	tocTimer(&T);
	//
	unsetConfig(myConfig);
}

void TestChew2(){
	Config *myConfig=createChew();
	setConfig(myConfig);
	saveConfig(myConfig);
	//
	int Ns=1000;
	double x_min=-3.0;
	double x_max=+3.0;
	double y=+1.0;
	double z=-0.3;
	double x_=+0.0;
	double y_=+0.0;
	double z_=-1.4;
	double theta0=deg2rad(20.0);
	double phi0=deg2rad(30.0);
	double tol=1.0E-3;
	//
	double Mx=sin(theta0)*cos(phi0);
	double My=sin(theta0)*sin(phi0);
	double Mz=cos(theta0);
	double dx=(x_max-x_min)/(Ns-1.0);
	double x;
	complex double Ex, Ey, Ez;
	Timer T;
	ticTimer(&T);
	FILE *file=fopen("Data/Field/Data2.dat", "w");
	assert(file!=NULL);
	for (int i=0; i<Ns; i++){
		
		x = x_min+i*dx;
		Ex = GEMxx(myConfig, x, x_, y, y_, z, z_, tol)*Mx+
		     GEMxy(myConfig, x, x_, y, y_, z, z_, tol)*My+
			 GEMxz(myConfig, x, x_, y, y_, z, z_, tol)*Mz;
		Ey = GEMyx(myConfig, x, x_, y, y_, z, z_, tol)*Mx+
		     GEMyy(myConfig, x, x_, y, y_, z, z_, tol)*My+
			 GEMyz(myConfig, x, x_, y, y_, z, z_, tol)*Mz;
		Ez = GEMzx(myConfig, x, x_, y, y_, z, z_, tol)*Mx+
		     GEMzy(myConfig, x, x_, y, y_, z, z_, tol)*My+
			 GEMzz(myConfig, x, x_, y, y_, z, z_, tol)*Mz;
		fprintf(file, "%21.14E ", x);
		fprintf(file, "%21.14E %21.14E %21.14E ", cabs(Ex), cabs(Ey), cabs(Ez));
		fprintf(file, "\n");
		progressBar(i, Ns);
	}
	fclose(file);
	tocTimer(&T);
	//
	unsetConfig(myConfig);
}

void TestFEKO1(){
	Config *myConfig=createFEKO();
	setConfig(myConfig);
	saveConfig(myConfig);
	//
	int Ns=1000;
	double x_min=-3.0;
	double x_max=+3.0;
	double y=+1.0;
	double z=-0.3;
	double x_=+0.0;
	double y_=+0.0;
	double z_=-0.4;
	double theta0=deg2rad(20.0);
	double phi0=deg2rad(30.0);
	double tol=1.0E-3;
	//
	double Jx=sin(theta0)*cos(phi0);
	double Jy=sin(theta0)*sin(phi0);
	double Jz=cos(theta0);
	double Mx=sin(theta0)*cos(phi0);
	double My=sin(theta0)*sin(phi0);
	double Mz=cos(theta0);
	double dx=(x_max-x_min)/(Ns-1.0);
	double x;
	complex double Ex, Ey, Ez;
	// Test J Dipole
	printf("Computing J Dipole:\n");
	Timer T;
	ticTimer(&T);
	FILE *file=fopen("Data/Validation/HCut/Data1.dat", "w");
	assert(file!=NULL);
	for (int i=0; i<Ns; i++){
		x = x_min+i*dx;
		fprintf(file, "%21.14E ", x);
		Ex = GEJxx(myConfig, x, x_, y, y_, z, z_, tol)*Jx+
		     GEJxy(myConfig, x, x_, y, y_, z, z_, tol)*Jy+
			 GEJxz(myConfig, x, x_, y, y_, z, z_, tol)*Jz;
		Ey = GEJyx(myConfig, x, x_, y, y_, z, z_, tol)*Jx+
		     GEJyy(myConfig, x, x_, y, y_, z, z_, tol)*Jy+
			 GEJyz(myConfig, x, x_, y, y_, z, z_, tol)*Jz;
		Ez = GEJzx(myConfig, x, x_, y, y_, z, z_, tol)*Jx+
		     GEJzy(myConfig, x, x_, y, y_, z, z_, tol)*Jy+
			 GEJzz(myConfig, x, x_, y, y_, z, z_, tol)*Jz;
		fprintf(file, "%21.14E %21.14E %21.14E ", cabs(Ex), cabs(Ey), cabs(Ez));
		Ex = GHJxx(myConfig, x, x_, y, y_, z, z_, tol)*Jx+
		     GHJxy(myConfig, x, x_, y, y_, z, z_, tol)*Jy+
			 GHJxz(myConfig, x, x_, y, y_, z, z_, tol)*Jz;
		Ey = GHJyx(myConfig, x, x_, y, y_, z, z_, tol)*Jx+
		     GHJyy(myConfig, x, x_, y, y_, z, z_, tol)*Jy+
			 GHJyz(myConfig, x, x_, y, y_, z, z_, tol)*Jz;
		Ez = GHJzx(myConfig, x, x_, y, y_, z, z_, tol)*Jx+
		     GHJzy(myConfig, x, x_, y, y_, z, z_, tol)*Jy+
			 GHJzz(myConfig, x, x_, y, y_, z, z_, tol)*Jz;
		fprintf(file, "%21.14E %21.14E %21.14E ", cabs(Ex), cabs(Ey), cabs(Ez));
		fprintf(file, "\n");
		progressBar(i, Ns);
	}
	fclose(file);
	tocTimer(&T);
	// Test M Dipole
	printf("Computing M Dipole:\n");
	ticTimer(&T);
	file = fopen("Data/Validation/HCut/Data2.dat", "w");
	assert(file!=NULL);
	for (int i=0; i<Ns; i++){
		x = x_min+i*dx;
		fprintf(file, "%21.14E ", x);
		Ex = GEMxx(myConfig, x, x_, y, y_, z, z_, tol)*Mx+
		     GEMxy(myConfig, x, x_, y, y_, z, z_, tol)*My+
			 GEMxz(myConfig, x, x_, y, y_, z, z_, tol)*Mz;
		Ey = GEMyx(myConfig, x, x_, y, y_, z, z_, tol)*Mx+
		     GEMyy(myConfig, x, x_, y, y_, z, z_, tol)*My+
			 GEMyz(myConfig, x, x_, y, y_, z, z_, tol)*Mz;
		Ez = GEMzx(myConfig, x, x_, y, y_, z, z_, tol)*Mx+
		     GEMzy(myConfig, x, x_, y, y_, z, z_, tol)*My+
			 GEMzz(myConfig, x, x_, y, y_, z, z_, tol)*Mz;
		fprintf(file, "%21.14E %21.14E %21.14E ", cabs(Ex), cabs(Ey), cabs(Ez));
		Ex = GHMxx(myConfig, x, x_, y, y_, z, z_, tol)*Mx+
		     GHMxy(myConfig, x, x_, y, y_, z, z_, tol)*My+
			 GHMxz(myConfig, x, x_, y, y_, z, z_, tol)*Mz;
		Ey = GHMyx(myConfig, x, x_, y, y_, z, z_, tol)*Mx+
		     GHMyy(myConfig, x, x_, y, y_, z, z_, tol)*My+
			 GHMyz(myConfig, x, x_, y, y_, z, z_, tol)*Mz;
		Ez = GHMzx(myConfig, x, x_, y, y_, z, z_, tol)*Mx+
		     GHMzy(myConfig, x, x_, y, y_, z, z_, tol)*My+
			 GHMzz(myConfig, x, x_, y, y_, z, z_, tol)*Mz;
		fprintf(file, "%21.14E %21.14E %21.14E ", cabs(Ex), cabs(Ey), cabs(Ez));
		fprintf(file, "\n");
		progressBar(i, Ns);
	}
	fclose(file);
	tocTimer(&T);
	//
	unsetConfig(myConfig);
}

void TestFEKO2(){
	Config *myConfig=createFEKO();
	setConfig(myConfig);
	saveConfig(myConfig);
	//
	int Ns=1000;
	double x=+0.0;
	double y=+1.0;
	double z_min=-1.0;
	double z_max=+1.0;
	double x_=+0.0;
	double y_=+0.0;
	double z_=-0.4;
	double theta0=deg2rad(20.0);
	double phi0=deg2rad(30.0);
	double tol=1.0E-3;
	//
	double Jx=sin(theta0)*cos(phi0);
	double Jy=sin(theta0)*sin(phi0);
	double Jz=cos(theta0);
	double Mx=sin(theta0)*cos(phi0);
	double My=sin(theta0)*sin(phi0);
	double Mz=cos(theta0);
	double dz=(z_max-z_min)/(Ns-1.0);
	double z;
	complex double Ex, Ey, Ez;
	// Test J Dipole
	printf("Computing J Dipole:\n");
	Timer T;
	ticTimer(&T);
	FILE *file=fopen("Data/Validation/VCut/Data1.dat", "w");
	assert(file!=NULL);
	for (int i=0; i<Ns; i++){
		z = z_min+i*dz;
		fprintf(file, "%21.14E ", z);
		Ex = GEJxx(myConfig, x, x_, y, y_, z, z_, tol)*Jx+
		     GEJxy(myConfig, x, x_, y, y_, z, z_, tol)*Jy+
			 GEJxz(myConfig, x, x_, y, y_, z, z_, tol)*Jz;
		Ey = GEJyx(myConfig, x, x_, y, y_, z, z_, tol)*Jx+
		     GEJyy(myConfig, x, x_, y, y_, z, z_, tol)*Jy+
			 GEJyz(myConfig, x, x_, y, y_, z, z_, tol)*Jz;
		Ez = GEJzx(myConfig, x, x_, y, y_, z, z_, tol)*Jx+
		     GEJzy(myConfig, x, x_, y, y_, z, z_, tol)*Jy+
			 GEJzz(myConfig, x, x_, y, y_, z, z_, tol)*Jz;
		fprintf(file, "%21.14E %21.14E %21.14E ", cabs(Ex), cabs(Ey), cabs(Ez));
		Ex = GHJxx(myConfig, x, x_, y, y_, z, z_, tol)*Jx+
		     GHJxy(myConfig, x, x_, y, y_, z, z_, tol)*Jy+
			 GHJxz(myConfig, x, x_, y, y_, z, z_, tol)*Jz;
		Ey = GHJyx(myConfig, x, x_, y, y_, z, z_, tol)*Jx+
		     GHJyy(myConfig, x, x_, y, y_, z, z_, tol)*Jy+
			 GHJyz(myConfig, x, x_, y, y_, z, z_, tol)*Jz;
		Ez = GHJzx(myConfig, x, x_, y, y_, z, z_, tol)*Jx+
		     GHJzy(myConfig, x, x_, y, y_, z, z_, tol)*Jy+
			 GHJzz(myConfig, x, x_, y, y_, z, z_, tol)*Jz;
		fprintf(file, "%21.14E %21.14E %21.14E ", cabs(Ex), cabs(Ey), cabs(Ez));
		fprintf(file, "\n");
		progressBar(i, Ns);
	}
	fclose(file);
	tocTimer(&T);
	// Test M Dipole
	printf("Computing M Dipole:\n");
	ticTimer(&T);
	file = fopen("Data/Validation/VCut/Data2.dat", "w");
	assert(file!=NULL);
	for (int i=0; i<Ns; i++){
		z = z_min+i*dz;
		fprintf(file, "%21.14E ", z);
		Ex = GEMxx(myConfig, x, x_, y, y_, z, z_, tol)*Mx+
		     GEMxy(myConfig, x, x_, y, y_, z, z_, tol)*My+
			 GEMxz(myConfig, x, x_, y, y_, z, z_, tol)*Mz;
		Ey = GEMyx(myConfig, x, x_, y, y_, z, z_, tol)*Mx+
		     GEMyy(myConfig, x, x_, y, y_, z, z_, tol)*My+
			 GEMyz(myConfig, x, x_, y, y_, z, z_, tol)*Mz;
		Ez = GEMzx(myConfig, x, x_, y, y_, z, z_, tol)*Mx+
		     GEMzy(myConfig, x, x_, y, y_, z, z_, tol)*My+
			 GEMzz(myConfig, x, x_, y, y_, z, z_, tol)*Mz;
		fprintf(file, "%21.14E %21.14E %21.14E ", cabs(Ex), cabs(Ey), cabs(Ez));
		Ex = GHMxx(myConfig, x, x_, y, y_, z, z_, tol)*Mx+
		     GHMxy(myConfig, x, x_, y, y_, z, z_, tol)*My+
			 GHMxz(myConfig, x, x_, y, y_, z, z_, tol)*Mz;
		Ey = GHMyx(myConfig, x, x_, y, y_, z, z_, tol)*Mx+
		     GHMyy(myConfig, x, x_, y, y_, z, z_, tol)*My+
			 GHMyz(myConfig, x, x_, y, y_, z, z_, tol)*Mz;
		Ez = GHMzx(myConfig, x, x_, y, y_, z, z_, tol)*Mx+
		     GHMzy(myConfig, x, x_, y, y_, z, z_, tol)*My+
			 GHMzz(myConfig, x, x_, y, y_, z, z_, tol)*Mz;
		fprintf(file, "%21.14E %21.14E %21.14E ", cabs(Ex), cabs(Ey), cabs(Ez));
		fprintf(file, "\n");
		progressBar(i, Ns);
	}
	fclose(file);
	tocTimer(&T);
	//
	unsetConfig(myConfig);
}

void TestReflection(Config *myConfig){
	int N=myConfig->N;
	int Ns=1000;
	double theta_min=0.0;
	double theta_max=pi/2.0;
	double d_theta=(theta_max-theta_min)/(Ns-1.0);
	FILE *file=fopen("Data/Reflection/Data.dat", "w");
	complex double k_rho;
	complex double Gamma;
	double theta;
	complex double k1, kN;
	k1 = myConfig->Layers[0].k;
	kN = myConfig->Layers[N-1].k;
	for (int i=0; i<Ns; i++){
		theta = theta_min+i*d_theta;
		fprintf(file, "%21.14E ", theta*180.0/pi);
		k_rho = k1*sin(theta);
		Gamma = Refl(myConfig, 0, 'e', 'd', k_rho);
		fprintf(file, "%21.14E ", cabs(Gamma));
		Gamma = Refl(myConfig, 0, 'h', 'd', k_rho);
		fprintf(file, "%21.14E ", cabs(Gamma));
		k_rho = kN*sin(theta);
		Gamma = Refl(myConfig, N-1, 'e', 'u', k_rho);
		fprintf(file, "%21.14E ", cabs(Gamma));
		Gamma = Refl(myConfig, N-1, 'h', 'u', k_rho);
		fprintf(file, "%21.14E ", cabs(Gamma));
		fprintf(file, "\n");
	}
	fclose(file);
}

void TestTLGF(Config *myConfig, double z_){
	int N=myConfig->N;
	int Ns=1000;
	double z_min=myConfig->Layers[N-1].z_d;
	double z_max=myConfig->Layers[0].z_u;
	double dz=(z_max-z_min)/(Ns-1.0);
	FILE *file=fopen("Data/TLGF/Data.dat", "w");
	complex double k_rho= 0.7*myConfig->k0;
	complex double V, I;
	double z;
	for (int i=0; i<Ns; i++){
		z = z_min+i*dz;
		fprintf(file, "%21.14E ", z);
		TLGF(myConfig, 'e', 'v', z, z_, k_rho, &V, &I);
		fprintf(file, "%21.14E %21.14E ", creal(V), cimag(V));
		fprintf(file, "%21.14E %21.14E ", creal(I), cimag(I));
		TLGF(myConfig, 'h', 'v', z, z_, k_rho, &V, &I);
		fprintf(file, "%21.14E %21.14E ", creal(V), cimag(V));
		fprintf(file, "%21.14E %21.14E ", creal(I), cimag(I));
		TLGF(myConfig, 'e', 'i', z, z_, k_rho, &V, &I);
		fprintf(file, "%21.14E %21.14E ", creal(V), cimag(V));
		fprintf(file, "%21.14E %21.14E ", creal(I), cimag(I));
		TLGF(myConfig, 'h', 'i', z, z_, k_rho, &V, &I);
		fprintf(file, "%21.14E %21.14E ", creal(V), cimag(V));
		fprintf(file, "%21.14E %21.14E ", creal(I), cimag(I));
		fprintf(file, "\n");
	}
	fclose(file);
}

void TestTLGFr(Config *myConfig, double z_){
	int N=myConfig->N;
	int Ns=1000;
	double z_min=myConfig->Layers[N-1].z_d;
	double z_max=myConfig->Layers[0].z_u;
	double dz=(z_max-z_min)/(Ns-1.0);
	FILE *file=fopen("Data/TLGFr/Data.dat", "w");
	complex double k_rho= 0.7*myConfig->k0;
	complex double V, I;
	double z;
	for (int i=0; i<Ns; i++){
		z = z_min+i*dz;
		fprintf(file, "%21.14E ", z);
		TLGFr(myConfig, 'e', 'v', z, z_, k_rho, &V, &I);
		fprintf(file, "%21.14E %21.14E ", creal(V), cimag(V));
		fprintf(file, "%21.14E %21.14E ", creal(I), cimag(I));
		TLGFr(myConfig, 'h', 'v', z, z_, k_rho, &V, &I);
		fprintf(file, "%21.14E %21.14E ", creal(V), cimag(V));
		fprintf(file, "%21.14E %21.14E ", creal(I), cimag(I));
		TLGFr(myConfig, 'e', 'i', z, z_, k_rho, &V, &I);
		fprintf(file, "%21.14E %21.14E ", creal(V), cimag(V));
		fprintf(file, "%21.14E %21.14E ", creal(I), cimag(I));
		TLGFr(myConfig, 'h', 'i', z, z_, k_rho, &V, &I);
		fprintf(file, "%21.14E %21.14E ", creal(V), cimag(V));
		fprintf(file, "%21.14E %21.14E ", creal(I), cimag(I));
		fprintf(file, "\n");
	}
	fclose(file);
}