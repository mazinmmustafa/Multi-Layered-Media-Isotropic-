//
#include "../include/Engine.h"

void spectralParameters(Config *myConfig, int n, 
	char pol, complex double *Z, complex double *kz, 
	complex double *Theta, complex double k_rho){
	assert(n>=0&&n<myConfig->N);
	complex double mu=myConfig->Layers[n].mu;
	complex double eps=myConfig->Layers[n].eps;
	complex double k=myConfig->Layers[n].k;
	double d=myConfig->Layers[n].d;
	double omega=myConfig->omega;
	(*kz) = csqrt(k*k-k_rho*k_rho);
	if (cimag(*kz)>0){(*kz)*=-1.0;}
	if (pol=='e'){
		(*Z) = (*kz)/(omega*eps0*eps);
	}
	if (pol=='h'){
		(*Z) = (omega*mu0*mu)/(*kz);
	}
	(*Theta) = (*kz)*d;
}

complex double Refl(Config *myConfig, int n, char pol, 
	char dir, complex double k_rho){
	int N=myConfig->N;
	complex double j=csqrt(-1.0);
	assert(pol=='e'||pol=='h');
	assert(dir=='d'||dir=='u');
	complex double Gamma=0.0;
	complex double Z_m, Z_n, Theta, kz;
	complex Gamma_m_n, sigma_s, Omega_m_n;
	complex double A, B, C, D;
	if (dir=='d'){
		Gamma = myConfig->Gamma_d;
		for (int i=N-2; i>=n; i--){
			sigma_s = myConfig->Layers[i].sigma_s;
			spectralParameters(myConfig, i, pol, &Z_n, &kz, &Theta, k_rho);
			spectralParameters(myConfig, i+1, pol, &Z_m, &kz, &Theta, k_rho);
			Gamma_m_n = (Z_m-Z_n)/(Z_m+Z_n);
			Omega_m_n = Z_m*Z_n/(Z_m+Z_n);
			A = Gamma_m_n-Omega_m_n*sigma_s;
			B = (1.0-Omega_m_n*sigma_s)*Gamma*cexp(-j*2.0*Theta);
			C = 1.0+Omega_m_n*sigma_s;
			D = (Gamma_m_n+Omega_m_n*sigma_s)*Gamma*cexp(-j*2.0*Theta);
			Gamma = (A+B)/(C+D);
		}
	}
	if (dir=='u'){
		Gamma = myConfig->Gamma_u;
		for (int i=1; i<=n; i++){
			sigma_s = myConfig->Layers[i-1].sigma_s;
			spectralParameters(myConfig, i, pol, &Z_n, &kz, &Theta, k_rho);
			spectralParameters(myConfig, i-1, pol, &Z_m, &kz, &Theta, k_rho);
			Gamma_m_n = (Z_m-Z_n)/(Z_m+Z_n);
			Omega_m_n = Z_m*Z_n/(Z_m+Z_n);
			A = Gamma_m_n-Omega_m_n*sigma_s;
			B = (1.0-Omega_m_n*sigma_s)*Gamma*cexp(-j*2.0*Theta);
			C = 1.0+Omega_m_n*sigma_s;
			D = (Gamma_m_n+Omega_m_n*sigma_s)*Gamma*cexp(-j*2.0*Theta);
			Gamma = (A+B)/(C+D);
		}
	}
	return Gamma;
}

void TLGF(Config *myConfig, char pol, char sor, 
	double z, double z_, complex double k_rho, 
	complex double *V, complex double *I){
	assert(pol=='e'||pol=='h');
	assert(sor=='v'||sor=='i');
	complex double j=csqrt(-1.0);
	double v, i;
	if (sor=='v'){v = 1.0; i = 0.0;}
	if (sor=='i'){v = 0.0; i = 1.0;}
	int m=findLayer(myConfig, z_);
	int n=findLayer(myConfig, z);
	// Source layer
	complex double Gamma_d, Gamma_u, D_n, kz, Z, Vp, Vm; 
	complex double Theta, Gamma_d_, Gamma_u_;
	complex double Gamma_d_m, Gamma_d_n, Gamma_u_m, Gamma_u_n;
	complex double Theta_m, Theta_n;
	Gamma_d = Refl(myConfig, m, pol, 'd', k_rho);
	Gamma_u = Refl(myConfig, m, pol, 'u', k_rho);
	complex double A, B, C, D;
	double z_n, z_m;
	z_n=myConfig->Layers[m].z_d;
	z_m=myConfig->Layers[m].z_u;
	spectralParameters(myConfig, m, pol, &Z, &kz, 
		&Theta, k_rho);
	D_n = 1.0-Gamma_d*Gamma_u*cexp(-j*2.0*Theta);
	A = 1.0-Gamma_d*cexp(-j*2.0*kz*(z_-z_n));
	B = 1.0+Gamma_d*cexp(-j*2.0*kz*(z_-z_n));
	C = 1.0+Gamma_u*cexp(-j*2.0*kz*(z_m-z_));
	D = 1.0-Gamma_u*cexp(-j*2.0*kz*(z_m-z_));
	Vp = (A*v+B*i*Z)/(2.0*D_n);
	Vm = (C*i*Z-D*v)/(2.0*D_n);
	// Case m=n
	if (m==n){
		Gamma_u_ = Gamma_u*cexp(-j*2.0*kz*(z_m-z_));
		Gamma_d_ = Gamma_d*cexp(-j*2.0*kz*(z_-z_n));
		if (z>=z_){
			A = cexp(-j*kz*(z-z_));
			B = Gamma_u_*cexp(+j*kz*(z-z_));
			(*V) = +Vp*(A+B);
			(*I) = +Vp*(A-B)/Z;
		}
		if (z<z_){
			A = cexp(+j*kz*(z-z_));
			B = Gamma_d_*cexp(-j*kz*(z-z_));
			(*V) = +Vm*(A+B);
			(*I) = -Vm*(A-B)/Z;
		}
	}
	// Case m<n
	if (m<n){
		Gamma_d_m = Gamma_d;
		Gamma_d_n = Refl(myConfig, m+1, pol, 'd', k_rho);
		A = (1.0+Gamma_d_m)*cexp(-j*kz*(z_-z_n));
		spectralParameters(myConfig, m+1, pol, &Z, &kz, 
		&Theta, k_rho);
		B = 1.0+Gamma_d_n*cexp(-j*2.0*Theta);
		Vm*=(A/B);
		for (int i=m+2; i<=n; i++){
			Gamma_d_m = Refl(myConfig, i-1, pol, 'd', k_rho);
			Gamma_d_n = Refl(myConfig, i, pol, 'd', k_rho);
			spectralParameters(myConfig, i-1, pol, &Z, &kz, 
				&Theta_m, k_rho);
			spectralParameters(myConfig, i, pol, &Z, &kz, 
				&Theta_n, k_rho);
			A = (1.0+Gamma_d_m)*cexp(-j*Theta_m);
			B = 1.0+Gamma_d_n*cexp(-j*2.0*Theta_n);
			Vm*=(A/B);
		}
		Gamma_d_n = Refl(myConfig, n, pol, 'd', k_rho);
		spectralParameters(myConfig, n, pol, &Z, &kz, 
			&Theta_n, k_rho);
		z_n=myConfig->Layers[n].z_d;
		A = cexp(+j*kz*(z-z_n));
		B = Gamma_d_n*cexp(-j*kz*(z-z_n));
		(*V) = +Vm*cexp(-j*Theta_n)*(A+B);
		(*I) = -Vm*cexp(-j*Theta_n)*(A-B)/Z;
	}
	// Case m>n
	if (m>n){
		Gamma_u_m = Gamma_u;
		Gamma_u_n = Refl(myConfig, m-1, pol, 'u', k_rho);
		A = (1.0+Gamma_u_m)*cexp(-j*kz*(z_m-z_));
		spectralParameters(myConfig, m-1, pol, &Z, &kz, 
		&Theta, k_rho);
		B = 1.0+Gamma_u_n*cexp(-j*2.0*Theta);
		Vp*=(A/B);
		for (int i=m-2; i>=n; i--){
			Gamma_u_m = Refl(myConfig, i+1, pol, 'u', k_rho);
			Gamma_u_n = Refl(myConfig, i, pol, 'u', k_rho);
			spectralParameters(myConfig, i+1, pol, &Z, &kz, 
				&Theta_m, k_rho);
			spectralParameters(myConfig, i, pol, &Z, &kz, 
				&Theta_n, k_rho);
			A = (1.0+Gamma_u_m)*cexp(-j*Theta_m);
			B = 1.0+Gamma_u_n*cexp(-j*2.0*Theta_n);
			Vp*=(A/B);
		}
		Gamma_u_n = Refl(myConfig, n, pol, 'u', k_rho);
		spectralParameters(myConfig, n, pol, &Z, &kz, 
			&Theta_n, k_rho);
		z_m=myConfig->Layers[n].z_u;
		A = cexp(-j*kz*(z-z_m));
		B = Gamma_u_n*cexp(+j*kz*(z-z_m));
		(*V) = +Vp*cexp(-j*Theta_n)*(A+B);
		(*I) = +Vp*cexp(-j*Theta_n)*(A-B)/Z;
	}
}

void TLGFr(Config *myConfig, char pol, char sor, 
	double z, double z_, complex double k_rho, 
	complex double *V, complex double *I){
	assert(pol=='e'||pol=='h');
	assert(sor=='v'||sor=='i');
	complex double j=csqrt(-1.0);
	double v, i;
	if (sor=='v'){v = 1.0; i = 0.0;}
	if (sor=='i'){v = 0.0; i = 1.0;}
	int m=findLayer(myConfig, z_);
	int n=findLayer(myConfig, z);
	// Source layer
	complex double Gamma_d, Gamma_u, D_n, kz, Z, Vp, Vm; 
	complex double Theta, Gamma_d_, Gamma_u_;
	complex double Gamma_d_m, Gamma_d_n, Gamma_u_m, Gamma_u_n;
	complex double Theta_m, Theta_n;
	Gamma_d = Refl(myConfig, m, pol, 'd', k_rho);
	Gamma_u = Refl(myConfig, m, pol, 'u', k_rho);
	complex double A, B, C, D;
	double z_n, z_m;
	z_n=myConfig->Layers[m].z_d;
	z_m=myConfig->Layers[m].z_u;
	spectralParameters(myConfig, m, pol, &Z, &kz, 
		&Theta, k_rho);
	D_n = 1.0-Gamma_d*Gamma_u*cexp(-j*2.0*Theta);
	A = 1.0-Gamma_d*cexp(-j*2.0*kz*(z_-z_n));
	B = 1.0+Gamma_d*cexp(-j*2.0*kz*(z_-z_n));
	C = 1.0+Gamma_u*cexp(-j*2.0*kz*(z_m-z_));
	D = 1.0-Gamma_u*cexp(-j*2.0*kz*(z_m-z_));
	Vp = (A*v+B*i*Z)/(2.0*D_n);
	Vm = (C*i*Z-D*v)/(2.0*D_n);
	// Case m=n
	if (m==n){
		complex double Vp0, Vm0;
		Vp0 = (v+i*Z)/2.0;
		Vm0 = (i*Z-v)/2.0;
		Gamma_u_ = Gamma_u*cexp(-j*2.0*kz*(z_m-z_));
		Gamma_d_ = Gamma_d*cexp(-j*2.0*kz*(z_-z_n));
		if (z>=z_){
			A = cexp(-j*kz*(z-z_));
			B = Gamma_u_*cexp(+j*kz*(z-z_));
			(*V) = (Vp-Vp0)*A+Vp*B;
			(*I) = ((Vp-Vp0)*A-Vp*B)/Z;
		}
		if (z<z_){
			A = cexp(+j*kz*(z-z_));
			B = Gamma_d_*cexp(-j*kz*(z-z_));
			(*V) = (Vm-Vm0)*A+Vm*B;
			(*I) = ((Vm0-Vm)*A+Vm*B)/Z;
		}
	}
	// Case m<n
	if (m<n){
		Gamma_d_m = Gamma_d;
		Gamma_d_n = Refl(myConfig, m+1, pol, 'd', k_rho);
		A = (1.0+Gamma_d_m)*cexp(-j*kz*(z_-z_n));
		spectralParameters(myConfig, m+1, pol, &Z, &kz, 
		&Theta, k_rho);
		B = 1.0+Gamma_d_n*cexp(-j*2.0*Theta);
		Vm*=(A/B);
		for (int i=m+2; i<=n; i++){
			Gamma_d_m = Refl(myConfig, i-1, pol, 'd', k_rho);
			Gamma_d_n = Refl(myConfig, i, pol, 'd', k_rho);
			spectralParameters(myConfig, i-1, pol, &Z, &kz, 
				&Theta_m, k_rho);
			spectralParameters(myConfig, i, pol, &Z, &kz, 
				&Theta_n, k_rho);
			A = (1.0+Gamma_d_m)*cexp(-j*Theta_m);
			B = 1.0+Gamma_d_n*cexp(-j*2.0*Theta_n);
			Vm*=(A/B);
		}
		Gamma_d_n = Refl(myConfig, n, pol, 'd', k_rho);
		spectralParameters(myConfig, n, pol, &Z, &kz, 
			&Theta_n, k_rho);
		z_n=myConfig->Layers[n].z_d;
		A = cexp(+j*kz*(z-z_n));
		B = Gamma_d_n*cexp(-j*kz*(z-z_n));
		(*V) = +Vm*cexp(-j*Theta_n)*(A+B);
		(*I) = -Vm*cexp(-j*Theta_n)*(A-B)/Z;
	}
	// Case m>n
	if (m>n){
		Gamma_u_m = Gamma_u;
		Gamma_u_n = Refl(myConfig, m-1, pol, 'u', k_rho);
		A = (1.0+Gamma_u_m)*cexp(-j*kz*(z_m-z_));
		spectralParameters(myConfig, m-1, pol, &Z, &kz, 
		&Theta, k_rho);
		B = 1.0+Gamma_u_n*cexp(-j*2.0*Theta);
		Vp*=(A/B);
		for (int i=m-2; i>=n; i--){
			Gamma_u_m = Refl(myConfig, i+1, pol, 'u', k_rho);
			Gamma_u_n = Refl(myConfig, i, pol, 'u', k_rho);
			spectralParameters(myConfig, i+1, pol, &Z, &kz, 
				&Theta_m, k_rho);
			spectralParameters(myConfig, i, pol, &Z, &kz, 
				&Theta_n, k_rho);
			A = (1.0+Gamma_u_m)*cexp(-j*Theta_m);
			B = 1.0+Gamma_u_n*cexp(-j*2.0*Theta_n);
			Vp*=(A/B);
		}
		Gamma_u_n = Refl(myConfig, n, pol, 'u', k_rho);
		spectralParameters(myConfig, n, pol, &Z, &kz, 
			&Theta_n, k_rho);
		z_m=myConfig->Layers[n].z_u;
		A = cexp(-j*kz*(z-z_m));
		B = Gamma_u_n*cexp(+j*kz*(z-z_m));
		(*V) = +Vp*cexp(-j*Theta_n)*(A+B);
		(*I) = +Vp*cexp(-j*Theta_n)*(A-B)/Z;
	}
}

int findLayer(Config *myConfig, double z){
	int n=0;
	for (int i=0; i<myConfig->N; i++){
		if ((z>=myConfig->Layers[i].z_d)&&
		    (z<=myConfig->Layers[i].z_u)){
			n = i;
			break;
		}
	}
	return n;
}

void findMax(Config *myConfig, double *a, double *b, double rho,
	double z, double z_){
	double max=cabs(csqrt(myConfig->Layers[0].mu*myConfig->Layers[0].eps));
	double temp;
	for (int i=1; i<myConfig->N; i++){
		temp = cabs(csqrt(myConfig->Layers[i].mu*myConfig->Layers[i].eps));
		if (temp>max){
			max = temp;
		}
	}
	(*a) = myConfig->k0*(1.0+max);
	if (rho>fabs(z-z_)){
		(*b) = myConfig->k0 < 1.0/rho ? myConfig->k0 : 1.0/rho;
	}else{
		(*b) = myConfig->k0;
	}
}

void spectralDetour(complex double t, complex double *k_rho, 
	complex double *dk_rho, double a, double b){
	complex double j=csqrt(-1.0);
	if (cabs(t)>=0.0&&cabs(t)<=a){
		(*k_rho) = t+j*b*csin(pi*t/a);
		(*dk_rho) = (1.0+j*pi*(b/a)*ccos(pi*t/a));
	}
	else{
		(*k_rho) = t+j*0.0;
		(*dk_rho) = 1.0;
	}
}

complex double GEJ0xx(double x, double x_, double y, double y_,
	double z, double z_, double omega, complex double mu, 
	complex double eps, double k0){
	complex double j=csqrt(-1.0);
	complex double k=k0*csqrt(mu*eps);
	double rho=sqrt((x-x_)*(x-x_)+(y-y_)*(y-y_));
	double r=sqrt(rho*rho+(z-z_)*(z-z_));
	complex double g=cexp(-j*k*r)/(4.0*pi*r);
	complex double dg1=-((1.0+j*k*r)/r)*g;
	complex double dg2=-((1.0+j*k*r)/r)*dg1+g/(r*r);
	if (r==0.0){
		return -j*omega*mu0*mu*(-1.0/(3.0*k*k));
	}else{
		return -j*omega*mu0*mu*(g+(1.0/(k*k))*(dg1/r+((x-x_)*(x-x_)/(r*r))*(dg2-dg1/r)));
	}
}

complex double GEJ0xy(double x, double x_, double y, double y_,
	double z, double z_, double omega, complex double mu, 
	complex double eps, double k0){
	complex double j=csqrt(-1.0);
	complex double k=k0*csqrt(mu*eps);
	double rho=sqrt((x-x_)*(x-x_)+(y-y_)*(y-y_));
	double r=sqrt(rho*rho+(z-z_)*(z-z_));
	complex double g=cexp(-j*k*r)/(4.0*pi*r);
	complex double dg1=-((1.0+j*k*r)/r)*g;
	complex double dg2=-((1.0+j*k*r)/r)*dg1+g/(r*r);
	if (r==0.0){
		return 0.0;
	}else{
		return -j*omega*mu0*mu*((1.0/(k*k))*(((x-x_)*(y-y_)/(r*r))*(dg2-dg1/r)));
	}
}

complex double GEJ0xz(double x, double x_, double y, double y_,
	double z, double z_, double omega, complex double mu, 
	complex double eps, double k0){
	complex double j=csqrt(-1.0);
	complex double k=k0*csqrt(mu*eps);
	double rho=sqrt((x-x_)*(x-x_)+(y-y_)*(y-y_));
	double r=sqrt(rho*rho+(z-z_)*(z-z_));
	complex double g=cexp(-j*k*r)/(4.0*pi*r);
	complex double dg1=-((1.0+j*k*r)/r)*g;
	complex double dg2=-((1.0+j*k*r)/r)*dg1+g/(r*r);
	if (r==0.0){
		return 0.0;
	}else{
		return -j*omega*mu0*mu*((1.0/(k*k))*(((x-x_)*(z-z_)/(r*r))*(dg2-dg1/r)));
	}
}

complex double GEJ0yx(double x, double x_, double y, double y_,
	double z, double z_, double omega, complex double mu, 
	complex double eps, double k0){
	complex double j=csqrt(-1.0);
	complex double k=k0*csqrt(mu*eps);
	double rho=sqrt((x-x_)*(x-x_)+(y-y_)*(y-y_));
	double r=sqrt(rho*rho+(z-z_)*(z-z_));
	complex double g=cexp(-j*k*r)/(4.0*pi*r);
	complex double dg1=-((1.0+j*k*r)/r)*g;
	complex double dg2=-((1.0+j*k*r)/r)*dg1+g/(r*r);
	if (r==0.0){
		return 0.0;
	}else{
		return -j*omega*mu0*mu*((1.0/(k*k))*(((y-y_)*(x-x_)/(r*r))*(dg2-dg1/r)));
	}
}

complex double GEJ0yy(double x, double x_, double y, double y_,
	double z, double z_, double omega, complex double mu, 
	complex double eps, double k0){
	complex double j=csqrt(-1.0);
	complex double k=k0*csqrt(mu*eps);
	double rho=sqrt((x-x_)*(x-x_)+(y-y_)*(y-y_));
	double r=sqrt(rho*rho+(z-z_)*(z-z_));
	complex double g=cexp(-j*k*r)/(4.0*pi*r);
	complex double dg1=-((1.0+j*k*r)/r)*g;
	complex double dg2=-((1.0+j*k*r)/r)*dg1+g/(r*r);
	if (r==0.0){
		return -j*omega*mu0*mu*(-1.0/(3.0*k*k));
	}else{
		return -j*omega*mu0*mu*(g+(1.0/(k*k))*(dg1/r+((y-y_)*(y-y_)/(r*r))*(dg2-dg1/r)));
	}
}

complex double GEJ0yz(double x, double x_, double y, double y_,
	double z, double z_, double omega, complex double mu, 
	complex double eps, double k0){
	complex double j=csqrt(-1.0);
	complex double k=k0*csqrt(mu*eps);
	double rho=sqrt((x-x_)*(x-x_)+(y-y_)*(y-y_));
	double r=sqrt(rho*rho+(z-z_)*(z-z_));
	complex double g=cexp(-j*k*r)/(4.0*pi*r);
	complex double dg1=-((1.0+j*k*r)/r)*g;
	complex double dg2=-((1.0+j*k*r)/r)*dg1+g/(r*r);
	if (r==0.0){
		return 0.0;
	}else{
		return -j*omega*mu0*mu*((1.0/(k*k))*(((y-y_)*(z-z_)/(r*r))*(dg2-dg1/r)));
	}
}

complex double GEJ0zx(double x, double x_, double y, double y_,
	double z, double z_, double omega, complex double mu, 
	complex double eps, double k0){
	complex double j=csqrt(-1.0);
	complex double k=k0*csqrt(mu*eps);
	double rho=sqrt((x-x_)*(x-x_)+(y-y_)*(y-y_));
	double r=sqrt(rho*rho+(z-z_)*(z-z_));
	complex double g=cexp(-j*k*r)/(4.0*pi*r);
	complex double dg1=-((1.0+j*k*r)/r)*g;
	complex double dg2=-((1.0+j*k*r)/r)*dg1+g/(r*r);
	if (r==0.0){
		return 0.0;
	}else{
		return -j*omega*mu0*mu*((1.0/(k*k))*(((z-z_)*(x-x_)/(r*r))*(dg2-dg1/r)));
	}
}

complex double GEJ0zy(double x, double x_, double y, double y_,
	double z, double z_, double omega, complex double mu, 
	complex double eps, double k0){
	complex double j=csqrt(-1.0);
	complex double k=k0*csqrt(mu*eps);
	double rho=sqrt((x-x_)*(x-x_)+(y-y_)*(y-y_));
	double r=sqrt(rho*rho+(z-z_)*(z-z_));
	complex double g=cexp(-j*k*r)/(4.0*pi*r);
	complex double dg1=-((1.0+j*k*r)/r)*g;
	complex double dg2=-((1.0+j*k*r)/r)*dg1+g/(r*r);
	if (r==0.0){
		return 0.0;
	}else{
		return -j*omega*mu0*mu*((1.0/(k*k))*(((z-z_)*(y-y_)/(r*r))*(dg2-dg1/r)));
	}
}

complex double GEJ0zz(double x, double x_, double y, double y_,
	double z, double z_, double omega, complex double mu, 
	complex double eps, double k0){
	complex double j=csqrt(-1.0);
	complex double k=k0*csqrt(mu*eps);
	double rho=sqrt((x-x_)*(x-x_)+(y-y_)*(y-y_));
	double r=sqrt(rho*rho+(z-z_)*(z-z_));
	complex double g=cexp(-j*k*r)/(4.0*pi*r);
	complex double dg1=-((1.0+j*k*r)/r)*g;
	complex double dg2=-((1.0+j*k*r)/r)*dg1+g/(r*r);
	if (r==0.0){
		return -j*omega*mu0*mu*(-1.0/(3.0*k*k));
	}else{
		return -j*omega*mu0*mu*(g+(1.0/(k*k))*(dg1/r+((z-z_)*(z-z_)/(r*r))*(dg2-dg1/r)));
	}
}

complex double GEJ1(complex double k_rho_in, void *args){
	GfuncArgs *myArgs=(GfuncArgs*) args;
	Config *myConfig=myArgs->myConfig;
	double z=myArgs->z;
	double z_=myArgs->z_;
	double rho=myArgs->rho;
	double a=myArgs->a;
	double b=myArgs->b;
	complex double k_rho, dk_rho;
	spectralDetour(k_rho_in, &k_rho, &dk_rho, a, b);
	complex double V_e_i, I_e_i, V_h_i, I_h_i;
	TLGFr(myConfig, 'e', 'i', z, z_, k_rho, &V_e_i, &I_e_i);
	TLGFr(myConfig, 'h', 'i', z, z_, k_rho, &V_h_i, &I_h_i);
	return (V_e_i+V_h_i)*besselj(0, k_rho*rho)*k_rho*dk_rho/(2.0*pi);
}

complex double GEJ2(complex double k_rho_in, void *args){
	GfuncArgs *myArgs=(GfuncArgs*) args;
	Config *myConfig=myArgs->myConfig;
	double z=myArgs->z;
	double z_=myArgs->z_;
	double rho=myArgs->rho;
	double a=myArgs->a;
	double b=myArgs->b;
	complex double k_rho, dk_rho;
	spectralDetour(k_rho_in, &k_rho, &dk_rho, a, b);
	complex double V_e_i, I_e_i, V_h_i, I_h_i;
	TLGFr(myConfig, 'e', 'i', z, z_, k_rho, &V_e_i, &I_e_i);
	TLGFr(myConfig, 'h', 'i', z, z_, k_rho, &V_h_i, &I_h_i);
	return (V_e_i-V_h_i)*besselj(2, k_rho*rho)*k_rho*dk_rho/(2.0*pi);	
}

complex double GEJ3(complex double k_rho_in, void *args){
	GfuncArgs *myArgs=(GfuncArgs*) args;
	Config *myConfig=myArgs->myConfig;
	double z=myArgs->z;
	double z_=myArgs->z_;
	double rho=myArgs->rho;
	double a=myArgs->a;
	double b=myArgs->b;
	complex double k_rho, dk_rho;
	spectralDetour(k_rho_in, &k_rho, &dk_rho, a, b);
	complex double V_e_v, I_e_v;
	TLGFr(myConfig, 'e', 'v', z, z_, k_rho, &V_e_v, &I_e_v);
	return k_rho*V_e_v*besselj(1, k_rho*rho)*k_rho*dk_rho/(2.0*pi);	
}

complex double GEJ4(complex double k_rho_in, void *args){
	GfuncArgs *myArgs=(GfuncArgs*) args;
	Config *myConfig=myArgs->myConfig;
	double z=myArgs->z;
	double z_=myArgs->z_;
	double rho=myArgs->rho;
	double a=myArgs->a;
	double b=myArgs->b;
	complex double k_rho, dk_rho;
	spectralDetour(k_rho_in, &k_rho, &dk_rho, a, b);
	complex double V_e_i, I_e_i;
	TLGFr(myConfig, 'e', 'i', z, z_, k_rho, &V_e_i, &I_e_i);
	return k_rho*I_e_i*besselj(1, k_rho*rho)*k_rho*dk_rho/(2.0*pi);	
}

complex double GEJ5(complex double k_rho_in, void *args){
	GfuncArgs *myArgs=(GfuncArgs*) args;
	Config *myConfig=myArgs->myConfig;
	double z=myArgs->z;
	double z_=myArgs->z_;
	double rho=myArgs->rho;
	double a=myArgs->a;
	double b=myArgs->b;
	complex double k_rho, dk_rho;
	spectralDetour(k_rho_in, &k_rho, &dk_rho, a, b);
	complex double V_e_v, I_e_v;
	TLGFr(myConfig, 'e', 'v', z, z_, k_rho, &V_e_v, &I_e_v);
	return k_rho*k_rho*I_e_v*besselj(0, k_rho*rho)*k_rho*dk_rho/(2.0*pi);	
}

complex double GEJxx(Config *myConfig, double x, double x_, 
	double y,  double y_, double z, double z_, double tol){
	double rho=sqrt((x-x_)*(x-x_)+(y-y_)*(y-y_));
	double phi=atan2(y-y_, x-x_);
	int n=findLayer(myConfig, z);
	int m=findLayer(myConfig, z_);
	double a, b;
	findMax(myConfig, &a, &b, rho, z, z_);
	complex double G1, G2;
	GfuncArgs myArgs={myConfig, z, z_, rho, a, b};
	G1 = QuadL_1D(GEJ1, &myArgs, 0.0, 1.5*a, tol);
	G2 = QuadL_1D(GEJ2, &myArgs, 0.0, 1.5*a, tol);
	complex double I=-0.5*G1+0.5*cos(2.0*phi)*G2;
	if (m==n){
		return I+GEJ0xx(x, x_, y, y_, z, z_,
			myConfig->omega, myConfig->Layers[n].mu, 
			myConfig->Layers[n].eps, myConfig->k0);
	}else{
		return I;
	}
}

complex double GEJxy(Config *myConfig, double x, double x_, 
	double y,  double y_, double z, double z_, double tol){
	double rho=sqrt((x-x_)*(x-x_)+(y-y_)*(y-y_));
	double phi=atan2(y-y_, x-x_);
	int n=findLayer(myConfig, z);
	int m=findLayer(myConfig, z_);
	double a, b;
	findMax(myConfig, &a, &b, rho, z, z_);
	complex double G1;
	GfuncArgs myArgs={myConfig, z, z_, rho, a, b};
	G1 = QuadL_1D(GEJ2, &myArgs, 0.0, 1.5*a, tol);
	complex double I=0.5*sin(2.0*phi)*G1;
	if (m==n){
		return I+GEJ0xy(x, x_, y, y_, z, z_,
			myConfig->omega, myConfig->Layers[n].mu, 
			myConfig->Layers[n].eps, myConfig->k0);
	}else{
		return I;
	}
}

complex double GEJyx(Config *myConfig, double x, double x_, 
	double y,  double y_, double z, double z_, double tol){
	double rho=sqrt((x-x_)*(x-x_)+(y-y_)*(y-y_));
	double phi=atan2(y-y_, x-x_);
	int n=findLayer(myConfig, z);
	int m=findLayer(myConfig, z_);
	double a, b;
	findMax(myConfig, &a, &b, rho, z, z_);
	complex double G1;
	GfuncArgs myArgs={myConfig, z, z_, rho, a, b};
	G1 = QuadL_1D(GEJ2, &myArgs, 0.0, 1.5*a, tol);
	complex double I=0.5*sin(2.0*phi)*G1;
	if (m==n){
		return I+GEJ0yx(x, x_, y, y_, z, z_,
			myConfig->omega, myConfig->Layers[n].mu, 
			myConfig->Layers[n].eps, myConfig->k0);
	}else{
		return I;
	}
}

complex double GEJyy(Config *myConfig, double x, double x_, 
	double y,  double y_, double z, double z_, double tol){
	double rho=sqrt((x-x_)*(x-x_)+(y-y_)*(y-y_));
	double phi=atan2(y-y_, x-x_);
	int n=findLayer(myConfig, z);
	int m=findLayer(myConfig, z_);
	double a, b;
	findMax(myConfig, &a, &b, rho, z, z_);
	complex double G1, G2;
	GfuncArgs myArgs={myConfig, z, z_, rho, a, b};
	G1 = QuadL_1D(GEJ1, &myArgs, 0.0, 1.5*a, tol);
	G2 = QuadL_1D(GEJ2, &myArgs, 0.0, 1.5*a, tol);
	complex double I=-0.5*G1-0.5*cos(2.0*phi)*G2;
	if (m==n){
		return I+GEJ0yy(x, x_, y, y_, z, z_,
			myConfig->omega, myConfig->Layers[n].mu, 
			myConfig->Layers[n].eps, myConfig->k0);
	}else{
		return I;
	}
}

complex double GEJxz(Config *myConfig, double x, double x_, 
	double y,  double y_, double z, double z_, double tol){
	double rho=sqrt((x-x_)*(x-x_)+(y-y_)*(y-y_));
	double phi=atan2(y-y_, x-x_);
	int n=findLayer(myConfig, z);
	int m=findLayer(myConfig, z_);
	double a, b;
	findMax(myConfig, &a, &b, rho, z, z_);
	complex double G1;
	GfuncArgs myArgs={myConfig, z, z_, rho, a, b};
	G1 = QuadL_1D(GEJ3, &myArgs, 0.0, 1.5*a, tol);
	complex double j=csqrt(-1.0);
	complex double eps_=myConfig->Layers[m].eps;
	double k0=myConfig->k0;
	complex double I=(eta0/(j*k0*eps_))*cos(phi)*G1;
	if (m==n){
		return I+GEJ0xz(x, x_, y, y_, z, z_,
			myConfig->omega, myConfig->Layers[n].mu, 
			myConfig->Layers[n].eps, myConfig->k0);
	}else{
		return I;
	}
}

complex double GEJyz(Config *myConfig, double x, double x_, 
	double y,  double y_, double z, double z_, double tol){
	double rho=sqrt((x-x_)*(x-x_)+(y-y_)*(y-y_));
	double phi=atan2(y-y_, x-x_);
	int n=findLayer(myConfig, z);
	int m=findLayer(myConfig, z_);
	double a, b;
	findMax(myConfig, &a, &b, rho, z, z_);
	complex double G1;
	GfuncArgs myArgs={myConfig, z, z_, rho, a, b};
	G1 = QuadL_1D(GEJ3, &myArgs, 0.0, 1.5*a, tol);
	complex double j=csqrt(-1.0);
	complex double eps_=myConfig->Layers[m].eps;
	double k0=myConfig->k0;
	complex double I=(eta0/(j*k0*eps_))*sin(phi)*G1;
	if (m==n){
		return I+GEJ0yz(x, x_, y, y_, z, z_,
			myConfig->omega, myConfig->Layers[n].mu, 
			myConfig->Layers[n].eps, myConfig->k0);
	}else{
		return I;
	}
}

complex double GEJzx(Config *myConfig, double x, double x_, 
	double y,  double y_, double z, double z_, double tol){
	double rho=sqrt((x-x_)*(x-x_)+(y-y_)*(y-y_));
	double phi=atan2(y-y_, x-x_);
	int n=findLayer(myConfig, z);
	int m=findLayer(myConfig, z_);
	double a, b;
	findMax(myConfig, &a, &b, rho, z, z_);
	complex double G1;
	GfuncArgs myArgs={myConfig, z, z_, rho, a, b};
	G1 = QuadL_1D(GEJ4, &myArgs, 0.0, 1.5*a, tol);
	complex double j=csqrt(-1.0);
	complex double eps=myConfig->Layers[n].eps;
	double k0=myConfig->k0;
	complex double I=(eta0/(j*k0*eps))*cos(phi)*G1;
	if (m==n){
		return I+GEJ0zx(x, x_, y, y_, z, z_,
			myConfig->omega, myConfig->Layers[n].mu, 
			myConfig->Layers[n].eps, myConfig->k0);
	}else{
		return I;
	}
}

complex double GEJzy(Config *myConfig, double x, double x_, 
	double y,  double y_, double z, double z_, double tol){
	double rho=sqrt((x-x_)*(x-x_)+(y-y_)*(y-y_));
	double phi=atan2(y-y_, x-x_);
	int n=findLayer(myConfig, z);
	int m=findLayer(myConfig, z_);
	double a, b;
	findMax(myConfig, &a, &b, rho, z, z_);
	complex double G1;
	GfuncArgs myArgs={myConfig, z, z_, rho, a, b};
	G1 = QuadL_1D(GEJ4, &myArgs, 0.0, 1.5*a, tol);
	complex double j=csqrt(-1.0);
	complex double eps=myConfig->Layers[n].eps;
	double k0=myConfig->k0;
	complex double I=(eta0/(j*k0*eps))*sin(phi)*G1;
	if (m==n){
		return I+GEJ0zy(x, x_, y, y_, z, z_,
			myConfig->omega, myConfig->Layers[n].mu, 
			myConfig->Layers[n].eps, myConfig->k0);
	}else{
		return I;
	}
}

complex double GEJzz(Config *myConfig, double x, double x_, 
	double y,  double y_, double z, double z_, double tol){
	double rho=sqrt((x-x_)*(x-x_)+(y-y_)*(y-y_));
	int n=findLayer(myConfig, z);
	int m=findLayer(myConfig, z_);
	double a, b;
	findMax(myConfig, &a, &b, rho, z, z_);
	complex double G1;
	GfuncArgs myArgs={myConfig, z, z_, rho, a, b};
	G1 = QuadL_1D(GEJ5, &myArgs, 0.0, 1.5*a, tol);
	complex double j=csqrt(-1.0);
	complex double eps_=myConfig->Layers[m].eps;
	complex double eps=myConfig->Layers[n].eps;
	double k0=myConfig->k0;
	complex double I=-(eta0*eta0/(k0*k0*eps_*eps))*G1;
	if (rho==0.0&&fabs(z-z_)==0.0){
		I-=(eta0/(j*k0*eps));
	}
	if (m==n){
		return I+GEJ0zz(x, x_, y, y_, z, z_,
			myConfig->omega, myConfig->Layers[n].mu, 
			myConfig->Layers[n].eps, myConfig->k0);
	}else{
		return I;
	}
}

complex double GEM0xy(double x, double x_, double y, double y_,
	double z, double z_, complex double mu,complex double eps, double k0){
	complex double j=csqrt(-1.0);
	complex double k=k0*csqrt(mu*eps);
	double rho=sqrt((x-x_)*(x-x_)+(y-y_)*(y-y_));
	double r=sqrt(rho*rho+(z-z_)*(z-z_));
	complex double g=cexp(-j*k*r)/(4.0*pi*r);
	complex double dg1=-((1.0+j*k*r)/r)*g;
	if (r==0.0){
		return 0.0;
	}else{
		return +(z-z_)*dg1/r;
	}
}

complex double GEM0xz(double x, double x_, double y, double y_,
	double z, double z_, complex double mu,complex double eps, double k0){
	complex double j=csqrt(-1.0);
	complex double k=k0*csqrt(mu*eps);
	double rho=sqrt((x-x_)*(x-x_)+(y-y_)*(y-y_));
	double r=sqrt(rho*rho+(z-z_)*(z-z_));
	complex double g=cexp(-j*k*r)/(4.0*pi*r);
	complex double dg1=-((1.0+j*k*r)/r)*g;
	if (r==0.0){
		return 0.0;
	}else{
		return -(y-y_)*dg1/r;
	}
}

complex double GEM0yx(double x, double x_, double y, double y_,
	double z, double z_, complex double mu,complex double eps, double k0){
	complex double j=csqrt(-1.0);
	complex double k=k0*csqrt(mu*eps);
	double rho=sqrt((x-x_)*(x-x_)+(y-y_)*(y-y_));
	double r=sqrt(rho*rho+(z-z_)*(z-z_));
	complex double g=cexp(-j*k*r)/(4.0*pi*r);
	complex double dg1=-((1.0+j*k*r)/r)*g;
	if (r==0.0){
		return 0.0;
	}else{
		return -(z-z_)*dg1/r;
	}
}

complex double GEM0yz(double x, double x_, double y, double y_,
	double z, double z_, complex double mu,complex double eps, double k0){
	complex double j=csqrt(-1.0);
	complex double k=k0*csqrt(mu*eps);
	double rho=sqrt((x-x_)*(x-x_)+(y-y_)*(y-y_));
	double r=sqrt(rho*rho+(z-z_)*(z-z_));
	complex double g=cexp(-j*k*r)/(4.0*pi*r);
	complex double dg1=-((1.0+j*k*r)/r)*g;
	if (r==0.0){
		return 0.0;
	}else{
		return +(x-x_)*dg1/r;
	}
}

complex double GEM0zx(double x, double x_, double y, double y_,
	double z, double z_, complex double mu,complex double eps, double k0){
	complex double j=csqrt(-1.0);
	complex double k=k0*csqrt(mu*eps);
	double rho=sqrt((x-x_)*(x-x_)+(y-y_)*(y-y_));
	double r=sqrt(rho*rho+(z-z_)*(z-z_));
	complex double g=cexp(-j*k*r)/(4.0*pi*r);
	complex double dg1=-((1.0+j*k*r)/r)*g;
	if (r==0.0){
		return 0.0;
	}else{
		return +(y-y_)*dg1/r;
	}
}

complex double GEM0zy(double x, double x_, double y, double y_,
	double z, double z_, complex double mu,complex double eps, double k0){
	complex double j=csqrt(-1.0);
	complex double k=k0*csqrt(mu*eps);
	double rho=sqrt((x-x_)*(x-x_)+(y-y_)*(y-y_));
	double r=sqrt(rho*rho+(z-z_)*(z-z_));
	complex double g=cexp(-j*k*r)/(4.0*pi*r);
	complex double dg1=-((1.0+j*k*r)/r)*g;
	if (r==0.0){
		return 0.0;
	}else{
		return -(x-x_)*dg1/r;
	}
}

complex double GEM1(complex double k_rho_in, void *args){
	GfuncArgs *myArgs=(GfuncArgs*) args;
	Config *myConfig=myArgs->myConfig;
	double z=myArgs->z;
	double z_=myArgs->z_;
	double rho=myArgs->rho;
	double a=myArgs->a;
	double b=myArgs->b;
	complex double k_rho, dk_rho;
	spectralDetour(k_rho_in, &k_rho, &dk_rho, a, b);
	complex double V_e_v, I_e_v, V_h_v, I_h_v;
	TLGFr(myConfig, 'e', 'v', z, z_, k_rho, &V_e_v, &I_e_v);
	TLGFr(myConfig, 'h', 'v', z, z_, k_rho, &V_h_v, &I_h_v);
	return (V_e_v-V_h_v)*besselj(2, k_rho*rho)*k_rho*dk_rho/(2.0*pi);
}

complex double GEM2(complex double k_rho_in, void *args){
	GfuncArgs *myArgs=(GfuncArgs*) args;
	Config *myConfig=myArgs->myConfig;
	double z=myArgs->z;
	double z_=myArgs->z_;
	double rho=myArgs->rho;
	double a=myArgs->a;
	double b=myArgs->b;
	complex double k_rho, dk_rho;
	spectralDetour(k_rho_in, &k_rho, &dk_rho, a, b);
	complex double V_e_v, I_e_v, V_h_v, I_h_v;
	TLGFr(myConfig, 'e', 'v', z, z_, k_rho, &V_e_v, &I_e_v);
	TLGFr(myConfig, 'h', 'v', z, z_, k_rho, &V_h_v, &I_h_v);
	return (V_e_v+V_h_v)*besselj(0, k_rho*rho)*k_rho*dk_rho/(2.0*pi);	
}

complex double GEM3(complex double k_rho_in, void *args){
	GfuncArgs *myArgs=(GfuncArgs*) args;
	Config *myConfig=myArgs->myConfig;
	double z=myArgs->z;
	double z_=myArgs->z_;
	double rho=myArgs->rho;
	double a=myArgs->a;
	double b=myArgs->b;
	complex double k_rho, dk_rho;
	spectralDetour(k_rho_in, &k_rho, &dk_rho, a, b);
	complex double V_h_i, I_h_i;
	TLGFr(myConfig, 'h', 'i', z, z_, k_rho, &V_h_i, &I_h_i);
	return k_rho*V_h_i*besselj(1, k_rho*rho)*k_rho*dk_rho/(2.0*pi);	
}

complex double GEM4(complex double k_rho_in, void *args){
	GfuncArgs *myArgs=(GfuncArgs*) args;
	Config *myConfig=myArgs->myConfig;
	double z=myArgs->z;
	double z_=myArgs->z_;
	double rho=myArgs->rho;
	double a=myArgs->a;
	double b=myArgs->b;
	complex double k_rho, dk_rho;
	spectralDetour(k_rho_in, &k_rho, &dk_rho, a, b);
	complex double V_e_v, I_e_v;
	TLGFr(myConfig, 'e', 'v', z, z_, k_rho, &V_e_v, &I_e_v);
	return k_rho*I_e_v*besselj(1, k_rho*rho)*k_rho*dk_rho/(2.0*pi);	
}

complex double GEMxx(Config *myConfig, double x, double x_, 
	double y,  double y_, double z, double z_, double tol){
	double rho=sqrt((x-x_)*(x-x_)+(y-y_)*(y-y_));
	double phi=atan2(y-y_, x-x_);
	double a, b;
	findMax(myConfig, &a, &b, rho, z, z_);
	complex double G1;
	GfuncArgs myArgs={myConfig, z, z_, rho, a, b};
	G1 = QuadL_1D(GEM1, &myArgs, 0.0, 1.5*a, tol);
	complex double I=-0.5*sin(2.0*phi)*G1;
	return I;
}

complex double GEMxy(Config *myConfig, double x, double x_, 
	double y,  double y_, double z, double z_, double tol){
	double rho=sqrt((x-x_)*(x-x_)+(y-y_)*(y-y_));
	double phi=atan2(y-y_, x-x_);
	int n=findLayer(myConfig, z);
	int m=findLayer(myConfig, z_);
	double a, b;
	findMax(myConfig, &a, &b, rho, z, z_);
	complex double G1, G2;
	GfuncArgs myArgs={myConfig, z, z_, rho, a, b};
	G1 = QuadL_1D(GEM1, &myArgs, 0.0, 1.5*a, tol);
	G2 = QuadL_1D(GEM2, &myArgs, 0.0, 1.5*a, tol);
	complex double I=-0.5*G2+0.5*cos(2.0*phi)*G1;
	if (m==n){
		return I+GEM0xy(x, x_, y, y_, z, z_, myConfig->Layers[n].mu, 
			myConfig->Layers[n].eps, myConfig->k0);
	}else{
		return I;
	}
}

complex double GEMyx(Config *myConfig, double x, double x_, 
	double y,  double y_, double z, double z_, double tol){
	double rho=sqrt((x-x_)*(x-x_)+(y-y_)*(y-y_));
	double phi=atan2(y-y_, x-x_);
	int n=findLayer(myConfig, z);
	int m=findLayer(myConfig, z_);
	double a, b;
	findMax(myConfig, &a, &b, rho, z, z_);
	complex double G1, G2;
	GfuncArgs myArgs={myConfig, z, z_, rho, a, b};
	G1 = QuadL_1D(GEM1, &myArgs, 0.0, 1.5*a, tol);
	G2 = QuadL_1D(GEM2, &myArgs, 0.0, 1.5*a, tol);
	complex double I=+0.5*G2+0.5*cos(2.0*phi)*G1;
	if (m==n){
		return I+GEM0yx(x, x_, y, y_, z, z_, myConfig->Layers[n].mu, 
			myConfig->Layers[n].eps, myConfig->k0);
	}else{
		return I;
	}
}

complex double GEMyy(Config *myConfig, double x, double x_, 
	double y,  double y_, double z, double z_, double tol){
	double rho=sqrt((x-x_)*(x-x_)+(y-y_)*(y-y_));
	double phi=atan2(y-y_, x-x_);
	double a, b;
	findMax(myConfig, &a, &b, rho, z, z_);
	complex double G1;
	GfuncArgs myArgs={myConfig, z, z_, rho, a, b};
	G1 = QuadL_1D(GEM1, &myArgs, 0.0, 1.5*a, tol);
	complex double I=+0.5*sin(2.0*phi)*G1;
	return I;
}

complex double GEMxz(Config *myConfig, double x, double x_, 
	double y,  double y_, double z, double z_, double tol){
	double rho=sqrt((x-x_)*(x-x_)+(y-y_)*(y-y_));
	double phi=atan2(y-y_, x-x_);
	int n=findLayer(myConfig, z);
	int m=findLayer(myConfig, z_);
	double a, b;
	findMax(myConfig, &a, &b, rho, z, z_);
	complex double G1;
	GfuncArgs myArgs={myConfig, z, z_, rho, a, b};
	G1 = QuadL_1D(GEM3, &myArgs, 0.0, 1.5*a, tol);
	complex double j=csqrt(-1.0);
	complex double mu_=myConfig->Layers[m].mu;
	double k0=myConfig->k0;
	complex double I=+(1.0/(j*eta0*k0*mu_))*sin(phi)*G1;
	if (m==n){
		return I+GEM0xz(x, x_, y, y_, z, z_, myConfig->Layers[n].mu, 
			myConfig->Layers[n].eps, myConfig->k0);
	}else{
		return I;
	}
}

complex double GEMyz(Config *myConfig, double x, double x_, 
	double y,  double y_, double z, double z_, double tol){
	double rho=sqrt((x-x_)*(x-x_)+(y-y_)*(y-y_));
	double phi=atan2(y-y_, x-x_);
	int n=findLayer(myConfig, z);
	int m=findLayer(myConfig, z_);
	double a, b;
	findMax(myConfig, &a, &b, rho, z, z_);
	complex double G1;
	GfuncArgs myArgs={myConfig, z, z_, rho, a, b};
	G1 = QuadL_1D(GEM3, &myArgs, 0.0, 1.5*a, tol);
	complex double j=csqrt(-1.0);
	complex double mu_=myConfig->Layers[m].mu;
	double k0=myConfig->k0;
	complex double I=-(1.0/(j*eta0*k0*mu_))*cos(phi)*G1;
	if (m==n){
		return I+GEM0yz(x, x_, y, y_, z, z_, myConfig->Layers[n].mu, 
			myConfig->Layers[n].eps, myConfig->k0);
	}else{
		return I;
	}
}

complex double GEMzx(Config *myConfig, double x, double x_, 
	double y,  double y_, double z, double z_, double tol){
	double rho=sqrt((x-x_)*(x-x_)+(y-y_)*(y-y_));
	double phi=atan2(y-y_, x-x_);
	int n=findLayer(myConfig, z);
	int m=findLayer(myConfig, z_);
	double a, b;
	findMax(myConfig, &a, &b, rho, z, z_);
	complex double G1;
	GfuncArgs myArgs={myConfig, z, z_, rho, a, b};
	G1 = QuadL_1D(GEM4, &myArgs, 0.0, 1.5*a, tol);
	complex double j=csqrt(-1.0);
	complex double eps=myConfig->Layers[n].eps;
	double k0=myConfig->k0;
	complex double I=(j*eta0/(k0*eps))*sin(phi)*G1;
	if (m==n){
		return I+GEM0zx(x, x_, y, y_, z, z_, myConfig->Layers[n].mu, 
			myConfig->Layers[n].eps, myConfig->k0);
	}else{
		return I;
	}
}

complex double GEMzy(Config *myConfig, double x, double x_, 
	double y,  double y_, double z, double z_, double tol){
	double rho=sqrt((x-x_)*(x-x_)+(y-y_)*(y-y_));
	double phi=atan2(y-y_, x-x_);
	int n=findLayer(myConfig, z);
	int m=findLayer(myConfig, z_);
	double a, b;
	findMax(myConfig, &a, &b, rho, z, z_);
	complex double G1;
	GfuncArgs myArgs={myConfig, z, z_, rho, a, b};
	G1 = QuadL_1D(GEM4, &myArgs, 0.0, 1.5*a, tol);
	complex double j=csqrt(-1.0);
	complex double eps=myConfig->Layers[n].eps;
	double k0=myConfig->k0;
	complex double I=-(j*eta0/(k0*eps))*cos(phi)*G1;
	if (m==n){
		return I+GEM0zy(x, x_, y, y_, z, z_, myConfig->Layers[n].mu, 
			myConfig->Layers[n].eps, myConfig->k0);
	}else{
		return I;
	}
}

complex double GEMzz(Config *myConfig, double x, double x_, 
	double y,  double y_, double z, double z_, double tol){
	double k0=myConfig->k0;
	return 0.0*x*x_*y*y_*z*z_*tol*k0;
}

complex double GHJ0xy(double x, double x_, double y, double y_,
	double z, double z_, complex double mu,complex double eps, double k0){
	complex double j=csqrt(-1.0);
	complex double k=k0*csqrt(mu*eps);
	double rho=sqrt((x-x_)*(x-x_)+(y-y_)*(y-y_));
	double r=sqrt(rho*rho+(z-z_)*(z-z_));
	complex double g=cexp(-j*k*r)/(4.0*pi*r);
	complex double dg1=-((1.0+j*k*r)/r)*g;
	if (r==0.0){
		return 0.0;
	}else{
		return -(z-z_)*dg1/r;
	}
}

complex double GHJ0xz(double x, double x_, double y, double y_,
	double z, double z_, complex double mu,complex double eps, double k0){
	complex double j=csqrt(-1.0);
	complex double k=k0*csqrt(mu*eps);
	double rho=sqrt((x-x_)*(x-x_)+(y-y_)*(y-y_));
	double r=sqrt(rho*rho+(z-z_)*(z-z_));
	complex double g=cexp(-j*k*r)/(4.0*pi*r);
	complex double dg1=-((1.0+j*k*r)/r)*g;
	if (r==0.0){
		return 0.0;
	}else{
		return +(y-y_)*dg1/r;
	}
}

complex double GHJ0yx(double x, double x_, double y, double y_,
	double z, double z_, complex double mu,complex double eps, double k0){
	complex double j=csqrt(-1.0);
	complex double k=k0*csqrt(mu*eps);
	double rho=sqrt((x-x_)*(x-x_)+(y-y_)*(y-y_));
	double r=sqrt(rho*rho+(z-z_)*(z-z_));
	complex double g=cexp(-j*k*r)/(4.0*pi*r);
	complex double dg1=-((1.0+j*k*r)/r)*g;
	if (r==0.0){
		return 0.0;
	}else{
		return +(z-z_)*dg1/r;
	}
}

complex double GHJ0yz(double x, double x_, double y, double y_,
	double z, double z_, complex double mu,complex double eps, double k0){
	complex double j=csqrt(-1.0);
	complex double k=k0*csqrt(mu*eps);
	double rho=sqrt((x-x_)*(x-x_)+(y-y_)*(y-y_));
	double r=sqrt(rho*rho+(z-z_)*(z-z_));
	complex double g=cexp(-j*k*r)/(4.0*pi*r);
	complex double dg1=-((1.0+j*k*r)/r)*g;
	if (r==0.0){
		return 0.0;
	}else{
		return -(x-x_)*dg1/r;
	}
}

complex double GHJ0zx(double x, double x_, double y, double y_,
	double z, double z_, complex double mu,complex double eps, double k0){
	complex double j=csqrt(-1.0);
	complex double k=k0*csqrt(mu*eps);
	double rho=sqrt((x-x_)*(x-x_)+(y-y_)*(y-y_));
	double r=sqrt(rho*rho+(z-z_)*(z-z_));
	complex double g=cexp(-j*k*r)/(4.0*pi*r);
	complex double dg1=-((1.0+j*k*r)/r)*g;
	if (r==0.0){
		return 0.0;
	}else{
		return -(y-y_)*dg1/r;
	}
}

complex double GHJ0zy(double x, double x_, double y, double y_,
	double z, double z_, complex double mu,complex double eps, double k0){
	complex double j=csqrt(-1.0);
	complex double k=k0*csqrt(mu*eps);
	double rho=sqrt((x-x_)*(x-x_)+(y-y_)*(y-y_));
	double r=sqrt(rho*rho+(z-z_)*(z-z_));
	complex double g=cexp(-j*k*r)/(4.0*pi*r);
	complex double dg1=-((1.0+j*k*r)/r)*g;
	if (r==0.0){
		return 0.0;
	}else{
		return +(x-x_)*dg1/r;
	}
}

complex double GHJ1(complex double k_rho_in, void *args){
	GfuncArgs *myArgs=(GfuncArgs*) args;
	Config *myConfig=myArgs->myConfig;
	double z=myArgs->z;
	double z_=myArgs->z_;
	double rho=myArgs->rho;
	double a=myArgs->a;
	double b=myArgs->b;
	complex double k_rho, dk_rho;
	spectralDetour(k_rho_in, &k_rho, &dk_rho, a, b);
	complex double V_h_i, I_h_i, V_e_i, I_e_i;
	TLGFr(myConfig, 'h', 'i', z, z_, k_rho, &V_h_i, &I_h_i);
	TLGFr(myConfig, 'e', 'i', z, z_, k_rho, &V_e_i, &I_e_i);
	return (I_h_i-I_e_i)*besselj(2, k_rho*rho)*k_rho*dk_rho/(2.0*pi);
}

complex double GHJ2(complex double k_rho_in, void *args){
	GfuncArgs *myArgs=(GfuncArgs*) args;
	Config *myConfig=myArgs->myConfig;
	double z=myArgs->z;
	double z_=myArgs->z_;
	double rho=myArgs->rho;
	double a=myArgs->a;
	double b=myArgs->b;
	complex double k_rho, dk_rho;
	spectralDetour(k_rho_in, &k_rho, &dk_rho, a, b);
	complex double V_h_i, I_h_i, V_e_i, I_e_i;
	TLGFr(myConfig, 'h', 'i', z, z_, k_rho, &V_h_i, &I_h_i);
	TLGFr(myConfig, 'e', 'i', z, z_, k_rho, &V_e_i, &I_e_i);
	return (I_h_i+I_e_i)*besselj(0, k_rho*rho)*k_rho*dk_rho/(2.0*pi);	
}

complex double GHJ3(complex double k_rho_in, void *args){
	GfuncArgs *myArgs=(GfuncArgs*) args;
	Config *myConfig=myArgs->myConfig;
	double z=myArgs->z;
	double z_=myArgs->z_;
	double rho=myArgs->rho;
	double a=myArgs->a;
	double b=myArgs->b;
	complex double k_rho, dk_rho;
	spectralDetour(k_rho_in, &k_rho, &dk_rho, a, b);
	complex double V_e_v, I_e_v;
	TLGFr(myConfig, 'e', 'v', z, z_, k_rho, &V_e_v, &I_e_v);
	return k_rho*I_e_v*besselj(1, k_rho*rho)*k_rho*dk_rho/(2.0*pi);	
}

complex double GHJ4(complex double k_rho_in, void *args){
	GfuncArgs *myArgs=(GfuncArgs*) args;
	Config *myConfig=myArgs->myConfig;
	double z=myArgs->z;
	double z_=myArgs->z_;
	double rho=myArgs->rho;
	double a=myArgs->a;
	double b=myArgs->b;
	complex double k_rho, dk_rho;
	spectralDetour(k_rho_in, &k_rho, &dk_rho, a, b);
	complex double V_h_i, I_h_i;
	TLGFr(myConfig, 'h', 'i', z, z_, k_rho, &V_h_i, &I_h_i);
	return k_rho*V_h_i*besselj(1, k_rho*rho)*k_rho*dk_rho/(2.0*pi);	
}

complex double GHJxx(Config *myConfig, double x, double x_, 
	double y,  double y_, double z, double z_, double tol){
	double rho=sqrt((x-x_)*(x-x_)+(y-y_)*(y-y_));
	double phi=atan2(y-y_, x-x_);
	double a, b;
	findMax(myConfig, &a, &b, rho, z, z_);
	complex double G1;
	GfuncArgs myArgs={myConfig, z, z_, rho, a, b};
	G1 = QuadL_1D(GHJ1, &myArgs, 0.0, 1.5*a, tol);
	complex double I=+0.5*sin(2.0*phi)*G1;
	return I;
}

complex double GHJxy(Config *myConfig, double x, double x_, 
	double y,  double y_, double z, double z_, double tol){
	double rho=sqrt((x-x_)*(x-x_)+(y-y_)*(y-y_));
	double phi=atan2(y-y_, x-x_);
	int n=findLayer(myConfig, z);
	int m=findLayer(myConfig, z_);
	double a, b;
	findMax(myConfig, &a, &b, rho, z, z_);
	complex double G1, G2;
	GfuncArgs myArgs={myConfig, z, z_, rho, a, b};
	G1 = QuadL_1D(GHJ1, &myArgs, 0.0, 1.5*a, tol);
	G2 = QuadL_1D(GHJ2, &myArgs, 0.0, 1.5*a, tol);
	complex double I=+0.5*G2-0.5*cos(2.0*phi)*G1;
	if (m==n){
		return I+GHJ0xy(x, x_, y, y_, z, z_, myConfig->Layers[n].mu, 
			myConfig->Layers[n].eps, myConfig->k0);
	}else{
		return I;
	}
}

complex double GHJyx(Config *myConfig, double x, double x_, 
	double y,  double y_, double z, double z_, double tol){
	double rho=sqrt((x-x_)*(x-x_)+(y-y_)*(y-y_));
	double phi=atan2(y-y_, x-x_);
	int n=findLayer(myConfig, z);
	int m=findLayer(myConfig, z_);
	double a, b;
	findMax(myConfig, &a, &b, rho, z, z_);
	complex double G1, G2;
	GfuncArgs myArgs={myConfig, z, z_, rho, a, b};
	G1 = QuadL_1D(GHJ1, &myArgs, 0.0, 1.5*a, tol);
	G2 = QuadL_1D(GHJ2, &myArgs, 0.0, 1.5*a, tol);
	complex double I=-0.5*G2-0.5*cos(2.0*phi)*G1;
	if (m==n){
		return I+GHJ0yx(x, x_, y, y_, z, z_, myConfig->Layers[n].mu, 
			myConfig->Layers[n].eps, myConfig->k0);
	}else{
		return I;
	}
}

complex double GHJyy(Config *myConfig, double x, double x_, 
	double y,  double y_, double z, double z_, double tol){
	double rho=sqrt((x-x_)*(x-x_)+(y-y_)*(y-y_));
	double phi=atan2(y-y_, x-x_);
	double a, b;
	findMax(myConfig, &a, &b, rho, z, z_);
	complex double G1;
	GfuncArgs myArgs={myConfig, z, z_, rho, a, b};
	G1 = QuadL_1D(GHJ1, &myArgs, 0.0, 1.5*a, tol);
	complex double I=-0.5*sin(2.0*phi)*G1;
	return I;
}

complex double GHJxz(Config *myConfig, double x, double x_, 
	double y,  double y_, double z, double z_, double tol){
	double rho=sqrt((x-x_)*(x-x_)+(y-y_)*(y-y_));
	double phi=atan2(y-y_, x-x_);
	int n=findLayer(myConfig, z);
	int m=findLayer(myConfig, z_);
	double a, b;
	findMax(myConfig, &a, &b, rho, z, z_);
	complex double G1;
	GfuncArgs myArgs={myConfig, z, z_, rho, a, b};
	G1 = QuadL_1D(GHJ3, &myArgs, 0.0, 1.5*a, tol);
	complex double j=csqrt(-1.0);
	complex double eps_=myConfig->Layers[m].eps;
	double k0=myConfig->k0;
	complex double I=-(eta0/(j*k0*eps_))*sin(phi)*G1;
	if (m==n){
		return I+GHJ0xz(x, x_, y, y_, z, z_, myConfig->Layers[n].mu, 
			myConfig->Layers[n].eps, myConfig->k0);
	}else{
		return I;
	}
}

complex double GHJyz(Config *myConfig, double x, double x_, 
	double y,  double y_, double z, double z_, double tol){
	double rho=sqrt((x-x_)*(x-x_)+(y-y_)*(y-y_));
	double phi=atan2(y-y_, x-x_);
	int n=findLayer(myConfig, z);
	int m=findLayer(myConfig, z_);
	double a, b;
	findMax(myConfig, &a, &b, rho, z, z_);
	complex double G1;
	GfuncArgs myArgs={myConfig, z, z_, rho, a, b};
	G1 = QuadL_1D(GHJ3, &myArgs, 0.0, 1.5*a, tol);
	complex double j=csqrt(-1.0);
	complex double eps_=myConfig->Layers[m].eps;
	double k0=myConfig->k0;
	complex double I=+(eta0/(j*k0*eps_))*cos(phi)*G1;
	if (m==n){
		return I+GHJ0yz(x, x_, y, y_, z, z_, myConfig->Layers[n].mu, 
			myConfig->Layers[n].eps, myConfig->k0);
	}else{
		return I;
	}
}

complex double GHJzx(Config *myConfig, double x, double x_, 
	double y,  double y_, double z, double z_, double tol){
	double rho=sqrt((x-x_)*(x-x_)+(y-y_)*(y-y_));
	double phi=atan2(y-y_, x-x_);
	int n=findLayer(myConfig, z);
	int m=findLayer(myConfig, z_);
	double a, b;
	findMax(myConfig, &a, &b, rho, z, z_);
	complex double G1;
	GfuncArgs myArgs={myConfig, z, z_, rho, a, b};
	G1 = QuadL_1D(GHJ4, &myArgs, 0.0, 1.5*a, tol);
	complex double j=csqrt(-1.0);
	complex double mu=myConfig->Layers[n].mu;
	double k0=myConfig->k0;
	complex double I=(1.0/(j*eta0*k0*mu))*sin(phi)*G1;
	if (m==n){
		return I+GHJ0zx(x, x_, y, y_, z, z_, myConfig->Layers[n].mu, 
			myConfig->Layers[n].eps, myConfig->k0);
	}else{
		return I;
	}
}

complex double GHJzy(Config *myConfig, double x, double x_, 
	double y,  double y_, double z, double z_, double tol){
	double rho=sqrt((x-x_)*(x-x_)+(y-y_)*(y-y_));
	double phi=atan2(y-y_, x-x_);
	int n=findLayer(myConfig, z);
	int m=findLayer(myConfig, z_);
	double a, b;
	findMax(myConfig, &a, &b, rho, z, z_);
	complex double G1;
	GfuncArgs myArgs={myConfig, z, z_, rho, a, b};
	G1 = QuadL_1D(GHJ4, &myArgs, 0.0, 1.5*a, tol);
	complex double j=csqrt(-1.0);
	complex double mu=myConfig->Layers[n].mu;
	double k0=myConfig->k0;
	complex double I=-(1.0/(j*eta0*k0*mu))*cos(phi)*G1;
	if (m==n){
		return I+GHJ0zy(x, x_, y, y_, z, z_, myConfig->Layers[n].mu, 
			myConfig->Layers[n].eps, myConfig->k0);
	}else{
		return I;
	}
}

complex double GHJzz(Config *myConfig, double x, double x_, 
	double y,  double y_, double z, double z_, double tol){
	double k0=myConfig->k0;
	return 0.0*x*x_*y*y_*z*z_*tol*k0;
}

complex double GHM0xx(double x, double x_, double y, double y_,
	double z, double z_, double omega, complex double mu, 
	complex double eps, double k0){
	complex double j=csqrt(-1.0);
	complex double k=k0*csqrt(mu*eps);
	double rho=sqrt((x-x_)*(x-x_)+(y-y_)*(y-y_));
	double r=sqrt(rho*rho+(z-z_)*(z-z_));
	complex double g=cexp(-j*k*r)/(4.0*pi*r);
	complex double dg1=-((1.0+j*k*r)/r)*g;
	complex double dg2=-((1.0+j*k*r)/r)*dg1+g/(r*r);
	if (r==0.0){
		return -j*omega*eps0*eps*(-1.0/(3.0*k*k));
	}else{
		return -j*omega*eps0*eps*(g+(1.0/(k*k))*(dg1/r+((x-x_)*(x-x_)/(r*r))*(dg2-dg1/r)));
	}
}

complex double GHM0xy(double x, double x_, double y, double y_,
	double z, double z_, double omega, complex double mu, 
	complex double eps, double k0){
	complex double j=csqrt(-1.0);
	complex double k=k0*csqrt(mu*eps);
	double rho=sqrt((x-x_)*(x-x_)+(y-y_)*(y-y_));
	double r=sqrt(rho*rho+(z-z_)*(z-z_));
	complex double g=cexp(-j*k*r)/(4.0*pi*r);
	complex double dg1=-((1.0+j*k*r)/r)*g;
	complex double dg2=-((1.0+j*k*r)/r)*dg1+g/(r*r);
	if (r==0.0){
		return 0.0;
	}else{
		return -j*omega*eps0*eps*((1.0/(k*k))*(((x-x_)*(y-y_)/(r*r))*(dg2-dg1/r)));
	}
}

complex double GHM0xz(double x, double x_, double y, double y_,
	double z, double z_, double omega, complex double mu, 
	complex double eps, double k0){
	complex double j=csqrt(-1.0);
	complex double k=k0*csqrt(mu*eps);
	double rho=sqrt((x-x_)*(x-x_)+(y-y_)*(y-y_));
	double r=sqrt(rho*rho+(z-z_)*(z-z_));
	complex double g=cexp(-j*k*r)/(4.0*pi*r);
	complex double dg1=-((1.0+j*k*r)/r)*g;
	complex double dg2=-((1.0+j*k*r)/r)*dg1+g/(r*r);
	if (r==0.0){
		return 0.0;
	}else{
		return -j*omega*eps0*eps*((1.0/(k*k))*(((x-x_)*(z-z_)/(r*r))*(dg2-dg1/r)));
	}
}

complex double GHM0yx(double x, double x_, double y, double y_,
	double z, double z_, double omega, complex double mu, 
	complex double eps, double k0){
	complex double j=csqrt(-1.0);
	complex double k=k0*csqrt(mu*eps);
	double rho=sqrt((x-x_)*(x-x_)+(y-y_)*(y-y_));
	double r=sqrt(rho*rho+(z-z_)*(z-z_));
	complex double g=cexp(-j*k*r)/(4.0*pi*r);
	complex double dg1=-((1.0+j*k*r)/r)*g;
	complex double dg2=-((1.0+j*k*r)/r)*dg1+g/(r*r);
	if (r==0.0){
		return 0.0;
	}else{
		return -j*omega*eps0*eps*((1.0/(k*k))*(((y-y_)*(x-x_)/(r*r))*(dg2-dg1/r)));
	}
}

complex double GHM0yy(double x, double x_, double y, double y_,
	double z, double z_, double omega, complex double mu, 
	complex double eps, double k0){
	complex double j=csqrt(-1.0);
	complex double k=k0*csqrt(mu*eps);
	double rho=sqrt((x-x_)*(x-x_)+(y-y_)*(y-y_));
	double r=sqrt(rho*rho+(z-z_)*(z-z_));
	complex double g=cexp(-j*k*r)/(4.0*pi*r);
	complex double dg1=-((1.0+j*k*r)/r)*g;
	complex double dg2=-((1.0+j*k*r)/r)*dg1+g/(r*r);
	if (r==0.0){
		return -j*omega*eps0*eps*(-1.0/(3.0*k*k));
	}else{
		return -j*omega*eps0*eps*(g+(1.0/(k*k))*(dg1/r+((y-y_)*(y-y_)/(r*r))*(dg2-dg1/r)));
	}
}

complex double GHM0yz(double x, double x_, double y, double y_,
	double z, double z_, double omega, complex double mu, 
	complex double eps, double k0){
	complex double j=csqrt(-1.0);
	complex double k=k0*csqrt(mu*eps);
	double rho=sqrt((x-x_)*(x-x_)+(y-y_)*(y-y_));
	double r=sqrt(rho*rho+(z-z_)*(z-z_));
	complex double g=cexp(-j*k*r)/(4.0*pi*r);
	complex double dg1=-((1.0+j*k*r)/r)*g;
	complex double dg2=-((1.0+j*k*r)/r)*dg1+g/(r*r);
	if (r==0.0){
		return 0.0;
	}else{
		return -j*omega*eps0*eps*((1.0/(k*k))*(((y-y_)*(z-z_)/(r*r))*(dg2-dg1/r)));
	}
}

complex double GHM0zx(double x, double x_, double y, double y_,
	double z, double z_, double omega, complex double mu, 
	complex double eps, double k0){
	complex double j=csqrt(-1.0);
	complex double k=k0*csqrt(mu*eps);
	double rho=sqrt((x-x_)*(x-x_)+(y-y_)*(y-y_));
	double r=sqrt(rho*rho+(z-z_)*(z-z_));
	complex double g=cexp(-j*k*r)/(4.0*pi*r);
	complex double dg1=-((1.0+j*k*r)/r)*g;
	complex double dg2=-((1.0+j*k*r)/r)*dg1+g/(r*r);
	if (r==0.0){
		return 0.0;
	}else{
		return -j*omega*eps0*eps*((1.0/(k*k))*(((z-z_)*(x-x_)/(r*r))*(dg2-dg1/r)));
	}
}

complex double GHM0zy(double x, double x_, double y, double y_,
	double z, double z_, double omega, complex double mu, 
	complex double eps, double k0){
	complex double j=csqrt(-1.0);
	complex double k=k0*csqrt(mu*eps);
	double rho=sqrt((x-x_)*(x-x_)+(y-y_)*(y-y_));
	double r=sqrt(rho*rho+(z-z_)*(z-z_));
	complex double g=cexp(-j*k*r)/(4.0*pi*r);
	complex double dg1=-((1.0+j*k*r)/r)*g;
	complex double dg2=-((1.0+j*k*r)/r)*dg1+g/(r*r);
	if (r==0.0){
		return 0.0;
	}else{
		return -j*omega*eps0*eps*((1.0/(k*k))*(((z-z_)*(y-y_)/(r*r))*(dg2-dg1/r)));
	}
}

complex double GHM0zz(double x, double x_, double y, double y_,
	double z, double z_, double omega, complex double mu, 
	complex double eps, double k0){
	complex double j=csqrt(-1.0);
	complex double k=k0*csqrt(mu*eps);
	double rho=sqrt((x-x_)*(x-x_)+(y-y_)*(y-y_));
	double r=sqrt(rho*rho+(z-z_)*(z-z_));
	complex double g=cexp(-j*k*r)/(4.0*pi*r);
	complex double dg1=-((1.0+j*k*r)/r)*g;
	complex double dg2=-((1.0+j*k*r)/r)*dg1+g/(r*r);
	if (r==0.0){
		return -j*omega*eps0*eps*(-1.0/(3.0*k*k));
	}else{
		return -j*omega*eps0*eps*(g+(1.0/(k*k))*(dg1/r+((z-z_)*(z-z_)/(r*r))*(dg2-dg1/r)));
	}
}

complex double GHM1(complex double k_rho_in, void *args){
	GfuncArgs *myArgs=(GfuncArgs*) args;
	Config *myConfig=myArgs->myConfig;
	double z=myArgs->z;
	double z_=myArgs->z_;
	double rho=myArgs->rho;
	double a=myArgs->a;
	double b=myArgs->b;
	complex double k_rho, dk_rho;
	spectralDetour(k_rho_in, &k_rho, &dk_rho, a, b);
	complex double V_h_v, I_h_v, V_e_v, I_e_v;
	TLGFr(myConfig, 'h', 'v', z, z_, k_rho, &V_h_v, &I_h_v);
	TLGFr(myConfig, 'e', 'v', z, z_, k_rho, &V_e_v, &I_e_v);
	return (I_h_v+I_e_v)*besselj(0, k_rho*rho)*k_rho*dk_rho/(2.0*pi);
}

complex double GHM2(complex double k_rho_in, void *args){
	GfuncArgs *myArgs=(GfuncArgs*) args;
	Config *myConfig=myArgs->myConfig;
	double z=myArgs->z;
	double z_=myArgs->z_;
	double rho=myArgs->rho;
	double a=myArgs->a;
	double b=myArgs->b;
	complex double k_rho, dk_rho;
	spectralDetour(k_rho_in, &k_rho, &dk_rho, a, b);
	complex double V_h_v, I_h_v, V_e_v, I_e_v;
	TLGFr(myConfig, 'h', 'v', z, z_, k_rho, &V_h_v, &I_h_v);
	TLGFr(myConfig, 'e', 'v', z, z_, k_rho, &V_e_v, &I_e_v);
	return (I_h_v-I_e_v)*besselj(2, k_rho*rho)*k_rho*dk_rho/(2.0*pi);	
}

complex double GHM3(complex double k_rho_in, void *args){
	GfuncArgs *myArgs=(GfuncArgs*) args;
	Config *myConfig=myArgs->myConfig;
	double z=myArgs->z;
	double z_=myArgs->z_;
	double rho=myArgs->rho;
	double a=myArgs->a;
	double b=myArgs->b;
	complex double k_rho, dk_rho;
	spectralDetour(k_rho_in, &k_rho, &dk_rho, a, b);
	complex double V_h_i, I_h_i;
	TLGFr(myConfig, 'h', 'i', z, z_, k_rho, &V_h_i, &I_h_i);
	return k_rho*I_h_i*besselj(1, k_rho*rho)*k_rho*dk_rho/(2.0*pi);	
}

complex double GHM4(complex double k_rho_in, void *args){
	GfuncArgs *myArgs=(GfuncArgs*) args;
	Config *myConfig=myArgs->myConfig;
	double z=myArgs->z;
	double z_=myArgs->z_;
	double rho=myArgs->rho;
	double a=myArgs->a;
	double b=myArgs->b;
	complex double k_rho, dk_rho;
	spectralDetour(k_rho_in, &k_rho, &dk_rho, a, b);
	complex double V_h_v, I_h_v;
	TLGFr(myConfig, 'h', 'v', z, z_, k_rho, &V_h_v, &I_h_v);
	return k_rho*V_h_v*besselj(1, k_rho*rho)*k_rho*dk_rho/(2.0*pi);	
}

complex double GHM5(complex double k_rho_in, void *args){
	GfuncArgs *myArgs=(GfuncArgs*) args;
	Config *myConfig=myArgs->myConfig;
	double z=myArgs->z;
	double z_=myArgs->z_;
	double rho=myArgs->rho;
	double a=myArgs->a;
	double b=myArgs->b;
	complex double k_rho, dk_rho;
	spectralDetour(k_rho_in, &k_rho, &dk_rho, a, b);
	complex double V_h_i, I_h_i;
	TLGFr(myConfig, 'h', 'i', z, z_, k_rho, &V_h_i, &I_h_i);
	return k_rho*k_rho*V_h_i*besselj(0, k_rho*rho)*k_rho*dk_rho/(2.0*pi);	
}

complex double GHMxx(Config *myConfig, double x, double x_, 
	double y,  double y_, double z, double z_, double tol){
	double rho=sqrt((x-x_)*(x-x_)+(y-y_)*(y-y_));
	double phi=atan2(y-y_, x-x_);
	int n=findLayer(myConfig, z);
	int m=findLayer(myConfig, z_);
	double a, b;
	findMax(myConfig, &a, &b, rho, z, z_);
	complex double G1, G2;
	GfuncArgs myArgs={myConfig, z, z_, rho, a, b};
	G1 = QuadL_1D(GHM1, &myArgs, 0.0, 1.5*a, tol);
	G2 = QuadL_1D(GHM2, &myArgs, 0.0, 1.5*a, tol);
	complex double I=-0.5*G1+0.5*cos(2.0*phi)*G2;
	if (m==n){
		return I+GHM0xx(x, x_, y, y_, z, z_,
			myConfig->omega, myConfig->Layers[n].mu, 
			myConfig->Layers[n].eps, myConfig->k0);
	}else{
		return I;
	}
}

complex double GHMxy(Config *myConfig, double x, double x_, 
	double y,  double y_, double z, double z_, double tol){
	double rho=sqrt((x-x_)*(x-x_)+(y-y_)*(y-y_));
	double phi=atan2(y-y_, x-x_);
	int n=findLayer(myConfig, z);
	int m=findLayer(myConfig, z_);
	double a, b;
	findMax(myConfig, &a, &b, rho, z, z_);
	complex double G1;
	GfuncArgs myArgs={myConfig, z, z_, rho, a, b};
	G1 = QuadL_1D(GHM2, &myArgs, 0.0, 1.5*a, tol);
	complex double I=0.5*sin(2.0*phi)*G1;
	if (m==n){
		return I+GHM0xy(x, x_, y, y_, z, z_,
			myConfig->omega, myConfig->Layers[n].mu, 
			myConfig->Layers[n].eps, myConfig->k0);
	}else{
		return I;
	}
}

complex double GHMyx(Config *myConfig, double x, double x_, 
	double y,  double y_, double z, double z_, double tol){
	double rho=sqrt((x-x_)*(x-x_)+(y-y_)*(y-y_));
	double phi=atan2(y-y_, x-x_);
	int n=findLayer(myConfig, z);
	int m=findLayer(myConfig, z_);
	double a, b;
	findMax(myConfig, &a, &b, rho, z, z_);
	complex double G1;
	GfuncArgs myArgs={myConfig, z, z_, rho, a, b};
	G1 = QuadL_1D(GHM2, &myArgs, 0.0, 1.5*a, tol);
	complex double I=0.5*sin(2.0*phi)*G1;
	if (m==n){
		return I+GHM0yx(x, x_, y, y_, z, z_,
			myConfig->omega, myConfig->Layers[n].mu, 
			myConfig->Layers[n].eps, myConfig->k0);
	}else{
		return I;
	}
}

complex double GHMyy(Config *myConfig, double x, double x_, 
	double y,  double y_, double z, double z_, double tol){
	double rho=sqrt((x-x_)*(x-x_)+(y-y_)*(y-y_));
	double phi=atan2(y-y_, x-x_);
	int n=findLayer(myConfig, z);
	int m=findLayer(myConfig, z_);
	double a, b;
	findMax(myConfig, &a, &b, rho, z, z_);
	complex double G1, G2;
	GfuncArgs myArgs={myConfig, z, z_, rho, a, b};
	G1 = QuadL_1D(GHM1, &myArgs, 0.0, 1.5*a, tol);
	G2 = QuadL_1D(GHM2, &myArgs, 0.0, 1.5*a, tol);
	complex double I=-0.5*G1-0.5*cos(2.0*phi)*G2;
	if (m==n){
		return I+GHM0yy(x, x_, y, y_, z, z_,
			myConfig->omega, myConfig->Layers[n].mu, 
			myConfig->Layers[n].eps, myConfig->k0);
	}else{
		return I;
	}
}

complex double GHMxz(Config *myConfig, double x, double x_, 
	double y,  double y_, double z, double z_, double tol){
	double rho=sqrt((x-x_)*(x-x_)+(y-y_)*(y-y_));
	double phi=atan2(y-y_, x-x_);
	int n=findLayer(myConfig, z);
	int m=findLayer(myConfig, z_);
	double a, b;
	findMax(myConfig, &a, &b, rho, z, z_);
	complex double G1;
	GfuncArgs myArgs={myConfig, z, z_, rho, a, b};
	G1 = QuadL_1D(GHM3, &myArgs, 0.0, 1.5*a, tol);
	complex double j=csqrt(-1.0);
	complex double mu_=myConfig->Layers[m].mu;
	double k0=myConfig->k0;
	complex double I=(1.0/(j*eta0*k0*mu_))*cos(phi)*G1;
	if (m==n){
		return I+GHM0xz(x, x_, y, y_, z, z_,
			myConfig->omega, myConfig->Layers[n].mu, 
			myConfig->Layers[n].eps, myConfig->k0);
	}else{
		return I;
	}
}

complex double GHMyz(Config *myConfig, double x, double x_, 
	double y,  double y_, double z, double z_, double tol){
	double rho=sqrt((x-x_)*(x-x_)+(y-y_)*(y-y_));
	double phi=atan2(y-y_, x-x_);
	int n=findLayer(myConfig, z);
	int m=findLayer(myConfig, z_);
	double a, b;
	findMax(myConfig, &a, &b, rho, z, z_);
	complex double G1;
	GfuncArgs myArgs={myConfig, z, z_, rho, a, b};
	G1 = QuadL_1D(GHM3, &myArgs, 0.0, 1.5*a, tol);
	complex double j=csqrt(-1.0);
	complex double mu_=myConfig->Layers[m].mu;
	double k0=myConfig->k0;
	complex double I=(1.0/(j*eta0*k0*mu_))*sin(phi)*G1;
	if (m==n){
		return I+GHM0yz(x, x_, y, y_, z, z_,
			myConfig->omega, myConfig->Layers[n].mu, 
			myConfig->Layers[n].eps, myConfig->k0);
	}else{
		return I;
	}
}

complex double GHMzx(Config *myConfig, double x, double x_, 
	double y,  double y_, double z, double z_, double tol){
	double rho=sqrt((x-x_)*(x-x_)+(y-y_)*(y-y_));
	double phi=atan2(y-y_, x-x_);
	int n=findLayer(myConfig, z);
	int m=findLayer(myConfig, z_);
	double a, b;
	findMax(myConfig, &a, &b, rho, z, z_);
	complex double G1;
	GfuncArgs myArgs={myConfig, z, z_, rho, a, b};
	G1 = QuadL_1D(GHM4, &myArgs, 0.0, 1.5*a, tol);
	complex double j=csqrt(-1.0);
	complex double mu=myConfig->Layers[n].mu;
	double k0=myConfig->k0;
	complex double I=(1.0/(j*eta0*k0*mu))*cos(phi)*G1;
	if (m==n){
		return I+GHM0zx(x, x_, y, y_, z, z_,
			myConfig->omega, myConfig->Layers[n].mu, 
			myConfig->Layers[n].eps, myConfig->k0);
	}else{
		return I;
	}
}

complex double GHMzy(Config *myConfig, double x, double x_, 
	double y,  double y_, double z, double z_, double tol){
	double rho=sqrt((x-x_)*(x-x_)+(y-y_)*(y-y_));
	double phi=atan2(y-y_, x-x_);
	int n=findLayer(myConfig, z);
	int m=findLayer(myConfig, z_);
	double a, b;
	findMax(myConfig, &a, &b, rho, z, z_);
	complex double G1;
	GfuncArgs myArgs={myConfig, z, z_, rho, a, b};
	G1 = QuadL_1D(GHM4, &myArgs, 0.0, 1.5*a, tol);
	complex double j=csqrt(-1.0);
	complex double mu=myConfig->Layers[n].mu;
	double k0=myConfig->k0;
	complex double I=(1.0/(j*eta0*k0*mu))*sin(phi)*G1;
	if (m==n){
		return I+GHM0zy(x, x_, y, y_, z, z_,
			myConfig->omega, myConfig->Layers[n].mu, 
			myConfig->Layers[n].eps, myConfig->k0);
	}else{
		return I;
	}
}

complex double GHMzz(Config *myConfig, double x, double x_, 
	double y,  double y_, double z, double z_, double tol){
	double rho=sqrt((x-x_)*(x-x_)+(y-y_)*(y-y_));
	int n=findLayer(myConfig, z);
	int m=findLayer(myConfig, z_);
	double a, b;
	findMax(myConfig, &a, &b, rho, z, z_);
	complex double G1;
	GfuncArgs myArgs={myConfig, z, z_, rho, a, b};
	G1 = QuadL_1D(GHM5, &myArgs, 0.0, 1.5*a, tol);
	complex double j=csqrt(-1.0);
	complex double mu_=myConfig->Layers[m].mu;
	complex double mu=myConfig->Layers[n].mu;
	double k0=myConfig->k0;
	complex double I=-(1.0/(eta0*eta0*k0*k0*mu_*mu))*G1;
	if (rho==0.0&&fabs(z-z_)==0.0){
		I-=(1.0/(j*eta0*k0*mu));
	}
	if (m==n){
		return I+GHM0zz(x, x_, y, y_, z, z_,
			myConfig->omega, myConfig->Layers[n].mu, 
			myConfig->Layers[n].eps, myConfig->k0);
	}else{
		return I;
	}
}