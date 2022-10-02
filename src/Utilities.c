//
#include "../include/Utilities.h"

void printComplex(complex double z){
	printf("(%21.14E, %21.14E)\n", creal(z), cimag(z));
}

double roundn(double x, int n){
	double r=pow(10.0, n);
	return round(x*r)/r;
}

void ticTimer(Timer *T){
	T->start = clock();
	T->flag = 1;
}

void tocTimer(Timer *T){
	if (T->flag!=1){
		printf("Timer was not set!\n");
	}else{
		T->stop = clock();
		T->elapsed=((T->stop)-(T->start))/CLOCKS_PER_SEC;
		if (T->elapsed<60.0){
			printf("Elapsed time is %0.2f seconds\n", T->elapsed);
		}else
		if (T->elapsed<3600.0){
			printf("Elapsed time is %0.2f minutes\n", T->elapsed/60.0);
		}else{
			printf("Elapsed time is %0.2f hours\n", T->elapsed/3600.0);
		}
		T->flag = 0;
	}
}

void progressBar(int n, int N){
	n++;
	int NBar=30;
	printf("\rProgress: [");
	for (int i=0; i<n*NBar/N; i++){
		printf("#");
	}
	for (int i=n*NBar/N; i<NBar; i++){
		printf("-");
	}
	printf("] %3d%%", 100*n/N);
	if (n==N){
		printf(" Done!\n");
	}
}

double deg2rad(double theta){
	return theta*pi/180.0;
}

double rad2deg(double theta){
	return theta*180.0/pi;
}