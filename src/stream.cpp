#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double var_fault_stream(double v);

extern long nodes_thermal;
extern double **xyz_thermal;
extern double **v_adv;

extern double seg_per_ano;
extern double V_meio;

extern double axis_stream;

double stream_x (double x,double z);
double stream_z (double x,double z);

extern long cont_falha;
extern long num_falha;
extern double **pos_falha_vec;

extern double tempo;

void stream()
{
    long i;
    double x,z;
    
    for (i=0;i<nodes_thermal;i++){
		z = xyz_thermal[i][2];
		x = xyz_thermal[i][0]-(axis_stream);
		v_adv[i][0]=stream_x(x,z);        
		v_adv[i][2]=stream_z(x,z);		
	}
	
}



double stream_x (double x,double z){
	
	double A,B,C,D;
    double pi=3.14159265;
    double U=V_meio/seg_per_ano;
    double Uz = U;
	
	double v_x;
    
    A = 0;
    B = (2*pi*Uz-pi*pi*U)/(pi*pi-4);
    C = (4*U-2*pi*Uz)/(pi*pi-4);
    D = (2*pi*U-4*Uz)/(pi*pi-4);
	
	double x0,z0;
	
	if (x==0.0) x0 = 0.000001;
	else x0=fabs(x);
	if (z==0) z0 = 0.000001;
	else z0=-z;
	
	if (x==0 && z==0){
		v_x = 0.0;
	}
	else {
		v_x = -B -D*(atan(z0/x0))+ (C*x0+D*z0)*(-x0/(x0*x0+z0*z0));
	}
	if (x<0) v_x*=-1.0;
	
	return(v_x+V_meio/seg_per_ano);
}



double stream_z (double x,double z){
	
	double A,B,C,D;
    double pi=3.14159265;
    double U=V_meio/seg_per_ano;
    double Uz = U;
	
	double v_z;
    
    A = 0;
    B = (2*pi*Uz-pi*pi*U)/(pi*pi-4);
    C = (4*U-2*pi*Uz)/(pi*pi-4);
    D = (2*pi*U-4*Uz)/(pi*pi-4);
	
	double x0,z0;
	
	if (x==0.0) x0 = 0.000001;
	else x0=fabs(x);
	if (z==0) z0 = 0.000001;
	else z0=-z;
	
	if (x==0 && z==0){
		v_z = 0.0;
	}
	else {
		v_z = -(A +C*(atan(z0/x0))+ (C*x0+D*z0)*(-z0/(x0*x0+z0*z0)));
	}
	
	return(v_z);
}