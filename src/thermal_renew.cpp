/*
 *  thermal_renew.cpp
 *  thermal1
 *
 *  Created by Victor Sacek on 15/03/10.
 *  Copyright 2010 IAG/USP. All rights reserved.
 *
 */

#include <math.h>
#include <stdio.h>

extern double *T_vec;
extern double *T_vec_fut;

extern double dt_calor_sec;
extern double seg_per_ano;

extern long nodes_thermal;
extern long *cond_thermal;

extern double **xyz_thermal;
extern double **xyz_thermal_fut;

extern double deltaz;

extern double kappa;

extern long nodes_flex;
extern long layers;

extern double *h_top;
extern double h_bot;

extern double *temper_q;
extern double *temper_q_a;
extern double alpha_exp_thermo;

extern double RHOM;

extern long nodes;

extern double **xy_flex;
extern double *h_q;

extern double minx;
extern double maxx;

extern double **peso;
extern long **peso_pos;

void thermal_renew()
{
	long i,j;
	double v_prov = 500.0/1.0E6/seg_per_ano;
	long p0,p1;
	double z0,z1,zz;
	double T0,T1,theta,A0,B0;
	long estranho=0;
	
	double h_aux,T_aux;
	
	for (i=0;i<nodes_flex;i++){
		temper_q_a[i]=0.0;
	}
	
	for (j=0;j<layers-1;j++){
		for (i=0;i<nodes_flex;i++){
			h_aux = xyz_thermal[j*nodes_flex+i][2]-xyz_thermal[(j+1)*nodes_flex+i][2];
			T_aux = (T_vec[j*nodes_flex+i]+T_vec[(j+1)*nodes_flex+i])/2;
			temper_q_a[i]+=h_aux*9.8*RHOM*(1.0-alpha_exp_thermo*T_aux);
			
			h_aux = xyz_thermal_fut[j*nodes_flex+i][2]-xyz_thermal_fut[(j+1)*nodes_flex+i][2];
			T_aux = (T_vec_fut[j*nodes_flex+i]+T_vec_fut[(j+1)*nodes_flex+i])/2;
			temper_q_a[i]-=h_aux*9.8*RHOM*(1.0-alpha_exp_thermo*T_aux);
		}
	}
	
	
	for (i=0;i<nodes_flex;i++){
		if (xy_flex[i][0]<minx || xy_flex[i][0]>maxx){
			temper_q[i]+=temper_q_a[i];
		}
	}
	
	
	
	for (j=0;j<layers;j++){
		for (i=0;i<nodes_flex;i++){
			xyz_thermal[j*nodes_flex+i][2]=h_top[i] + j*(h_bot-h_top[i])/(layers-1);
		}
	}
	
	for (i=0;i<nodes_thermal;i++){
		if (cond_thermal[i]==0){
			
			if (xyz_thermal[i][2]>xyz_thermal_fut[i][2]){
				p0 = i-nodes_flex;
				p1 = i;
			}
			else {
				p0 = i;
				p1 = i+nodes_flex;
			}
			
			//T_vec[i]=T_vec_fut[i]*(ll-v_prov*dt_calor_sec)/ll+T_vec_fut[i+36]*v_prov*dt_calor_sec/ll;
			
			z0 = xyz_thermal_fut[p0][2];
			z1 = xyz_thermal_fut[p1][2];
			zz = xyz_thermal[i][2]-z0;
			
			theta=exp((z1-z0)*v_prov/kappa);
			
			T0=T_vec_fut[p0];
			T1=T_vec_fut[p1];
			
			A0 = (T1-T0)/(theta-1);
			B0 = T0-A0;
			T_vec[i]=A0*exp(zz*v_prov/kappa)+B0;
			
			if (T_vec[i]>0);
			else estranho=1;
		}
	}
	printf("T_vec:%ld\n",estranho);

}