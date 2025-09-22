/*
 *  thermal_interp.cpp
 *  Thermal_hist3
 *
 *  Created by Victor Sacek on 23/03/10.
 *  Copyright 2010 IAG/USP. All rights reserved.
 *
 */
#include <math.h>
#include <stdio.h>

extern long layers;
extern double **xyz_thermal;
extern double *T_vec;
extern long nodes_thermal;
extern long nodes_flex;

extern double kappa;
extern double seg_per_ano;

double thermal_interp(long p, double TT)
{
	long i,p0,p1,cond=0;
	double zz=0;
	double z0,z1;
	//printf("teper");
	
	for (i=p;i<nodes_thermal-nodes_flex && cond==0;i+=nodes_flex){
		if (TT>T_vec[i] && TT<=T_vec[i+nodes_flex]){
			p0=i;
			p1=i+nodes_flex;
			cond=1;
		}
	}
	
	if (cond==1){
		z1 = xyz_thermal[p1][2];
		z0 = xyz_thermal[p0][2];
		double v_prov = 500.0/1.0E6/seg_per_ano;
		double T1=T_vec[p1], T0=T_vec[p0];
		double theta=exp((z1-z0)*v_prov/kappa);
		double A0 = (T1-T0)/(theta-1),B0 = T0-A0;
		
		zz = log((TT-B0)/A0)*kappa/v_prov+z0;
		
	}
	return (zz);
	
}
