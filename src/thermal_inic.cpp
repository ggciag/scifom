/*
 *  thermal_inic.cpp
 *  thermal1
 *
 *  Created by Victor Sacek on 15/03/10.
 *  Copyright 2010 IAG/USP. All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>

void thermal_solv(double **Kf, double *Kd, double *b_vec, double *x_vec);

extern double **Kthermal1;
extern double *Kthermal1_diag;
extern double *fthermal;

extern long nodes_thermal;
extern double **xyz_thermal;
extern double *T_vec;
extern double *T_vec_fut;
extern double *T_vec0;

extern double **xy_flex;
extern double maxy;


extern long layers;
extern long nodes_flex;

void thermal_inic()
{
	
	long i,j;
	for (j=0;j<layers;j++){
		for (i=0;i<nodes_flex;i++){
			T_vec[j*nodes_flex+i]=j*1300.0/(layers-1);
		}
	}
	
	/*j=layers-1;
	for (i=0;i<nodes_flex;i++){
		if (xy_flex[i][0]>600000.0 && xy_flex[i][0]<800000.0 && xy_flex[i][1]<maxy/2){
			T_vec[j*nodes_flex+i]=1500.0;
		}
	}*/
	
	thermal_solv(Kthermal1,Kthermal1_diag,fthermal,T_vec);
	
	for (j=0;j<layers;j++){
		for (i=0;i<nodes_flex;i++){
			T_vec_fut[j*nodes_flex+i]=T_vec[j*nodes_flex+i];
			T_vec0[j*nodes_flex+i]=T_vec[j*nodes_flex+i];
		}
	}
	
	
}