/*
 *  fault3D.cpp
 *  Flex_Erod_2_6_1
 *
 *  Created by Victor Sacek on 11/06/12.
 *  Copyright 2012 IAG/USP. All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

extern double *h_topo;
extern double **xy;

extern double *moho;
extern double *h_crust_sup;
extern double *h_flex_cumulat;

extern double *h_bed;

extern double RHOC;

extern double dt;

extern long nodes;

extern long n_falhas;
extern double **falhas_param;

extern double tempo;

void fault3D()
{
	double g = 9.8;

	double w_min;
	double w_max;
	double Larg;
	double transX;
	double transY;
	double theta;	
	
	double Te_uc;
	double nu = 0.25;
	double Elast = 1E11;
	double D;
	
	long n;
	
	for (n=0;n<n_falhas;n++){
		if (tempo>=falhas_param[n][7] && tempo<=falhas_param[n][8]){
			w_min = falhas_param[n][0]*dt;
			w_max = falhas_param[n][1]*dt;
			Larg = falhas_param[n][2];
			transX = falhas_param[n][3];
			transY = falhas_param[n][4];
			theta = falhas_param[n][5]*3.14159/180.;
			Te_uc = falhas_param[n][6];
			
			D = Elast*pow(Te_uc,3.0)/(12*(1-nu*nu));
			
			double alpha = pow(D/(RHOC*g),0.25);
			
			double x,y,x1,y1;
			
			double h_aux;
			
			for (long i=0;i<nodes;i++){
				x = xy[i][0];
				y = xy[i][1];
				
				x = x-transX;
				y = y-transY;
				
				x1 = x*cos(theta)-y*sin(theta);
				y1 = x*sin(theta)+y*cos(theta);
				
				x=x1;
				y=y1;
				
				if (x>=0){
					h_aux=w_min*exp(-0.701*x/alpha)*cos(0.701*x/alpha)*exp(-pow(y/Larg,4));
				}
				else{
					x=-x;
					h_aux=w_max*exp(-0.701*x/alpha)*cos(0.701*x/alpha)*exp(-pow(y/Larg,4));
				}
				
				h_topo[i]+=h_aux;
				h_bed[i]+=h_aux;
				h_crust_sup[i]+=h_aux;
				if (moho[i]>h_crust_sup[i]) moho[i]=h_crust_sup[i];
				h_flex_cumulat[i]+=h_aux;
				
			}
		}
	}
}