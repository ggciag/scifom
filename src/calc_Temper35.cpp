/*
 *  calc_Temper35.cpp
 *  FlexErod2_7
 *
 *  Created by Victor Sacek on 15/06/12.
 *  Copyright 2012 IAG/USP. All rights reserved.
 *
 */


extern long nodes;

extern double *Temper35;

extern double *h_temper35;

extern long layers;
extern double **xyz_thermal;
extern double *T_vec;
extern long nodes_thermal;
extern long nodes_flex;

extern double kappa;
extern double seg_per_ano;


extern long **peso_pos;
extern double **peso;


void calc_Temper35()
{
	
	long i,p0,p1,cond=0;
	double zz=-35000.0;
	double z0,z1;
	double TT,T0,T1;
	
	long p;
	
	for (p=0;p<nodes_flex;p++){
		for (i=p,cond=0;i<nodes_thermal-nodes_flex && cond==0;i+=nodes_flex){
			if (zz<xyz_thermal[i][2] && zz>=xyz_thermal[i+nodes_flex][2]){
				p0=i;
				p1=i+nodes_flex;
				cond=1;
				
				z1 = xyz_thermal[p1][2];
				z0 = xyz_thermal[p0][2];
				
				T1=T_vec[p1], 
				T0=T_vec[p0];
				
				TT = T1*(zz-z0)/(z1-z0)+ T0*(z1-zz)/(z1-z0);
				
				Temper35[p]=TT;
				
				cond=1;
			}
		}
	}
	
	
	for (i=0;i<nodes;i++){
		TT =peso[i][0]*Temper35[peso_pos[i][0]];
		TT+=peso[i][1]*Temper35[peso_pos[i][1]];
		TT+=peso[i][2]*Temper35[peso_pos[i][2]];
		
		h_temper35[i]=TT;
	}
	
}
