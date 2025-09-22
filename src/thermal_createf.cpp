/*
 *  thermal_createf.cpp
 *  thermal1
 *
 *  Created by Victor Sacek on 15/03/10.
 *  Copyright 2010 IAG/USP. All rights reserved.
 *
 */

extern double *fthermal;
extern long nodes_thermal;

extern double **Mthermal;
extern double *Mthermal_diag;
extern double **Kthermal1;
extern double *Kthermal1_diag;

extern long **Kthermal_conec;
extern long *Kthermal_posconec;

extern double dt_calor_sec;

extern double comp_alpha_thermal;

extern double *T_vec;



void thermal_createf()
{
	long i,j,pos_aux;
	for (i=0;i<nodes_thermal;i++){
		fthermal[i]=(Mthermal_diag[i]-comp_alpha_thermal*dt_calor_sec*Kthermal1_diag[i])*T_vec[i];
		pos_aux = Kthermal_posconec[i];
		for (j=0;j<pos_aux;j++){
			fthermal[i]+=(Mthermal[i][j]-comp_alpha_thermal*dt_calor_sec*Kthermal1[i][j])*T_vec[Kthermal_conec[i][j]];
		}
	}
	
}