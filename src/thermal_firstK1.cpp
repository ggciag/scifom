/*
 *  thermal_firstK1.cpp
 *  thermal5
 *
 *  Created by Victor Sacek on 19/03/10.
 *  Copyright 2010 IAG/USP. All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>

extern double **Kthermal1;
extern double *Kthermal1_diag;

extern long **Kthermal_conec;
extern long *Kthermal_posconec;

extern double **TKe;
extern double **TMe;

extern long tetra_thermal;
extern long **Tet_thermal;
extern long nodes_thermal;
extern double **xyz_thermal;

extern double dt_calor_sec;

extern double seg_per_ano;

extern double kappa;

extern double alpha_thermal;

extern double soma_volume;


void montaKeThermal(double **Ke, double **Me, long **Tet_thermal, double **xyz_thermal, long t);

void thermal_firstK1()
{
	long t,i,ii,jj,iaux,pos_aux;
	
	
	for (i=0;i<nodes_thermal;i++){
		Kthermal1_diag[i]=0.0;
		pos_aux = Kthermal_posconec[i];
		for (iaux=0;iaux<pos_aux;iaux++){
			Kthermal1[i][iaux]=0.0;
		}			
	}
	
	soma_volume=0.0;
	
	for (t=0; t<tetra_thermal; t++){
		montaKeThermal(TKe,TMe,Tet_thermal,xyz_thermal,t);		
		for (ii=0;ii<4;ii++){
			pos_aux = Tet_thermal[t][ii];
			for (jj=0;jj<4;jj++){				
				if (ii==jj){
					Kthermal1_diag[pos_aux]+=TKe[ii][jj];
				}   
				else {
					for (iaux=0;iaux<Kthermal_posconec[pos_aux];iaux++){
						if (Kthermal_conec[pos_aux][iaux]==Tet_thermal[t][jj]){
							Kthermal1[pos_aux][iaux]+=TKe[ii][jj];
							iaux=4000;//numero bem grande para sair do for
						}
					}
				}                                  
			}
		}
	}	
	
	
}
