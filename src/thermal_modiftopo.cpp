/*
 *  thermal_modiftopo.cpp
 *  Thermal_hist5
 *
 *  Created by Victor Sacek on 29/03/10.
 *  Copyright 2010 IAG/USP. All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>

extern long nodes_flex;

extern long layers;

extern double *h_top;
extern double *v_adv2D;

extern double **v_adv;

extern double seg_per_ano;

extern double **xyz_thermal;
extern double **xyz_thermal_fut;

void thermal_modiftopo()
{
	long i,j;
	
	for (j=0;j<layers;j++){
		for (i=0;i<nodes_flex;i++){
			//v_adv[j*nodes_flex+i][2]=100.0/((tempo_fut-tempo_ant)*seg_per_ano);
			//v_adv[j*nodes_flex+i][0]=500.0/1.0E6/seg_per_ano;
			//v_adv[j*nodes_flex+i][2]=500.0/1.0E6/seg_per_ano;
		}
	}
	
	for (i=0;i<nodes_flex;i++){
		h_top[i]=0;
	}
	
	for (i=0;i<nodes_flex;i++){
		xyz_thermal_fut[i][2]=h_top[i];
	}
	
	for (j=1;j<layers;j++){
		for (i=0;i<nodes_flex;i++){
			xyz_thermal_fut[j*nodes_flex+i][2]=xyz_thermal[j*nodes_flex+i][2];
		}
	}
	
}
