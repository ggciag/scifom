/*
 *  lateral.cpp
 *  Thermal_hist3
 *
 *  Created by Victor Sacek on 22/03/10.
 *  Copyright 2010 IAG/USP. All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>

long **Aloc_matrix_long (long p, long n);

extern long layers;

extern long ntri_lat;

extern long **tri_lat;

extern long nodes_thermal;

void lateral(double **xyz, long ntet, long **tetra){
	long t,i;
	double x[4],y[4];
	long verif;
	long conttri=0;
	long aux=0;
	free(tri_lat);
	tri_lat = Aloc_matrix_long(20000, 3);
	
	double xmin,xmax,ymin,ymax;
    
    xmin=xyz[0][0]; xmax=xyz[0][0];
    ymin=xyz[0][1]; ymax=xyz[0][1];    
    
    for (i=1;i<nodes_thermal;i++){
        if (xmin>xyz[i][0]) xmin=xyz[i][0];
        if (xmax<xyz[i][0]) xmax=xyz[i][0];
        if (ymin>xyz[i][1]) ymin=xyz[i][1];
        if (ymax<xyz[i][1]) ymax=xyz[i][1];
    }
	
	
	for (t=0,conttri=0;t<ntet;t++){
		
		x[0]=xyz[tetra[t][0]][0];
		x[1]=xyz[tetra[t][1]][0];
		x[2]=xyz[tetra[t][2]][0];
		x[3]=xyz[tetra[t][3]][0];
		
		y[0]=xyz[tetra[t][0]][1];
		y[1]=xyz[tetra[t][1]][1];
		y[2]=xyz[tetra[t][2]][1];
		y[3]=xyz[tetra[t][3]][1];
		
		for (i=0,verif=0;i<4;i++){
			if (x[i]==xmin) verif++;
		}
		if (verif==3){
			aux=0;
			for (i=0;i<4;i++){
				if (x[i]==xmin){
					tri_lat[conttri][aux]=tetra[t][i];
					aux++;
				}
			}
			conttri++;
		}
		
		for (i=0,verif=0;i<4;i++){
			if (y[i]==ymin) verif++;
		}
		if (verif==3){
			aux=0;
			for (i=0;i<4;i++){
				if (y[i]==ymin){
					tri_lat[conttri][aux]=tetra[t][i];
					aux++;
				}
			}
			conttri++;
		}/**/
		
	}
	
	
	ntri_lat=conttri;
	//return(conttri);
	
	
	
}