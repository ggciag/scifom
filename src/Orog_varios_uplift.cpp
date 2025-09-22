//
//  Orog.cpp
//  Orogenia1.3
//
//  Created by LabTectonofisica on 12/11/12.
//  Copyright (c) 2012 Lab Tectonofísica. All rights reserved.
//

#include <math.h>
#include <stdio.h>

extern long nodes;
extern double *h_topo;
extern double *h_bed;

extern double minx;
extern double maxx;
extern double miny;
extern double maxy;


extern double dt;
extern double **xy;

extern double RHOC;

extern double *h_q;

extern double *tempos_uplift_min;
extern double *tempos_uplift_max;

extern long numero_uplift;
extern double **uplift_map;
extern double *uplift_factor;

extern double tempo;


extern int lith_flag;
extern int nsr;
extern double **h_sr;

void orog()
{

    long i,j;
    double dh;

    for (j=0;j<nodes;j++){
        dh=0;
        for (i=0;i<numero_uplift;i++){
            if (tempo>=tempos_uplift_min[i] && tempo<tempos_uplift_max[i]){
                    dh = dh + uplift_factor[i]*dt*uplift_map[j][i];
            }
        }
        if (h_topo[j]>6000.0){
            dh *= 1.0 - (h_topo[j]-6000.0)/1000.0;
            if (dh<0.0) dh=0;
        }

        h_topo[j]+=dh;
        h_bed[j]+=dh;
        h_q[j]+=-dh*RHOC;

        if (lith_flag==1){
            for (int l=0;l<nsr;l++){
                h_sr[j][l]+=dh;
            }
        }
        
    }
        //if (int(tempo)%100000==0){printf("u=%lf\n", uplift_scale);}

}

