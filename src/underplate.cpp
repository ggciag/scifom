//
//  UnderPlate.cpp
//
//
//  Created by LabTectonofisica on 02/11/20.
//  Copyright (c) 2020 Lab Tectonofísica. All rights reserved.
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
extern double RHOM;


extern double *h_q;

extern double *tempos_underp_min;
extern double *tempos_underp_max;

extern long numero_underp;
extern double **underp_map;
extern double *underp_factor;

extern double tempo;

void underplate()
{

    long i,j;
    double dh;

    for (j=0;j<nodes;j++){
        dh=0;
        for (i=0;i<numero_underp;i++){
            if (tempo>=tempos_underp_min[i] && tempo<tempos_underp_max[i]){
                    dh = dh + underp_factor[i]*dt*underp_map[j][i]/1.0E6;
            }
        }
        h_q[j]+=dh*(RHOM-RHOC);
        //if (h_q[j]!=0){
            //printf("Carga underplate: %f\n", h_q[j]);
        //}
    }
        //if (int(tempo)%100000==0){printf("u=%lf\n", uplift_scale);}

}

