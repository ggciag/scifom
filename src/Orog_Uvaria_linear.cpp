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

extern double uplift_scale;

//extern double time_ofchangeu;

extern double uplift_scale2;

extern double tempo;

extern double *uplift_map;

void orog()
{

    long i;

    double dh;
    float tmax=5.0E7;		
    float uplift_scale_real;
    uplift_scale_real=((uplift_scale2 - uplift_scale)*tempo/tmax) + uplift_scale;
    for (i=0;i<nodes;i++){
        dh = uplift_scale_real*dt*uplift_map[i]/1.0E6;
        if (h_topo[i]>6000.0){
            dh *= 1.0 - (h_topo[i]-6000.0)/1000.0;
            if (dh<0.0) dh=0;
            printf("entrou aqui");
        }

        h_topo[i]+=dh;
        h_bed[i]+=dh;

        h_q[i]+=-dh*RHOC;
        //if (int(tempo)%100000==0){printf("u=%lf\n", uplift_scale);}
	
	}
}
