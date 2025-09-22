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

extern double time_ofchangeu;

extern double uplift_scale2;

extern double tempo;

extern double tempo_max;

extern double *uplift_map;

void orog()
{

    long i;

    double dh;
    float uplift_scale_real;
    for (i=0;i<nodes;i++){

        if (time_ofchangeu != 0){
            if (tempo<=time_ofchangeu){
                dh = uplift_scale*dt*uplift_map[i]/1.0E6;

                h_topo[i]+=dh;
                h_bed[i]+=dh;

                h_q[i]+=-dh*RHOC;
                //if (int(tempo)%100000==0){printf("u=%lf\n", uplift_scale);}
            }else{
                dh = uplift_scale2*dt*uplift_map[i]/1.0E6;

                h_topo[i]+=dh;
                h_bed[i]+=dh;

                h_q[i]+=-dh*RHOC;
                //if (int(tempo)%100000==0){printf("u=%lf\n", uplift_scale2);}
            }
        }else{
            uplift_scale_real=((uplift_scale2 - uplift_scale)*tempo/tempo_max) + uplift_scale;
            dh = uplift_scale_real*dt*uplift_map[i]/1.0E6;

            h_topo[i]+=dh;
            h_bed[i]+=dh;

            h_q[i]+=-dh*RHOC;
        }
    }


}
