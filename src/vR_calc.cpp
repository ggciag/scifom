/*
 *  vR_calc.cpp
 *  Orogenia1.17_Te_variavel_3_m0.5
 *
 *  Created by Victor Sacek on 08/04/13.
 *  Copyright 2013 IAG/USP. All rights reserved.
 *
 */


//vR constante!!!

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

extern double **uplift_map;

extern long *vR_ordem;
extern double *vR_aux;

extern double *vR_map;
extern double *vR_flow;

extern double *h_topo;

extern long **conec;
extern long *pos_conec;

extern double **aresta_vor;
extern double **dist_vor;
extern double *area_vor;

extern long nodes;
extern long nivel;

extern double vR;
extern double time_ofchangevR;
extern double vR2;

extern double vRandes;
extern double time_ofchangevRandes;
extern double vR2andes;


extern double tempo;
extern double tempo_max;

extern double **xy;

extern double minx;
extern double maxx;
extern double miny;
extern double maxy;

extern double windx;
extern double windy;

extern double *Lf_vec;

void vR_calc()
{


	long i;

	for (i=0;i<nodes;i++){
		vR_map[i]=0.0;
		vR_flow[i]=0.0;
	}
	
    double vR_efetivo;
    double vR_efetivoandes;
	for (i=0;i<nodes;i++){
		if (uplift_map[i][0]==0){
            		if (time_ofchangevR!=0){
                		if (tempo<=time_ofchangevR){
                    			vR_map[i]=vR*area_vor[i]; //// vR_map = 1 * Area // vR constante
                		}else{
                    			vR_map[i]=vR2*area_vor[i]; //// vR_map = 1 * Area // vR constante
                		}
            		}else{
                		vR_efetivo = ((vR2 - vR)*tempo/tempo_max) + vR;
                		vR_map[i]=vR_efetivo*area_vor[i]; //// vR_map = 1 * Area // vR constante
		    	}
		}else{
            		if (time_ofchangevRandes!=0){
                		if (tempo<=time_ofchangevRandes){
                    			vR_map[i]=vRandes*area_vor[i]; //// vR_map = 1 * Area // vR constante
                		}else{
                    			vR_map[i]=vR2andes*area_vor[i]; //// vR_map = 1 * Area // vR constante
                		}
            		}else{
                		vR_efetivoandes = ((vR2andes - vRandes)*tempo/tempo_max) + vRandes;
                		vR_map[i]=vR_efetivoandes*area_vor[i];
            		}
		}
	}
}

void vR_calc_orography()
{


	double vapour = 2.0E6; //m^2/yr
	double l_oro = 100.0E3;
	double h_oro = 1000.0;

	long i,ii,j,jj;

	double dist_aux;



	for (i=0;i<nodes;i++){
		if (h_topo[i]<nivel){
			vR_flow[i]=vapour;
		}
		else {
			if (xy[i][0]==maxx && windx<0.0) vR_flow[i]=vapour;
			if (xy[i][0]==minx && windx>0.0) vR_flow[i]=vapour;
			if (xy[i][1]==maxy && windy<0.0) vR_flow[i]=vapour;
			if (xy[i][1]==miny && windy>0.0) vR_flow[i]=vapour;
		}

	}

	double xxn;

	xxn = (minx+0.32*(maxx-minx));



	for (ii=nodes-1;ii>=0;ii--){
		i=vR_ordem[ii];

		vR_map[i]=3.0*area_vor[i];//!!!

		if (h_topo[i]>nivel){
			jj=i;


			for (j=0;j<pos_conec[i];j++){
				if (vR_aux[conec[i][j]]>vR_aux[jj]){
					if (aresta_vor[i][j]>0){
						jj = conec[i][j];
						dist_aux = dist_vor[i][j];
					}
				}
			}

			if (vR_flow[i]==0.0 && jj!=i){
				//if (xy[i][0]<(maxx-minx)/3.0){
					vR_flow[i]=vR_flow[jj]-vR_map[jj]*dist_aux/area_vor[jj];
					if (vR_flow[i]<0.0) vR_flow[i]=0.0;
				//}
				//else vR_flow[i]=vapour;

			}


			//if (h_topo[i]>500 && xy[i][0]<xxn)
			if (xy[i][0]<xxn){
				if (h_topo[i]>500 && uplift_map[i][0]>0)
					vR_map[i]=vR_flow[i]*(h_topo[i]-450.)*area_vor[i]/(l_oro*h_oro);
				else
					vR_map[i]=vR_flow[i]*(50.)*area_vor[i]/(l_oro*h_oro);
			}
			else{
				vR_map[i]=1.0*area_vor[i];
				vR_flow[i]=vapour;
			}

			vR_map[i]=1.0*area_vor[i];//!!!

			//printf("%ld %ld %ld vR_map: %lg | vR_flow: %lg\n",ii,i,jj,vR_map[i]/area_vor[i],vR_flow[i]);
		}
	}
	

}


extern double **vR_maps;
extern double *h_vR_external;
extern int nvR_maps;


void calc_vR_external()
{
	long i,t;
	double h_mean=0.0;
	long cont_h_mean=0;
	double fac;
	for (i=0;i<nodes;i++){
		if (uplift_map[i][0]>600.0 && xy[i][1]<500.0E3){
			//if (h_max<h_topo[i]) h_max=h_topo[i];
			h_mean+=h_topo[i];
			cont_h_mean+=1;
		}
	}
	h_mean=h_mean/cont_h_mean;
	printf("h_mean = %.3f\n",h_mean);
	if (h_mean<h_vR_external[0]){
		for (i=0;i<nodes;i++){
			vR_map[i] = vR_maps[0][i];
		}
	}
	for (t=1;t<nvR_maps;t++){
		if (h_mean>=h_vR_external[t-1] && h_mean<h_vR_external[t]){
			fac = (h_mean-h_vR_external[t-1])/(h_vR_external[t]-h_vR_external[t-1]);
			for (i=0;i<nodes;i++){
				vR_map[i] = (1.0-fac)*vR_maps[t-1][i] + fac*vR_maps[t][i];
			}
		}
	}
	t = nvR_maps-1;
	if (h_mean>=h_vR_external[t]){
		for (i=0;i<nodes;i++){
			vR_map[i] = vR_maps[t][i];
		}
	}

	for (i=0;i<nodes;i++){
		vR_map[i] = vR_map[i]*area_vor[i];;
	}
	

	//for (i=0;i<nodes;i++){
	//	printf("%f ",vR_map[i]);
	//}
}