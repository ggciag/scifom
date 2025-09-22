//
//  topo_din.cpp
//  incorpora a topografia dinamica lida de um arquivo (topo_din_map.txt)
//	para cada intervalo de tempo a topografia dinamica pode variar,
//os intervals de tempo estao no arquivo param_topo_din.txt (a primeira linha � o numero de intvalos de tempo
//  Created by LabTectonofisica on 30/04/18.
//  Copyright (c) 2018 Lab Tectonofísica. All rights reserved.
//

#include <math.h>
#include <stdio.h>
extern long nodes;
extern double *h_topo;
extern double *h_bed;
extern double *h_crust_sup;
extern double *moho;
extern double *h_flex_cumulat;
extern double *h_topo_din_cumulat;
extern double fac_topo_din;

extern double dt;

extern double RHOC;

//extern double *h_q;

extern double tempo;

extern double **topo_din_map;
extern double *tempo_topo_din;
extern long numero_intervalos;
extern double tempo_inicial;


extern int lith_flag;
extern int nsr;
extern double **h_sr;

void topo_din()
{

    long i = 0;
	//long verif = 0;
    double dh;
    double tempo_aux;
    double tempo_aux2;

	long n_tempo;

    for (i=0;i<numero_intervalos;i++){
        tempo_aux=tempo_inicial-tempo_topo_din[i];
    	tempo_aux2=tempo_inicial-tempo_topo_din[i+1];
		if (tempo/1000000>=tempo_aux && tempo/1000000<tempo_aux2){
			n_tempo=i;
        }
    }

    for (i=0;i<nodes;i++){
        dh = fac_topo_din*(topo_din_map[n_tempo][i])*(dt/1000000);
        h_topo[i]+=dh;
        h_bed[i]+=dh;
		h_crust_sup[i]+=dh;
		moho[i]+=dh;
        h_flex_cumulat[i]+=dh;
        h_topo_din_cumulat[i]+=dh;
		//h_q[i]+=-dh*RHOC; deve-se colocar a carga h[q]?
        //if (int(tempo)%100000==0){printf("u=%lf\n", uplift_scale);}
        if (lith_flag==1){
            for (int l=0;l<nsr;l++){
                h_sr[i][l]+=dh;
            }
        }
	}
}


extern double *tempo_sea;
extern double *h_sea;
extern long n_sea_levels;

void sea_level_change()
{
    long t=0;
    double dh;
    double sea1=0.0,sea2=0.0,aux;

    for (t=0;t<n_sea_levels-1;t++){
        if (tempo-dt>tempo_sea[t] && tempo-dt<=tempo_sea[t+1]){
            aux = (tempo-dt-tempo_sea[t])/(tempo_sea[t+1]-tempo_sea[t]);
            sea1 = h_sea[t]*(1.0-aux) + h_sea[t+1]*aux;
        }
        if (tempo>tempo_sea[t] && tempo<=tempo_sea[t+1]){
            aux = (tempo-tempo_sea[t])/(tempo_sea[t+1]-tempo_sea[t]);
            sea2 = h_sea[t]*(1.0-aux) + h_sea[t+1]*aux;
        }
    }
    dh = sea1-sea2;
    printf("sea1=%lf sea2=%lf dh=%lf\n",sea1,sea2,dh);

    for (long i=0;i<nodes;i++){
        h_topo[i]+=dh;
        h_bed[i]+=dh;
		h_crust_sup[i]+=dh;
		moho[i]+=dh;
        h_flex_cumulat[i]+=dh;
        h_topo_din_cumulat[i]+=dh;
		//h_q[i]+=-dh*RHOC; deve-se colocar a carga h[q]?
        //if (int(tempo)%100000==0){printf("u=%lf\n", uplift_scale);}
        if (lith_flag==1){
            for (int l=0;l<nsr;l++){
                h_sr[i][l]+=dh;
            }
        }
	}

}
