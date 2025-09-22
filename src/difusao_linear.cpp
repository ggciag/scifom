#include <stdio.h>
#include <stdlib.h>
#include <math.h>

extern long **conec;
extern long *pos_conec;
extern long nodes;

extern double *area_vor;
extern double **aresta_vor;
extern double **dist_vor;

extern double *h_topo;

extern double *h_bed;

extern double *Ds;

extern double K_d;
extern double K_m;

extern double dt;

extern double nivel;

void difusao()
{
    long i,j;
    double h_var,K_aux,sed_t;
    double zl = 1.0;

	double vol;

    for (i=0;i<nodes;i++){
        Ds[i]=0;
    }

    for (i=0;i<nodes;i++){
        if (h_topo[i]<nivel){
            if (h_topo[i]>nivel-zl) K_aux=(nivel-h_topo[i])*(K_m-K_d)/zl+K_d;
            else {
                 K_aux = K_m;
                 //printf("K_aux: %f %f\n",K_aux,h_topo[i]);
            }
            sed_t = h_topo[i]-h_bed[i];
            /*if (sed_t<1.0){
                K_aux=K_aux*sed_t/1.0;
            }
            else {
                K_aux=K_aux;
            } */
			for (j=0;j<pos_conec[i];j++){
				h_var = h_topo[i]-h_topo[conec[i][j]];
				if (h_var>0){
					vol = K_aux*aresta_vor[i][j]*(h_var/dist_vor[i][j])*dt;
					if (vol<sed_t*area_vor[i]+Ds[i]){
						Ds[i]-= vol;
						Ds[conec[i][j]]+=vol;
					}
					else {
						vol = sed_t*area_vor[i]+Ds[i];
						Ds[i]-=vol;
						Ds[conec[i][j]]+=vol;
					}

				}
			}

        }
        else{
            K_aux=K_d;
			for (j=0;j<pos_conec[i];j++){
				h_var = h_topo[i]-h_topo[conec[i][j]];
				if (h_var>0){
					Ds[i]-=K_aux*aresta_vor[i][j]*(h_var/dist_vor[i][j])*dt;
					Ds[conec[i][j]]+=K_aux*aresta_vor[i][j]*(h_var/dist_vor[i][j])*dt;
				}
			}
        }/**/
    }


    for (i=0;i<nodes;i++){
        Ds[i]/=area_vor[i];
    }


}
