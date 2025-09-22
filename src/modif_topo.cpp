#include <stdio.h>
#include <stdlib.h>

extern double *h_topo;
extern double **xy;

extern double maxx;
extern double minx;


extern double Sed_esq;
extern double Sed_dir;
extern double Sed_in;

extern double *h_bed;
extern double *h_w;
extern double *Ds;
extern double *Df;
extern double *area_vor;
extern long nodes;

extern double *h_foot;

//extern double *h_crust_sup;

extern double *h_q;

extern long *cond_topo_modif;


extern double RHOW;
extern double RHOS;
extern double RHOC;
extern double RHOM;

extern double nivel;


extern double *moho;
extern double *moho_flex;

extern double **peso;
extern long **peso_pos;

extern double *h_crust_sup;

extern double *lsr;
extern double *depth_lsr;
extern int nsr;
extern double **h_sr;

void modif_topo()
{
    long i; 
    double aux_bed;  
	double aux_topo;
	double aux_dh;
	double dh;
	
    
    for (i=0;i<nodes;i++){   
        
		//if (cond_topo_modif[i]=0) h_topo[i]=0.0;
		
		aux_dh=cond_topo_modif[i]*(Ds[i]+Df[i]);
		aux_topo=h_topo[i];
		aux_bed=h_bed[i];
		
		if (aux_dh>0){
			h_q[i]+=-aux_dh*RHOS;//sedimentacao
		}
		else {
			if (aux_topo==aux_bed) {
				h_q[i]+=-aux_dh*RHOC; //erosao embasamento
			}
			else {
				if (aux_topo+aux_dh<aux_bed) {
					dh = aux_topo-aux_bed;
					h_q[i]+=(dh)*RHOS-(aux_dh+dh)*RHOC;  //erosao embasamento + sedimento
				}
				else {
					h_q[i]+=-aux_dh*RHOS; //erosao sedimento
				}

			}

		}
		
		if (aux_topo<=nivel && aux_topo+aux_dh<=nivel){
			h_q[i]+=aux_dh*RHOW;
		}
		else {
			if (aux_topo<=nivel)		h_q[i]+=(nivel-aux_topo)*RHOW;
			if (aux_topo+aux_dh<=nivel) h_q[i]+=-(nivel-aux_topo-aux_dh)*RHOW;			
		}

			

        h_topo[i]+=aux_dh;

		for (int lit=0; lit<nsr; lit++){
			if (h_topo[i]-depth_lsr[lit] > h_sr[i][lit]) h_sr[i][lit] = h_topo[i]-depth_lsr[lit];
			if (h_sr[i][lit]>h_topo[i]) h_sr[i][lit]=h_topo[i];
		}

		
        aux_bed=h_bed[i];
        if (aux_bed>h_topo[i]){
            h_bed[i]=h_topo[i];
            //if (h_foot[i]>h_bed[i]) h_foot[i]=h_bed[i];
        }
	}
	
	
	
	Sed_in=0;
	for (i=0;i<nodes;i++){           
		Sed_in+=(h_topo[i]-h_bed[i])*area_vor[i]/1.0E9;
		if (cond_topo_modif[i]==0){
			if (xy[i][0]<(minx+maxx)/2){
				if (Ds[i]+Df[i]>0){
					Sed_esq+=(Ds[i]+Df[i])*area_vor[i]/1.0E9;
				}
			}
			else {
				if (Ds[i]+Df[i]>0){
					Sed_dir+=(Ds[i]+Df[i])*area_vor[i]/1.0E9;
				}
			}
		}
		
	}
}
