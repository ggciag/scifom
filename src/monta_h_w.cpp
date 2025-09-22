extern double **peso;
extern long **peso_pos;

extern long nodes;
extern long nodes_flex;

extern double *h_topo;

extern double *h_bed;
extern double *h_flex_cumulat;
extern double *h_w;

extern double *h_foot;

extern double *moho;
extern double *h_crust_sup;

extern double *wflex_aux;

extern double **falha_plot; //posi��o da falha em y = 0
extern long **falha_plot_pos;
extern long num_falha;


extern double *moho_flex;
extern double **xy_flex;

extern double *h_q;
extern double *h_q_a;
extern double depth;
extern double nivel;

extern double RHOM;
extern double RHOC;
extern double RHOS;
extern double RHOW;

extern double TeConstante;
extern double TeConstante2;

extern double minx;
extern double maxx;


extern int lith_flag;
extern int nsr;
extern double **h_sr;

void monta_h_w()
{
    long i,j;
	            
	for (i=0;i<nodes;i++){
		h_w[i] =peso[i][0]*wflex_aux[peso_pos[i][0]];
		h_w[i]+=peso[i][1]*wflex_aux[peso_pos[i][1]];
		h_w[i]+=peso[i][2]*wflex_aux[peso_pos[i][2]];
	}
	
	for (i=0;i<nodes_flex;i++){
		if (xy_flex[i][0]>=minx && xy_flex[i][0]<=maxx)
			moho_flex[i]+=wflex_aux[i];
	}
	
	double aux_dh;
	double aux_topo;
	
	
    for (i=0;i<nodes;i++){
		
		aux_dh=h_w[i];
		aux_topo=h_topo[i];
		
		if (aux_topo<=nivel && aux_topo+aux_dh<=nivel){
			h_q_a[i]=aux_dh*RHOW;
		}
		else {
			if (aux_topo<=nivel)		h_q_a[i]=(nivel-aux_topo)*RHOW;
			if (aux_topo+aux_dh<=nivel) h_q_a[i]=-(nivel-aux_topo-aux_dh)*RHOW;			
		}
		if (aux_topo>nivel && aux_topo+aux_dh>nivel){
			h_q_a[i]=0;
		}/**/
		
		
        h_topo[i]+=h_w[i];
        h_bed[i]+=h_w[i];
        h_foot[i]+=h_w[i];
        moho[i]+=h_w[i];//!!!! Ver como interpolar com moho_flex !!!!
		h_crust_sup[i]+=h_w[i];
        
        h_flex_cumulat[i]+=h_w[i];

    }

	if (lith_flag==1){
		for (i=0;i<nodes;i++){
			for (int l=0;l<nsr;l++){
				h_sr[i][l]+=h_w[i];
			}
		}
	}
    
    
    long aux;
    
    for (i=0;i<num_falha;i++){
        aux = falha_plot_pos[i][0];
        for (j=1; j<=aux; j++){
            falha_plot[i][j]+=h_w[falha_plot_pos[i][j]];
            if (falha_plot[i][j]>h_bed[falha_plot_pos[i][j]])
                falha_plot[i][j]=h_bed[falha_plot_pos[i][j]];
        }
    }
    
}
