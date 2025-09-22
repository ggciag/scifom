#include <stdio.h>
#include <stdlib.h>

double **Aloc_matrix_real (long p, long n);
long **Aloc_matrix_long (long p, long n);
double *Aloc_vector_real (long n);
long *Aloc_vector_long (long n);

extern long nodes_max_aloca;

extern double *Df;
extern long *ordem_fluvi;
extern long *direc_fluvi;
extern double *dist_fluvi;
extern double *Qr;
extern double *Qr_prov;
extern double *Qf;
extern long *basins;
extern long **global_basin;
extern long **global_index;
extern double **h_global_basin;
extern long *n_global;
extern long *global_flag;
extern double *h_min;
extern long *pos_min;

extern long *stack_fluvi;

extern double *h_topo_prov;

extern long *lagos;

extern double *lsr;
extern double *depth_lsr;
extern int nsr;
extern double **h_sr;
extern double *lsr_map;

extern long nodes;

extern double *h_topo;

extern int lith_flag;

void aloca_fluvi()
{
     Df = Aloc_vector_real (nodes_max_aloca);
     
     ordem_fluvi = Aloc_vector_long (nodes_max_aloca);
     direc_fluvi = Aloc_vector_long (nodes_max_aloca);
     dist_fluvi = Aloc_vector_real (nodes_max_aloca);
     
     Qr = Aloc_vector_real (nodes_max_aloca);
     Qr_prov = Aloc_vector_real (nodes_max_aloca);	
     Qf = Aloc_vector_real (nodes_max_aloca); 
     
     basins = Aloc_vector_long (nodes_max_aloca);
     global_basin = Aloc_matrix_long (nodes_max_aloca,200);
	 global_index = Aloc_matrix_long (nodes_max_aloca,200);
	 n_global = Aloc_vector_long (nodes_max_aloca);	
	 global_flag = Aloc_vector_long (nodes_max_aloca);	
	 h_global_basin = Aloc_matrix_real (nodes_max_aloca,200);
	
	stack_fluvi = Aloc_vector_long(nodes_max_aloca);
	
	lagos = Aloc_vector_long(nodes_max_aloca);
	
	h_topo_prov = Aloc_vector_real(nodes_max_aloca);
	
     h_min = Aloc_vector_real (nodes_max_aloca);
     pos_min = Aloc_vector_long (nodes_max_aloca);

     FILE *f;
     nsr=0;
     f = fopen("lithification.txt","r");
     if (f!=NULL){
          lith_flag=1;
          fscanf(f,"%d",&nsr);
          
          printf("\n\nLithification on\n\n");
          depth_lsr = Aloc_vector_real(nsr);
          lsr = Aloc_vector_real(nsr);
          lsr_map = Aloc_vector_real(nodes_max_aloca);
          h_sr = Aloc_matrix_real(nodes_max_aloca,nsr);
          for (int j=0;j<nodes;j++){
               for (int i=0;i<nsr;i++){
                    h_sr[j][i]=h_topo[j];
               }
          }
          for (int i=0;i<nsr;i++){
               fscanf(f,"%lf %lf",&lsr[i],&depth_lsr[i]);
               printf("%lf %lf\n",lsr[i],depth_lsr[i]);
          }
          fclose(f);
     }
     else {
          lith_flag=0;
          printf("\n\nLithification off\n\n");
          lsr_map = Aloc_vector_real(nodes_max_aloca);
     }
     
}
