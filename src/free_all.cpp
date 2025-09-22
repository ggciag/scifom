#include <stdio.h>
#include <stdlib.h>


extern long **Tri;
extern long **conec;
extern long *pos_conec;
extern double **xy;

extern double *area_vor;
extern double **aresta_vor;
extern double **dist_vor;

extern double *h_topo;
extern double *h_bed;
extern double *water;
extern double *h_q;
extern double *h_w;

extern double *Ds;
extern double *Df;

extern long *ordem_fluvi;
extern long *direc_fluvi;
extern double *dist_fluvi;
extern double *Qr;
extern double *Qf;
//extern FILE *f_h_topo;
//////////////////////////////////////////
//Flexura

extern long **Tri_flex;
extern double **xy_flex;
extern long **cond_c;
extern double **peso;
extern long **peso_pos;
extern double **Kflex;
extern long **Kconec;
extern long *Kposconec;
extern double **Ke;
extern double *Te;
extern double *props;
extern double *Kdiag;
//extern double *qflex;
//extern double *qflex_a;
extern double *bflex;
extern double *wflex;
extern double *wflex_aux;
extern double *wflex_cumula;

void free_all()
{
    free( Tri);
    free( conec);
    free( pos_conec);
    free( xy);
    
    free( area_vor);
    free( aresta_vor);
    free( dist_vor);
    
    free( h_topo);
    free( h_bed);
    free( water);
    free( h_q);
    free( h_w);
    
    free( Ds);
    free( Df);
    
    free( ordem_fluvi);
    free( direc_fluvi);
    free( dist_fluvi);
    
    free( Qr);
    free( Qf);
	
    
    //free(f_h_topo);
    
    
    
    //////////////////////////////////////////
    //Flexura
    
    free( Tri_flex);
    free( xy_flex);
    
    free( cond_c);
    
    
    free( peso);
    free( peso_pos);
    
    free( Kflex);
    free( Kconec);
    free( Kposconec);
    
    free( Ke);
    free( Te);
    free( props);
    free( Kdiag);
    //free( qflex);
    //free( qflex_a);
    free( bflex);
    free( wflex);
    free( wflex_aux);
    free( wflex_cumula);
    
}
