
//#include <stdio.h>
#include <stdlib.h>

double *Aloc_vector_real (long n);

extern long nodes;
extern long nodes_max_aloca;

extern double *Ds;


void aloca_difusao()
{
     Ds = Aloc_vector_real (nodes_max_aloca);     
}
     
