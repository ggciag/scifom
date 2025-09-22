/*
 *  aloca_falha_3D.cpp
 *  Flex_Erod_2_6_2
 *
 *  Created by Victor Sacek on 14/06/12.
 *  Copyright 2012 IAG/USP. All rights reserved.
 *
 */

double **Aloc_matrix_real (long p, long n);

#include <stdio.h>
#include <stdlib.h>

extern long n_falhas;
extern double **falhas_param;

void aloca_falha_3D()
{
	FILE *f_falha;
	
	f_falha = fopen("falhas.txt","r");
	
	fscanf(f_falha, "%ld",&n_falhas);
	
	long i,j;
	
	if (n_falhas>0){
		
		falhas_param = Aloc_matrix_real(n_falhas, 9);
		
		for (i=0; i<n_falhas; i++){
			for (j=0;j<9;j++){
				fscanf(f_falha, "%lf",&falhas_param[i][j]);
			}
		}
	}
	
	fclose(f_falha);

}