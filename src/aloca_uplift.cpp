//Aloca espa�o e le os parametros da tpografia dinamica
#include <stdio.h>
#include <stdlib.h>

extern double **Aloc_matrix_real (long p, long n);
extern double *Aloc_vector_real (long n);

extern double *tempos_uplift_min;
extern double *tempos_uplift_max;

extern long numero_uplift;
extern double **uplift_map;
extern double *uplift_factor;

extern long nodes;
extern long nodes_max_aloca;

void aloca_uplift()
{
	int i,j;


	FILE *parametros_uplift;
	parametros_uplift = fopen("param_uplift.txt", "r");


	fscanf(parametros_uplift,"%ld\n", &numero_uplift);

	tempos_uplift_min=Aloc_vector_real(numero_uplift);
    tempos_uplift_max=Aloc_vector_real(numero_uplift);
    uplift_factor=Aloc_vector_real(numero_uplift);

	for (i=0;i<numero_uplift;i++){
        fscanf(parametros_uplift,"%lf %lf %lf\n", &uplift_factor[i],&tempos_uplift_min[i], &tempos_uplift_max[i]);
	}
    fclose(parametros_uplift);

	printf("numero de fases de uplifts:%ld \n", numero_uplift);
    printf("nodes_max_aloca:%ld \n", nodes_max_aloca);

    uplift_map=Aloc_matrix_real(nodes_max_aloca, numero_uplift);

	printf("alocado espaco para matriz uplift_map");

	FILE *arquivo_uplift_map;
	arquivo_uplift_map = fopen("uplift_map.txt", "r");

    printf("lendo arquivo uplift_map.txt");

    for (j=0;j<nodes;j++){
        for (i=0;i<numero_uplift;i++){
            fscanf(arquivo_uplift_map,"%lf ",&uplift_map[j][i]);
        }
        fscanf(arquivo_uplift_map,"\n");
    }
    fclose(arquivo_uplift_map);
    printf("arquivo uplift_map.txt lido");
    //printf("%lf %lf \n %lf %lf \n",uplift_map[18][0],uplift_map[18][1],uplift_map[19][0],uplift_map[19][1]);
}

