//Aloca espaço e le os parametros relacionados ao processo de underplate
#include <stdio.h>
#include <stdlib.h>

extern double **Aloc_matrix_real (long p, long n);
extern double *Aloc_vector_real (long n);

extern double *tempos_underp_min;
extern double *tempos_underp_max;

extern long numero_underp;
extern double **underp_map;
extern double *underp_factor;

extern long nodes;
extern long nodes_max_aloca;

void aloca_underp()
{
	int i,j,aux;


	FILE *parametros_underp;
	parametros_underp = fopen("param_underp.txt", "r");


	fscanf(parametros_underp,"%ld\n", &numero_underp);

	tempos_underp_min=Aloc_vector_real(numero_underp);
    tempos_underp_max=Aloc_vector_real(numero_underp);
    underp_factor=Aloc_vector_real(numero_underp);

	for (i=0;i<numero_underp;i++){
        fscanf(parametros_underp,"%lf %lf %lf\n", &underp_factor[i],&tempos_underp_min[i], &tempos_underp_max[i]);
	}
    fclose(parametros_underp);

	printf("numero de fases de underps:%ld \n", numero_underp);
    printf("nodes_max_aloca:%ld \n", nodes_max_aloca);

    underp_map=Aloc_matrix_real(nodes_max_aloca, numero_underp);

	printf("alocado espaco para matriz underp_map \n");

	FILE *arquivo_underp_map;
	arquivo_underp_map = fopen("underp_map.txt", "r");

    printf("lendo arquivo underp_map.txt\n");

    for (j=0;j<nodes;j++){
        for (i=0;i<numero_underp;i++){
            fscanf(arquivo_underp_map,"%lf ",&underp_map[j][i]);
        }
        fscanf(arquivo_underp_map,"\n");
        aux=j+1;
        //printf("linha %d = %f, %f \n", aux, underp_map[j][0], underp_map[j][1]);

    }
    fclose(arquivo_underp_map);
    printf("arquivo underp_map.txt lido\n");
    printf("underpmap linha 2390 = %f, %f \n", underp_map[2389][0], underp_map[2389][1]);
    printf("underpmap linha 1004 = %f, %f \n", underp_map[1004][0], underp_map[1004][1]);
    //printf("%lf %lf \n %lf %lf \n",underp_map[18][0],underp_map[18][1],underp_map[19][0],underp_map[19][1]);
}

