//Aloca espaço e le os parametros da tpografia dinamica
#include <stdio.h>
#include <stdlib.h>

extern double **Aloc_matrix_real (long p, long n);
extern double *Aloc_vector_real (long n);

extern double *h_topo_din_cumulat;
extern double **topo_din_map;
extern double *tempo_topo_din;
extern long numero_intervalos;
extern double fac_topo_din;
extern double tempo_inicial;
extern long nodes;
extern long nodes_max_aloca;

void aloca_topo_din()
{
	int i,j;

	h_topo_din_cumulat=Aloc_vector_real(nodes_max_aloca);

	FILE *topografia_dinamica;
	topografia_dinamica = fopen("param_topo_din.txt", "r");
	fscanf(topografia_dinamica,"%ld", &numero_intervalos);
	fscanf(topografia_dinamica,"%lf", &fac_topo_din);
	fscanf(topografia_dinamica,"%lf", &tempo_inicial);
	printf("topografia dinamica:%ld intervalos de tempo \n", numero_intervalos);

	topo_din_map=Aloc_matrix_real(numero_intervalos,nodes);
	tempo_topo_din=Aloc_vector_real(numero_intervalos+1);

	for (i=0;i<=numero_intervalos;i++){
		fscanf(topografia_dinamica,"%lf",&tempo_topo_din[i]);
		printf("topografia tempos (Ma): %lf \n", tempo_topo_din[i]);
	}
	fclose(topografia_dinamica);

	FILE *topografia_dinamica_map;
	topografia_dinamica_map = fopen("topo_din_map.txt", "r");
	printf("(teste) valores duas primeiras linhas: \n");
	for(j=0;j<nodes;j++){
		for (i=0;i<numero_intervalos;i++){
			fscanf(topografia_dinamica_map, "%lf", &topo_din_map[i][j]);
			if (j<2){printf("%lf ", topo_din_map[i][j]);}
		}

	}
    fclose(topografia_dinamica_map);
}

