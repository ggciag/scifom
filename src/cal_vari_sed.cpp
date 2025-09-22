#include <stdio.h>
#include <stdlib.h>

void IO1 (float *maxy, float *miny, int *nx, int *ny, float *Termic, float *Te, float *Vr, float *Kf, float *Km, float *ls, float *lb, float *U){
	FILE *entrada;
	entrada = fopen ("param_OrogSedFlex_1.1.txt", "r");
	fscanf (entrada,"%f\n%f\n%d\n%d%f\n%f\n%f\n%f\n%f\n%f\n%f\n%f", maxy, miny, ny, nx, Termic, Te, Vr, Kf, Km, ls, lb, U);
	fclose (entrada);
	return;
}

int IOpontos (int **labels, int nx, int ny, double x1, double y1, double x2, double y2){
	double x, y;
				
	FILE *entrada;
	entrada = fopen ("pontos.txt", "r");
	
	int cont=0, dotsinside=0;
	while (cont < nx*ny){
		fscanf (entrada,"%lf %lf\n", &x, &y);
		if (x>x1 && x<x2 && y>y1 && y<y2){
			labels[0][cont]=1;
			dotsinside++;
		}
		cont++;
	}
	fclose (entrada);
	return(dotsinside);
}

double IOgeologia (double tempo, int *labels, int nx, int ny){
	int cont=0, teste=0;
	char nome[30];
	double htopo, hbed, diftopbed;
	float aux1, aux2, aux3, aux4, aux5, aux6, aux7, aux8;
		
	FILE *entrada;

	sprintf(nome, "m_Topo_%.3lf.txt", tempo);
	printf("%s \n", nome);
	entrada = fopen(nome, "r");
	
	diftopbed=0;
	while (cont < nx*ny){
		fscanf (entrada,"%lf %lf %f %f %f %f %f %f %f %f\n", &htopo, &hbed, &aux1, &aux2, &aux3, &aux4, &aux5, &aux6, &aux7, &aux8);
		if(labels[cont]==1){
			diftopbed=diftopbed + (htopo-hbed);
		}
		cont++;
	}
	fclose (entrada);
return (diftopbed);
}

int main(int argc, char **argv){

	int nx, ny, cont, pontosmarcados;
	double x1, y1, x2, y2, tempo, dif, areamedia, volume;;
	float  maxy, miny, Termic, Te, Vr, Kf, Km, ls, lb, U;
	int *labels;
	/*
	if (argc != 5){
		printf("Deve passar os 4 argumentos: vertice inferior esquerdo (x1, y1) e vertice superior direito (x2, y2)");
		return 0;
	}else{
	x1=atof(argv[1]);
	y1=atof(argv[2]);
	x2=atof(argv[3]);
	y2=atof(argv[4]);
	}
	*/
	//printf("x1=%lf y1=%lf x2=%lf y2=%lf \n", x1, y1, x2, y2);


	
	

			
	IO1 (&maxy, &miny, &nx, &ny, &Termic, &Te, &Vr, &Kf, &Km, &ls, &lb, &U);	

	//Abre e salva o arquivo de saída
	char nome[100];
	sprintf(nome, "../../t_sed_Te%.1f_U%.1f.txt", Te, U);
	FILE *saida;
	saida=fopen(nome, "w");
	
	//imprime cabeçalho
	fprintf(saida, "%f\n%f\n%d\n%d\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n", maxy, miny, ny, nx, Termic, Te, Vr, Kf, Km, ls, lb, U);
	
	
	//Area que corresponde ao leque do amazonas
	x1=3500000;
	y1=1600000;
	x2=4900000;
	y2=2800000;
	
	labels=(int*)calloc(nx*ny,sizeof(int));
	pontosmarcados=IOpontos (&labels, nx, ny, x1, y1, x2, y2);
	
	areamedia=(x2-x1)*(y2-y1)/(pontosmarcados);	
	
	//Define o tempo inicial em M.a. após o início da simulação 
	tempo=1;
	
	// Tempo Máximo
	while (tempo <= 80){
		dif=IOgeologia (tempo, labels, nx, ny);
		volume=dif*areamedia;
		fprintf(saida, "%lf %lf\n", tempo, volume);
		// passo de tempo
		tempo = tempo + 0.2;
	}


	free(labels);
		
	// dimensões bacia do solimões aproximada
	
	x1=1500000;
	x2=2600000;
	y1=1000000;
	y2=2200000;
	
	labels=(int*)calloc(nx*ny,sizeof(int));
	pontosmarcados=IOpontos (&labels, nx, ny, x1, y1, x2, y2);
	
	areamedia=(x2-x1)*(y2-y1)/(pontosmarcados);	
	
	fprintf(saida, "#bacia do solimões\n");
	//Define o tempo inicial em M.a. após o início da simulação 
	tempo=1;

	// Tempo Máximo
	while (tempo <= 80){
		dif=IOgeologia (tempo, labels, nx, ny);
		volume=dif*areamedia;
		fprintf(saida, "%lf %lf\n", tempo, volume);
		// passo de tempo
		tempo = tempo + 0.2;
	}


	free(labels);
	
	// bacias de ante-país
	
	x1=800000;
	x2=1450000;
	y1=200000;
	y2=2900000;	
	
	labels=(int*)calloc(nx*ny,sizeof(int));
	pontosmarcados=IOpontos (&labels, nx, ny, x1, y1, x2, y2);
	
	areamedia=(x2-x1)*(y2-y1)/(pontosmarcados);	
	
	fprintf(saida, "#bacias de ante-país\n");
	//Define o tempo inicial em M.a. após o início da simulação 
	tempo=1;

	// Tempo Máximo
	while (tempo <= 80){
		dif=IOgeologia (tempo, labels, nx, ny);
		volume=dif*areamedia;
		fprintf(saida, "%lf %lf\n", tempo, volume);
		// passo de tempo
		tempo = tempo + 0.2;
	}


	free(labels);

	//	
	fclose(saida);	
	return (0);
}
