//
//
// Sem DIFUSAO no CONTINENTE!!!
//
// Fluvi_min:  Area^m com m=1.0
//
//
////////////////////////
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "header.h"
#include <sys/time.h>

/*//Função que lê o arquivo de entrada e salva dados dos sismos lidos. Retorna o numero de sismos lidos.
void IO (*nx, *ny){

	char linha[20];
	int cont=0;
	FILE *entrada;

	entrada = fopen ("param_OrogSedFlex_1.1.txt", "r");

	fscanf (entrada, "*\n*\n%ld\n%ld", ny, nx);
	fclose (entrada);
	return ();
}*/

void malha();
void malha_regular(double minx,double maxx,double miny,double maxy,long Nx, long Ny);

void aloca_difusao();
void aloca_topo();

void aloca_topo_din();
void topo_din();

void aloca_fluvi();
void aloca_uplift();
void aloca_underp();

void underplate();

void print_topo(double tempo,long muda_ponto);


void para_V3D();

void difusao();
void modif_topo();

void fluvi();
void aloca_fluvi();
void fluvi_min();

void flexuraK();
void flexura_load();
void flexura_solv(double **Kflex, double *Kdiag, double *bflex, double *wflex);

void aloca_calor(double minx,double maxx,double miny,double maxy,double depth, long Nx, long Ny, long Nz);
void calor(double dt_calor);
void perturb_thermal(double tempo);

void monta_h_w();


void modif_moho();

void free_all();

void calcula_IsoT();

long gera_malha(long n, double **xy, long **Tri);

void desloca(long n, double **xy,double pos_falha,double despl,long tri,long **Tri);

void re_voronoi(long n, double **xy, long tri, long **Tri);

long gera_malha_parc(long n, long n_old, double **xy, long **Tri,long tri);

double var_fault(double v);

void aloca_falhas();
void monta_falha();

void thermal_modiftopo();
void thermal_monta_K();
void thermal_createf();
void thermal_solv(double **Kf, double *Kd, double *b_vec, double *x_vec);
void thermal_renew();
void thermal_firstK1();
void thermal_plot_all();

void stream();
void desloca_tudo(long n, double **xy, long tri, long **Tri);

void fault3D();

double **Aloc_matrix_real (long p, long n);
long **Aloc_matrix_long (long p, long n);
double *Aloc_vector_real (long n);
long *Aloc_vector_long (long n);

void printf_profile();

void aloca_falha_3D();

void orog();

void vR_sort();
void vR_calc();
void vR_calc_orography();
void calc_vR_external();
void read_vR_external();
void read_sea_level();

void sea_level_change();

int main()
{

    long i,j;
    double soma_h;

    nivel=0.0;

	FILE *entra_var;
	entra_var = fopen("param_OrogSedFlex_1.1.txt", "r");

	fscanf(entra_var,"%lf",&maxy);
	fscanf(entra_var,"%lf",&miny);

	fscanf(entra_var,"%ld",&n_lat);
	fscanf(entra_var,"%ld",&n_latx);


	fscanf(entra_var, "%lf",&axis_stream);
	fscanf(entra_var, "%lf",&Terigida);
	fscanf(entra_var, "%lf",&Teoffshore);
	//TeConstante2=TeConstante;
	//fscanf(entra_var, "%lf",&Telitho);
	//Telitho = TeConstante;

	fscanf(entra_var,"%lf",&vR);
	fscanf(entra_var,"%lf",&time_ofchangevR);
	fscanf(entra_var,"%lf",&vR2);

	fscanf(entra_var,"%lf",&vRandes);
	fscanf(entra_var,"%lf",&time_ofchangevRandes);
	fscanf(entra_var,"%lf",&vR2andes);

	fscanf(entra_var, "%lf",&Kf);
	K_d = 0.0;
	fscanf(entra_var, "%lf",&K_m);

	fscanf(entra_var, "%lf",&ls);
	fscanf(entra_var, "%lf",&lb);
	fscanf(entra_var, "%lf",&lb2);

	fscanf(entra_var, "%lf",&uplift_scale);

    fscanf(entra_var, "%lf",&time_ofchangeu);
    fscanf(entra_var, "%lf",&uplift_scale2);
    fscanf(entra_var, "%lf",&tempo_max);
    fscanf(entra_var, "%lf",&dt);

	aloca_falhas();


	fclose(entra_var);


	//char nome[80];
	//FILE *f_perfil;


    time_t t1,t2;


    tri_p = Aloc_matrix_long (2000,4);
    cond_p = Aloc_vector_long (2000);
    aresta_p = Aloc_matrix_long (2000,5);


    malha();
    aloca_topo();
    aloca_falhas();
    aloca_difusao();
    aloca_fluvi();
	aloca_topo_din();
	
	printf("aloca topo_din ok\n");
    aloca_uplift();
	printf("aloca uplift ok\n");
	aloca_underp();
	printf("aloca underplate ok\n");
	read_vR_external();
	read_sea_level();
	vR_sort();


    aloca_falha_3D();

    flexuraK();


	///// Provisório
	for (i=0;i<nodes;i++){
		h_flex_cumulat[i]=0;
	}
	for (i=0;i<tri_flex;i++){
		h_flex_cumulat[Tri_flex[i][0]]=Te[i];
		h_flex_cumulat[Tri_flex[i][1]]=Te[i];
		h_flex_cumulat[Tri_flex[i][2]]=Te[i];
	}

	for (i=0;i<nodes;i++){
		h_flex_cumulat[i]=0;
	}
	//exit(-1);
	/////

    (void) time(&t1);


    for (soma_h=0,i=0;i<nodes;i++){
        soma_h+=h_topo[i]*area_vor[i];
    }
    printf("Area: %f\n",soma_h);
    printf("%f %f %f\n",pos_falha,despl,Edge);

    printf("\n");


    tempo=0;
    print_topo(tempo,0);

	printf("print_topo ok\n");

    monta_falha();

	printf("monta_falha ok\n");

	struct timeval start, stop;
	double secs;

	vR_calc();

	if (vR_external_flag==1) calc_vR_external();
	gettimeofday(&start, NULL);

    for (tempo=dt;tempo<=tempo_max;tempo+=dt){

		if (cont_falha<num_falha){
			//if (tempo-pos_falha_vec[cont_falha][2]<=0.4E6) Te_uc_on=1;
			if (tempo<=2000.0E6) Te_uc_on=1;
			else Te_uc_on=0;
		}

		if (n_sea_levels>0) sea_level_change();

        difusao();    //printf("d");
        fluvi_min();  //printf("f");
        modif_topo(); //printf("M");


        orog();			//printf("orog ok\n");
		underplate();
		topo_din();

        //monta_falha();

		if (long(tempo)%20000==0){
			if (vR_external_flag==-1) vR_calc_orography();
			if (vR_external_flag==0) vR_calc();
			if (vR_external_flag==1) calc_vR_external();
		}

		


        if (long(tempo)%(2000)==0){

			gettimeofday(&stop, NULL);
			secs = (double)(stop.tv_usec - start.tv_usec) / 1000000 + (double)(stop.tv_sec - start.tv_sec);
			printf("time taken %f\n",secs);

			flexura_load();
			flexura_solv(Kflex, Kdiag, bflex, wflex);

			monta_h_w();


			printf("Tempo: %f flexura \n",tempo);

			for (soma_h=0,i=0,j=0;i<nodes;i++){
                if (soma_h>h_bed[i] && h_topo[i]<nivel){
					soma_h=h_bed[i];
					j=i;
				}
            }
            printf("Tempo: %f Minimo : %g %g \n",tempo,h_bed[j],h_topo[j]);
			gettimeofday(&start, NULL);
        }

		fault3D();

		print_topo(tempo,0);
        /*
		if (long(tempo)%1000000==0){
			sprintf(nome,"perfil_FMHE_%.2fMa.txt",tempo/1E6); //topografia
			f_perfil = fopen(nome,"w");
			for (i=0; i<nodes; i++){
				if (xy[i][1]==0)
					fprintf(f_perfil,"%f %f %f\n", xy[i][0], h_topo[i], h_bed[i]);
			}
			fclose(f_perfil);
		}
		*/

    }
    (void) time(&t2);
    printf("\nTempo: %ld",t2-t1);
    free_all();
    return(0);
}

double **Aloc_matrix_real (long p, long n)
{
	double **v;
	long i;
	if(p < 1 || n<1){
		printf("** Erro: Parametro invalido **\n");
		return(NULL);
	}
	v = (double **) calloc(p, sizeof(double *));
	for(i=0; i<p; i++){
		v[i] = (double *) calloc (n, sizeof(double));
		if (v[i] == NULL) {
			printf("** Erro: Memoria Insuficiente **");
			return(NULL);
		}
	}
	return(v);
}
long **Aloc_matrix_long (long p, long n)
{
	long **v;
	long i;
	if(p < 1 || n<1){
		printf("** Erro: Parametro invalido **\n");
		return(NULL);
	}
	v = (long **) calloc(p, sizeof(long *));
	for(i=0; i<p; i++){
		v[i] = (long *) calloc (n, sizeof(long));
		if (v[i] == NULL) {
			printf("** Erro: Memoria Insuficiente **");
			return(NULL);
		}
	}
	return(v);
}
double *Aloc_vector_real (long n)
{
	double *v;
	v = (double *) calloc(n, sizeof(double));
	if (v == NULL){
		printf("* Erro: Memoria Insuficiente *");
		return(NULL);
	}
	return(v);
}

long *Aloc_vector_long (long n)
{
	long *v;
	v = (long *) calloc(n, sizeof(long));
	if (v == NULL){
		printf("* Erro: Memoria Insuficiente *");
		return(NULL);
	}
	return(v);
}

