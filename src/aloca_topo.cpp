//modif: ver moho flex!!!

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

extern FILE *f_h_topo;
extern FILE *f_h_sed;

double **Aloc_matrix_real (long p, long n);
double *Aloc_vector_real (long n);
long *Aloc_vector_long (long n);

double var_fault_stream(double v);

extern double nivel;

extern double **xy;

extern double *h_topo;

extern double *f_difusao;
extern double *h_bed;
extern double *h_flex_cumulat;
extern double *hu_flex_cumulat;
extern double *water;

extern double *h_temper35;


extern double *moho;
extern double *load_Temper;

extern long *mar;


//double *h_isost;
extern double *h_q;
extern double *h_q_a;
extern double *h_w;
extern double *h_u;

extern long nodes;
extern long nodes_max_aloca;

extern double *h_foot;

extern double *h_crust_sup;

extern double *moho_aux;

extern double *Lf_vec;


extern double pos_falha;
extern double inclina;

extern double tgdip;

extern double H_C;
extern double H_brit;

extern double RHOW;
extern double RHOS;
extern double RHOC;
extern double RHOM;

extern double depth;

extern long *side;

extern double topo_seno;

extern double minx;
extern double maxx;
extern double miny;

extern double axis_stream;


double var_fault(double v);


extern long nodes_flex;
extern double **xy_flex;
extern double *moho_flex;

extern long **peso_pos;
extern double **peso;


extern double *temper_q_a;

extern long layers;
extern double **xyz_thermal;
extern double *T_vec_fut;

extern double alpha_exp_thermo;
/*
extern double *Te;
extern long tri_flex;
*/
extern double *Te_map;

void aloca_topo()
{
	printf("nodes_max_aloca: %ld, nodes: %ld\n",nodes_max_aloca,nodes);
	long i,j;
	//double dist2,dist_max=40000;
    //double aux;
	h_topo = Aloc_vector_real (nodes_max_aloca);
	h_crust_sup = Aloc_vector_real (nodes_max_aloca);

	Te_map = Aloc_vector_real (nodes_max_aloca);
	//Te = Aloc_vector_real (tri_flex);

	h_bed = Aloc_vector_real (nodes_max_aloca);
	h_flex_cumulat = Aloc_vector_real (nodes_max_aloca);
	hu_flex_cumulat = Aloc_vector_real (nodes_max_aloca);
	h_q = Aloc_vector_real (nodes_max_aloca);
	h_q_a = Aloc_vector_real (nodes_max_aloca);
	h_w = Aloc_vector_real (nodes_max_aloca);
	h_u = Aloc_vector_real (nodes_max_aloca);
	moho = Aloc_vector_real (nodes_max_aloca);
	load_Temper = Aloc_vector_real (nodes_max_aloca);

	water = Aloc_vector_real (nodes_max_aloca);

	h_temper35 = Aloc_vector_real (nodes_max_aloca);

	h_foot = Aloc_vector_real (nodes_max_aloca);

	Lf_vec = Aloc_vector_real(nodes_max_aloca);


	moho_aux = Aloc_vector_real (nodes_max_aloca);

	side = Aloc_vector_long(nodes_max_aloca);

	//double x_aux,y_aux;


	//double fundo=-7000.0;

	mar = Aloc_vector_long(nodes_max_aloca);




	f_h_topo = fopen("topo_historia.txt","w");
	f_h_sed = fopen("sed_historia.txt","w");


	/*double xh[20][2];
	long cont_xh=7;

	xh[0][0] =-200000.0+50000.0; xh[0][1]=0.0;
	xh[1][0] = 100000.0+50000.0; xh[1][1]=1.0;
	xh[2][0] = 110000.0+50000.0; xh[2][1]=1500.0;
	xh[3][0] = 300000.0+50000.0; xh[3][1]=1600.0;
	xh[4][0] = 490000.0+50000.0; xh[4][1]=1500.0;
	xh[5][0] = 500000.0+50000.0; xh[5][1]=1.0;
	xh[6][0] = 1000000.0+50000.0; xh[6][1]=0.0;*/

	FILE *f_topo;

	f_topo = fopen("topo_moho_lito_Te.txt","r");


	for (i=0;i<nodes;i++){

		//h_topo[i]=1000.0;

		/*if (xy[i][0]<=xh[0][0])
			h_topo[i]=xh[0][1];

		if (xy[i][0]>xh[cont_xh-1][0])
			h_topo[i]=xh[cont_xh-1][1];

		for (j=0;j<cont_xh-1;j++){
			if (xy[i][0]>xh[j][0] && xy[i][0]<=xh[j+1][0]){
				h_topo[i] = xh[j][1]*(xh[j+1][0]-xy[i][0])/(xh[j+1][0]-xh[j][0]);
				h_topo[i]+= xh[j+1][1]*(xy[i][0]-xh[j][0])/(xh[j+1][0]-xh[j][0]);
			}
		}*/

		fscanf(f_topo,"%lf %lf %lf %lf",&h_topo[i],&moho[i],&Lf_vec[i],&Te_map[i]);

		h_bed[i]=h_topo[i];

		if (h_topo[i]<nivel){
			mar[i]=1;
		}


		//if (h_topo[i]>200.0) Lf_vec[i]=4000000.0;


		//moho[i]=-H_C;
		h_crust_sup[i] = h_topo[i]-H_brit;

		if (h_crust_sup[i]<moho[i]) h_crust_sup[i]=moho[i];

		if (h_topo[i]<nivel) water[i]=nivel-h_topo[i];
	}

	fclose(f_topo);

	printf("topo_moho_litho_up: done\n");

	double dist;
	double dx;
	double dy;

	for (i=0;i<nodes_flex;i++){

		dist=1E11;
		for (j=0;j<nodes;j++){
			dx = xy[j][0]-xy_flex[i][0];
			dy = xy[j][1]-xy_flex[i][1];
			if (dist>dx*dx+dy*dy){
				moho_flex[i]=moho[j];
				dist=dx*dx+dy*dy;
			}
		}
	}

	/*for (i=0;i<nodes;i++){
		h_w[i] =peso[i][0]*moho_flex[peso_pos[i][0]];
		h_w[i]+=peso[i][1]*moho_flex[peso_pos[i][1]];
		h_w[i]+=peso[i][2]*moho_flex[peso_pos[i][2]];

		moho[i]=h_w[i];
	}*/

	double h_aux;
	double T_aux;
	double z;

	for (i=0;i<nodes_flex;i++){
		temper_q_a[i]=0.0;
	}

	for (j=0;j<layers-1;j++){
		for (i=0;i<nodes_flex;i++){
			h_aux = xyz_thermal[j*nodes_flex+i][2]-xyz_thermal[(j+1)*nodes_flex+i][2];
			T_aux = (T_vec_fut[j*nodes_flex+i]+T_vec_fut[(j+1)*nodes_flex+i])/2;
			temper_q_a[i]+=h_aux*RHOM*(1.0-alpha_exp_thermo*T_aux);
		}
	}

	for (i=0;i<nodes;i++){
		z =peso[i][0]*temper_q_a[peso_pos[i][0]];
		z+=peso[i][1]*temper_q_a[peso_pos[i][1]];
		z+=peso[i][2]*temper_q_a[peso_pos[i][2]];

		load_Temper[i]=z;
	}




}
