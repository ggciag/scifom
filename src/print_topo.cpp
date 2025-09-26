void calc_Temper35();

#include <stdio.h>
#include <stdlib.h>

extern FILE *f_h_topo;
extern FILE *f_h_sed;


extern double Sed_esq;
extern double Sed_dir;
extern double Sed_in;

extern long nodes;
extern double **xy;
extern long tri;
extern long **Tri;

extern double *h_topo;
extern double *h_bed;
extern double *h_foot;
extern double *moho;

extern double *h_crust_sup;

extern double *h_flex_cumulat;
extern double *hu_flex_cumulat;

extern double *h_temper35;

extern long nodes_thermal;
extern double **xyz_thermal;
extern double *Temper;

extern double axis_stream;

extern double tgdip;
extern double V_meio;
extern double seg_per_ano;

extern double **pos_falha_vec;
extern long cont_falha;
extern long num_falha;

extern long *lagos;

extern long n_lat;
extern long *in_tri;

extern double tempo_max_stream;

extern double tempo_despl;

extern long Te_uc_on;


extern long *direc_fluvi;

extern double *Qr;

extern long n_sub_dt;

extern double dt;

extern double *vR_map;

extern double *area_vor;

extern int vR_external_flag;

void calc_vR_external(long t);

void fluvi_min();


void print_topo(double tempo,long muda_ponto)
{

	calc_Temper35();

	long i,t;

	char nome[30];

	if (long(tempo)%10000==0){

		FILE *Ftopo;

		sprintf(nome, "m_Topo_%.3f.txt",tempo/1e6);

		Ftopo = fopen(nome, "w");

		//for (i=0;i<nodes;i++) fprintf(Ftopo,"%10.3f %10.3f %10.3f %10.3f %10.4f %6.1f %6.1f %ld %1ld %.4lf\n",h_topo[i],h_bed[i],h_crust_sup[i],moho[i],h_flex_cumulat[i],h_temper35[i],(Qr[i]*n_sub_dt)/dt,direc_fluvi[i],lagos[i],vR_map[i]/area_vor[i]);
		for (i=0;i<nodes;i++) fprintf(Ftopo,"%10.3f %10.3f %6.1f %ld %1ld %.4lf\n",h_topo[i],h_bed[i],(Qr[i]*n_sub_dt)/dt,direc_fluvi[i],lagos[i],vR_map[i]/area_vor[i]);

		fclose(Ftopo);


		if (vR_external_flag==1){
			calc_vR_external(1);
			fluvi_min();
			sprintf(nome, "Prec_max_%.3f.txt",tempo/1e6);
			Ftopo = fopen(nome, "w");
			for (i=0;i<nodes;i++) fprintf(Ftopo,"%6.1f %.4lf\n",(Qr[i]*n_sub_dt)/dt,vR_map[i]/area_vor[i]);
			fclose(Ftopo);

			calc_vR_external(2);
			fluvi_min();
			sprintf(nome, "Prec_min_%.3f.txt",tempo/1e6);
			Ftopo = fopen(nome, "w");
			for (i=0;i<nodes;i++) fprintf(Ftopo,"%6.1f %.4lf\n",(Qr[i]*n_sub_dt)/dt,vR_map[i]/area_vor[i]);
			fclose(Ftopo);


			calc_vR_external(0);
		}
		fluvi_min();


		if (tempo==0){
			Ftopo = fopen("m_Tri.txt", "w");
			fprintf(Ftopo, "%ld\n", tri);
			for (t=0;t<tri;t++) fprintf(Ftopo,"%ld %ld %ld\n",Tri[t][0],Tri[t][1],Tri[t][2]);
			fclose(Ftopo);

			Ftopo = fopen("m_xy.txt","w");
			fprintf(Ftopo, "%ld\n", nodes);
			for (i=0;i<nodes;i++) fprintf(Ftopo,"%10.3f %10.3f\n",xy[i][0],xy[i][1]);
			fclose(Ftopo);
		}


		for (i=0;i<nodes;i++){
			h_flex_cumulat[i]=0;
			hu_flex_cumulat[i]=0;
		}

	}

}
