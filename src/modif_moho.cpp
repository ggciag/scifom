#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double var_fault_stream(double v);

extern double *moho;
extern double *h_crust_sup;

extern double *h_bed;

extern long **conec_Tri;
extern long *pos_conec_Tri;

extern long nodes;

extern double **xy;

extern double tempo_despl;
extern double V_meio;

extern double maxx;
extern double minx;
extern double maxy;
extern double miny;

extern double Edge;

extern double seg_per_ano;

extern long **Tri;

extern double axis_stream;

extern double *h_q;

extern double RHOC;
extern double RHOM;

double stream_x (double x,double z);
double stream_z (double x,double z);

int intriangle(double **pontos, long i, long j, long k, double xx, double yy);

extern long Nx;
extern long Ny;

extern double *moho_flex;
extern double *moho_flex_aux;
extern double **xy_flex;

extern double dx_flex;
extern double dy_flex;

extern double **peso;
extern long **peso_pos;

extern double *temper_q;



void modif_moho()
{
    
	
	long i;
	double z;
	
	double dx,dz;
    
	double x_mean,y_mean,z_mean;
	
	for (i=0;i<nodes;i++){
		
		x_mean = xy[i][0];
		y_mean = xy[i][1];
		z_mean = moho[i];
		
		//if (moho[i]>moho_max_antes) moho_max_antes=moho[i];

		
		dx=stream_x(x_mean-(axis_stream+var_fault_stream(y_mean)),z_mean)*tempo_despl*seg_per_ano;
		dz=stream_z(x_mean-(axis_stream+var_fault_stream(y_mean)),z_mean)*tempo_despl*seg_per_ano;
		
		//if (fabs(dz)>zmax) zmax=fabs(dz);
		
		z = moho[i]+dz;
		
		if (z>h_crust_sup[i]) z=h_crust_sup[i];
		h_q[i]+=-(z-moho[i])*(RHOM-RHOC);
		moho[i]=z;
		
		//if (moho[i]>moho_max) moho_max=moho[i];
	}
	
	//printf("dz_max = %lf  %lf   %lf\n",zmax,moho_max_antes,moho_max);
      
}
