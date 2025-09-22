#include <stdio.h>
#include <stdlib.h>
#define NODE_POR_ELE 3
#define GLNODE 3

void montaf (double *load,double xi,double xj,double xk,double yi,double yj,double yk, 
             double q1, long ni, long nj, long nk,double L1,double L2,double L3);

void montaf3(double *load,double xi,double xj,double xk,double yi,
			 double yj,double yk, double q1, double q2, double q3, long ni, long nj, long nk);

extern long **Tri;
extern long tri;
extern double **xy;

extern double nivel;

extern double *area_vor;
extern long nodes;

extern double *h_topo;
extern double *h_bed;
extern double *water;
//extern double *h_original;


extern double *h_q;
extern double *h_q_a;

extern double *bflex;

extern long tri_flex;
extern long **Tri_flex;
extern long nodes_flex;
extern double **xy_flex;

extern double area_ele_flex;

extern double **peso;
extern long **peso_pos;

extern double RHOW;
extern double RHOS;
extern double RHOC;
extern double RHOM;
extern double alpha_thermal;

extern double *Temper;

extern long *pos_temper;
extern double **xyz_thermal;

extern double *temper_q;
extern double *temper_q_a;

extern double *moho;

extern long nodes_thermal;

extern double H_C;

extern double depth;

extern long Nx;
extern long Ny;

extern double *load_Temper;
extern long layers;

extern double *T_vec;
extern double alpha_exp_thermo;

extern double minx;
extern double maxx;


void flexura_load()
{
    long i,fi,j;

	double xi,xj,xk,yi,yj,yk,tq1;
    
    for (fi=0;fi<nodes_flex;fi++){
        for (j=0;j<3;j++){
            bflex[fi*3+j]=0;
        }
    }
	
	
	
	double h_aux;
	double T_aux;
	double z;
	
	for (i=0;i<nodes_flex;i++){
		temper_q_a[i]=0.0;
	}
	
	for (j=0;j<layers-1;j++){
		for (i=0;i<nodes_flex;i++){
			h_aux = xyz_thermal[j*nodes_flex+i][2]-xyz_thermal[(j+1)*nodes_flex+i][2];
			T_aux = (T_vec[j*nodes_flex+i]+T_vec[(j+1)*nodes_flex+i])/2;
			temper_q_a[i]+=h_aux*RHOM*(1.0-alpha_exp_thermo*T_aux);
		}
	}
	
	/////Calcula diferença
	for (i=0;i<nodes;i++){
		z =peso[i][0]*temper_q_a[peso_pos[i][0]];
		z+=peso[i][1]*temper_q_a[peso_pos[i][1]];
		z+=peso[i][2]*temper_q_a[peso_pos[i][2]];
		
		load_Temper[i]=load_Temper[i]-z;
	}
	
	
	//////Loads
	for (i=0;i<nodes;i++){
        xi=xy_flex[peso_pos[i][0]][0];
		xj=xy_flex[peso_pos[i][1]][0];
		xk=xy_flex[peso_pos[i][2]][0];
		yi=xy_flex[peso_pos[i][0]][1];
		yj=xy_flex[peso_pos[i][1]][1];
		yk=xy_flex[peso_pos[i][2]][1];
		tq1=0;
		tq1=(h_q[i]+h_q_a[i])*9.8*area_vor[i];
		tq1+=(load_Temper[i])*9.8*area_vor[i];

		montaf (bflex,xi,xj,xk,yi,yj,yk,tq1,
				peso_pos[i][0],peso_pos[i][1],peso_pos[i][2],
				peso[i][0],peso[i][1],peso[i][2]);
    }
	
	for (i=0;i<nodes;i++){
		h_q[i]=0;
	}
	
	for (i=0;i<nodes;i++){
		z =peso[i][0]*temper_q_a[peso_pos[i][0]];
		z+=peso[i][1]*temper_q_a[peso_pos[i][1]];
		z+=peso[i][2]*temper_q_a[peso_pos[i][2]];
		
		load_Temper[i]=z;
		if (i==100) printf("load_Temper = %lf\n",load_Temper[i]);
	}
	
	///// Calcula carga fora
	
	long verif,i_min=0,i_max=Nx-1,ii;
	j=0;
	for (i=0,verif=0;i<Nx && verif==0;i++){
		if (xy_flex[i*Ny+j][0]>minx){
			i_min=i;
			verif=1;
			if (i-1>=0) i_min=i-1;
		}	
	}
	for (i=Nx-1,verif=0;i>=0 && verif==0;i--){
		if (xy_flex[i*Ny+j][0]<maxx){
			i_max=i;
			if (i+1<Nx) i_max=i+1;
			verif=1;
		}	
	}
	
	for (i=i_min-1;i>=0;i--){
		ii=i_min;
		//ii=(i_min-i)+i_min-1;
		if (ii<Nx){
			bflex[(i*Ny+j)*3]=bflex[(ii*Ny+j)*3];
			bflex[(i*Ny+j)*3+1]=bflex[(ii*Ny+j)*3+1];
			bflex[(i*Ny+j)*3+2]=bflex[(ii*Ny+j)*3+2];
		}
	}
	
	for (i=i_max+1;i<Nx;i++){
		ii=i_max;
		//ii=(i_max-i)+i_max+1;
		if (ii>=0){
			bflex[(i*Ny+j)*3]=bflex[(ii*Ny+j)*3];
			bflex[(i*Ny+j)*3+1]=bflex[(ii*Ny+j)*3+1];
			bflex[(i*Ny+j)*3+2]=bflex[(ii*Ny+j)*3+2];
		}
	}
	
	
	for (fi=0;fi<nodes_flex;fi++){
		temper_q[fi]=0;
	}
	
    
    
    
    /**/
    
}

void montaf (double *load,double xi,double xj,double xk,double yi,double yj,double yk, 
             double q1, long ni, long nj, long nk,double L1,double L2,double L3)
{
	double N[9];
	double b1,b2,b3,c1,c2,c3,l1,l2,l3,mi1,mi2,mi3,delta;
	b1=yj-yk; c1=xk-xj;
	b2=yk-yi; c2=xi-xk;
	b3=yi-yj; c3=xj-xi;
	l1=(xk-xj)*(xk-xj)+(yk-yj)*(yk-yj);
	l2=(xk-xi)*(xk-xi)+(yk-yi)*(yk-yi);
	l3=(xi-xj)*(xi-xj)+(yi-yj)*(yi-yj);
	mi1=(l3-l2)/l1;
	mi2=(l1-l3)/l2;
	mi3=(l2-l1)/l3;
	delta=0.5*((xj*yk-yj*xk)+xi*(yj-yk)+yi*(xk-xj));
	if (delta==0) {
        printf("Problemas! delta=0!"); 
        printf("%f %f %f %f %f %f\n",xi,xj,xk, yi,yj,yk);
        printf("%ld %ld %ld",ni,nj,nk);
        exit(-1);
    }
	if (delta<0) {printf("Que estranho"); delta=delta*(-1); }


	N[0]=2*(-L1*L3*L3+L1*L2*L3*((3*mi3+1)*L3-(3*mi3+1)*L2+3*(1-mi3)*L1)/2-L1*L2*L3*(3*(1-mi2)*L3+(3*mi2+1)*L2-(3*mi2+1)*L1)/2+L1*L1*L2)+L1*L3-L1*L2+L1;
	N[1]=-b2*(L1*L3*L3+L1*L2*L3*(3*(1-mi2)*L3+(3*mi2+1)*L2-(3*mi2+1)*L1)/2-L1*L3)-b3*(L1*L2*L3*((3*mi3+1)*L3-(3*mi3+1)*L2+3*(1-mi3)*L1)/2+L1*L1*L2);
	N[2]=-c2*(L1*L3*L3+L1*L2*L3*(3*(1-mi2)*L3+(3*mi2+1)*L2-(3*mi2+1)*L1)/2-L1*L3)-c3*(L1*L2*L3*((3*mi3+1)*L3-(3*mi3+1)*L2+3*(1-mi3)*L1)/2+L1*L1*L2);
	N[3]=2*(-L1*L2*L3*((3*mi3+1)*L3-(3*mi3+1)*L2+3*(1-mi3)*L1)/2+L1*L2*L3*(-(3*mi1+1)*L3+3*(1-mi1)*L2+(3*mi1+1)*L1)/2+L2*L2*L3-L1*L1*L2)-L2*L3+L1*L2+L2;
	N[4]=-b3*(L1*L2*L3*((3*mi3+1)*L3-(3*mi3+1)*L2+3*(1-mi3)*L1)/2+L1*L1*L2-L1*L2)-b1*(L1*L2*L3*(-(3*mi1+1)*L3+3*(1-mi1)*L2+(3*mi1+1)*L1)/2+L2*L2*L3);
	N[5]=-c3*(L1*L2*L3*((3*mi3+1)*L3-(3*mi3+1)*L2+3*(1-mi3)*L1)/2+L1*L1*L2-L1*L2)-c1*(L1*L2*L3*(-(3*mi1+1)*L3+3*(1-mi1)*L2+(3*mi1+1)*L1)/2+L2*L2*L3);
	N[6]=2*(L1*L3*L3+L1*L2*L3*(3*(1-mi2)*L3+(3*mi2+1)*L2-(3*mi2+1)*L1)/2-L1*L2*L3*(-(3*mi1+1)*L3+3*(1-mi1)*L2+(3*mi1+1)*L1)/2-L2*L2*L3)+L2*L3-L1*L3+L3;
	N[7]=-b2*(L1*L3*L3+L1*L2*L3*(3*(1-mi2)*L3+(3*mi2+1)*L2-(3*mi2+1)*L1)/2)-b1*(L1*L2*L3*(-(3*mi1+1)*L3+3*(1-mi1)*L2+(3*mi1+1)*L1)/2+L2*L2*L3-L2*L3);
	N[8]=-c2*(L1*L3*L3+L1*L2*L3*(3*(1-mi2)*L3+(3*mi2+1)*L2-(3*mi2+1)*L1)/2)-c1*(L1*L2*L3*(-(3*mi1+1)*L3+3*(1-mi1)*L2+(3*mi1+1)*L1)/2+L2*L2*L3-L2*L3);
	
	load[ni*3+0]+= N[0]*q1;
	load[ni*3+1]+= N[1]*q1;
	load[ni*3+2]+= N[2]*q1;
	load[nj*3+0]+= N[3]*q1;
	load[nj*3+1]+= N[4]*q1;
	load[nj*3+2]+= N[5]*q1;
	load[nk*3+0]+= N[6]*q1;
	load[nk*3+1]+= N[7]*q1;
	load[nk*3+2]+= N[8]*q1; /**/
  	
	
	/*for(cont=0,q1=q1/delta; cont<3; cont++){
		
		if(cont==0) { L1=0; 	L2=0.5; 	L3=0.5; }
		if(cont==1) { L1=0.5; 	L2=0; 		L3=0.5; }
		if(cont==2) { L1=0.5; 	L2=0.5; 	L3=0; }
		
		N[0]=2*(-L1*L3*L3+L1*L2*L3*((3*mi3+1)*L3-(3*mi3+1)*L2+3*(1-mi3)*L1)/2-L1*L2*L3*(3*(1-mi2)*L3+(3*mi2+1)*L2-(3*mi2+1)*L1)/2+L1*L1*L2)+L1*L3-L1*L2+L1;
    	N[1]=-b2*(L1*L3*L3+L1*L2*L3*(3*(1-mi2)*L3+(3*mi2+1)*L2-(3*mi2+1)*L1)/2-L1*L3)-b3*(L1*L2*L3*((3*mi3+1)*L3-(3*mi3+1)*L2+3*(1-mi3)*L1)/2+L1*L1*L2);
		N[2]=-c2*(L1*L3*L3+L1*L2*L3*(3*(1-mi2)*L3+(3*mi2+1)*L2-(3*mi2+1)*L1)/2-L1*L3)-c3*(L1*L2*L3*((3*mi3+1)*L3-(3*mi3+1)*L2+3*(1-mi3)*L1)/2+L1*L1*L2);
		N[3]=2*(-L1*L2*L3*((3*mi3+1)*L3-(3*mi3+1)*L2+3*(1-mi3)*L1)/2+L1*L2*L3*(-(3*mi1+1)*L3+3*(1-mi1)*L2+(3*mi1+1)*L1)/2+L2*L2*L3-L1*L1*L2)-L2*L3+L1*L2+L2;
		N[4]=-b3*(L1*L2*L3*((3*mi3+1)*L3-(3*mi3+1)*L2+3*(1-mi3)*L1)/2+L1*L1*L2-L1*L2)-b1*(L1*L2*L3*(-(3*mi1+1)*L3+3*(1-mi1)*L2+(3*mi1+1)*L1)/2+L2*L2*L3);
		N[5]=-c3*(L1*L2*L3*((3*mi3+1)*L3-(3*mi3+1)*L2+3*(1-mi3)*L1)/2+L1*L1*L2-L1*L2)-c1*(L1*L2*L3*(-(3*mi1+1)*L3+3*(1-mi1)*L2+(3*mi1+1)*L1)/2+L2*L2*L3);
		N[6]=2*(L1*L3*L3+L1*L2*L3*(3*(1-mi2)*L3+(3*mi2+1)*L2-(3*mi2+1)*L1)/2-L1*L2*L3*(-(3*mi1+1)*L3+3*(1-mi1)*L2+(3*mi1+1)*L1)/2-L2*L2*L3)+L2*L3-L1*L3+L3;
		N[7]=-b2*(L1*L3*L3+L1*L2*L3*(3*(1-mi2)*L3+(3*mi2+1)*L2-(3*mi2+1)*L1)/2)-b1*(L1*L2*L3*(-(3*mi1+1)*L3+3*(1-mi1)*L2+(3*mi1+1)*L1)/2+L2*L2*L3-L2*L3);
		N[8]=-c2*(L1*L3*L3+L1*L2*L3*(3*(1-mi2)*L3+(3*mi2+1)*L2-(3*mi2+1)*L1)/2)-c1*(L1*L2*L3*(-(3*mi1+1)*L3+3*(1-mi1)*L2+(3*mi1+1)*L1)/2+L2*L2*L3-L2*L3);
		
		load[ni*3+0] = load[ni*3+0] + N[0]*(q1)*delta/3;
		load[ni*3+1] = load[ni*3+1] + N[1]*(q1)*delta/3;
		load[ni*3+2] = load[ni*3+2] + N[2]*(q1)*delta/3;
		load[nj*3+0] = load[nj*3+0] + N[3]*(q1)*delta/3;
		load[nj*3+1] = load[nj*3+1] + N[4]*(q1)*delta/3;
		load[nj*3+2] = load[nj*3+2] + N[5]*(q1)*delta/3;
		load[nk*3+0] = load[nk*3+0] + N[6]*(q1)*delta/3;
		load[nk*3+1] = load[nk*3+1] + N[7]*(q1)*delta/3;
		load[nk*3+2] = load[nk*3+2] + N[8]*(q1)*delta/3;
	}*/



}


void montaf3(double *load,double xi,double xj,double xk,double yi,
			double yj,double yk, double q1, double q2, double q3, long ni, long nj, long nk)
{
	long cont;
	double N[9], L1, L2, L3;
	double b1,b2,b3,c1,c2,c3,l1,l2,l3,mi1,mi2,mi3,delta;
	b1=yj-yk; c1=xk-xj;
	b2=yk-yi; c2=xi-xk;
	b3=yi-yj; c3=xj-xi;
	l1=(xk-xj)*(xk-xj)+(yk-yj)*(yk-yj);
	l2=(xk-xi)*(xk-xi)+(yk-yi)*(yk-yi);
	l3=(xi-xj)*(xi-xj)+(yi-yj)*(yi-yj);
	mi1=(l3-l2)/l1;
	mi2=(l1-l3)/l2;
	mi3=(l2-l1)/l3;
	delta=0.5*((xj*yk-yj*xk)+xi*(yj-yk)+yi*(xk-xj));
	if (delta==0) exit(-1);
	if (delta<0) delta=delta*(-1);
	
	for(cont=0; cont<3; cont++){
		
		if(cont==0) { L1=0; 	L2=0.5; 	L3=0.5; }
		if(cont==1) { L1=0.5; 	L2=0; 		L3=0.5; }
		if(cont==2) { L1=0.5; 	L2=0.5; 	L3=0; }
		
		N[0]=2*(-L1*L3*L3+L1*L2*L3*((3*mi3+1)*L3-(3*mi3+1)*L2+3*(1-mi3)*L1)/2-L1*L2*L3*(3*(1-mi2)*L3+(3*mi2+1)*L2-(3*mi2+1)*L1)/2+L1*L1*L2)+L1*L3-L1*L2+L1;
    	N[1]=-b2*(L1*L3*L3+L1*L2*L3*(3*(1-mi2)*L3+(3*mi2+1)*L2-(3*mi2+1)*L1)/2-L1*L3)-b3*(L1*L2*L3*((3*mi3+1)*L3-(3*mi3+1)*L2+3*(1-mi3)*L1)/2+L1*L1*L2);
		N[2]=-c2*(L1*L3*L3+L1*L2*L3*(3*(1-mi2)*L3+(3*mi2+1)*L2-(3*mi2+1)*L1)/2-L1*L3)-c3*(L1*L2*L3*((3*mi3+1)*L3-(3*mi3+1)*L2+3*(1-mi3)*L1)/2+L1*L1*L2);
		N[3]=2*(-L1*L2*L3*((3*mi3+1)*L3-(3*mi3+1)*L2+3*(1-mi3)*L1)/2+L1*L2*L3*(-(3*mi1+1)*L3+3*(1-mi1)*L2+(3*mi1+1)*L1)/2+L2*L2*L3-L1*L1*L2)-L2*L3+L1*L2+L2;
		N[4]=-b3*(L1*L2*L3*((3*mi3+1)*L3-(3*mi3+1)*L2+3*(1-mi3)*L1)/2+L1*L1*L2-L1*L2)-b1*(L1*L2*L3*(-(3*mi1+1)*L3+3*(1-mi1)*L2+(3*mi1+1)*L1)/2+L2*L2*L3);
		N[5]=-c3*(L1*L2*L3*((3*mi3+1)*L3-(3*mi3+1)*L2+3*(1-mi3)*L1)/2+L1*L1*L2-L1*L2)-c1*(L1*L2*L3*(-(3*mi1+1)*L3+3*(1-mi1)*L2+(3*mi1+1)*L1)/2+L2*L2*L3);
		N[6]=2*(L1*L3*L3+L1*L2*L3*(3*(1-mi2)*L3+(3*mi2+1)*L2-(3*mi2+1)*L1)/2-L1*L2*L3*(-(3*mi1+1)*L3+3*(1-mi1)*L2+(3*mi1+1)*L1)/2-L2*L2*L3)+L2*L3-L1*L3+L3;
		N[7]=-b2*(L1*L3*L3+L1*L2*L3*(3*(1-mi2)*L3+(3*mi2+1)*L2-(3*mi2+1)*L1)/2)-b1*(L1*L2*L3*(-(3*mi1+1)*L3+3*(1-mi1)*L2+(3*mi1+1)*L1)/2+L2*L2*L3-L2*L3);
		N[8]=-c2*(L1*L3*L3+L1*L2*L3*(3*(1-mi2)*L3+(3*mi2+1)*L2-(3*mi2+1)*L1)/2)-c1*(L1*L2*L3*(-(3*mi1+1)*L3+3*(1-mi1)*L2+(3*mi1+1)*L1)/2+L2*L2*L3-L2*L3);
		
		load[ni*3+0] = load[ni*3+0] + N[0]*(q1*L1+q2*L2+q3*L3)*delta/3;
		load[ni*3+1] = load[ni*3+1] + N[1]*(q1*L1+q2*L2+q3*L3)*delta/3;
		load[ni*3+2] = load[ni*3+2] + N[2]*(q1*L1+q2*L2+q3*L3)*delta/3;
		load[nj*3+0] = load[nj*3+0] + N[3]*(q1*L1+q2*L2+q3*L3)*delta/3;
		load[nj*3+1] = load[nj*3+1] + N[4]*(q1*L1+q2*L2+q3*L3)*delta/3;
		load[nj*3+2] = load[nj*3+2] + N[5]*(q1*L1+q2*L2+q3*L3)*delta/3;
		load[nk*3+0] = load[nk*3+0] + N[6]*(q1*L1+q2*L2+q3*L3)*delta/3;
		load[nk*3+1] = load[nk*3+1] + N[7]*(q1*L1+q2*L2+q3*L3)*delta/3;
		load[nk*3+2] = load[nk*3+2] + N[8]*(q1*L1+q2*L2+q3*L3)*delta/3;
	}
	
}
