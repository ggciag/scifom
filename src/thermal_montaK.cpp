/*
 *  thermal_montaK.cpp
 *  thermal1
 *
 *  Created by Victor Sacek on 15/03/10.
 *  Copyright 2010 IAG/USP. All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

extern double **Kthermal;
extern double *Kthermal_diag;
extern double **Mthermal;
extern double *Mthermal_diag;
extern double **Kthermal2;
extern double *Kthermal2_diag;

extern long **Kthermal_conec;
extern long *Kthermal_posconec;

extern double **TKe;
extern double **TMe;

extern long tetra_thermal;
extern long **Tet_thermal;
extern long nodes_thermal;
extern double **xyz_thermal_fut;

extern double dt_calor_sec;

extern double seg_per_ano;

extern double kappa;

extern double alpha_thermal;

extern double soma_volume;

extern double **v_adv;

void montaKeThermal(double **Ke, double **Me, long **Tet_thermal, double **xyz_thermal, long t);

void thermal_monta_K()
{
	long t,i,ii,jj,iaux,pos_aux;
	
	for (i=0;i<nodes_thermal;i++){
		Mthermal_diag[i]=0.0; 
		Kthermal2_diag[i]=0.0;
		pos_aux = Kthermal_posconec[i];
		for (iaux=0;iaux<pos_aux;iaux++){
			Mthermal[i][iaux]=0.0;
			Kthermal2[i][iaux]=0.0;
		}			
	}
	
	soma_volume=0.0;
	
	for (t=0; t<tetra_thermal; t++){
		montaKeThermal(TKe,TMe,Tet_thermal,xyz_thermal_fut,t);			
		for (ii=0;ii<4;ii++){
			pos_aux = Tet_thermal[t][ii];
			for (jj=0;jj<4;jj++){				
				if (ii==jj){
					Kthermal2_diag[pos_aux]+=TKe[ii][jj];
					Mthermal_diag[pos_aux]+=TMe[ii][jj];
				}   
				else {
					for (iaux=0;iaux<Kthermal_posconec[pos_aux];iaux++){
						if (Kthermal_conec[pos_aux][iaux]==Tet_thermal[t][jj]){
							Kthermal2[pos_aux][iaux]+=TKe[ii][jj];
							Mthermal[pos_aux][iaux]+=TMe[ii][jj];
							iaux=4000;//numero bem grande para sair do for
						}
					}
				}                                  
			}
		}
	}	
	
	long estranho=0;
	for (i=0;i<nodes_thermal;i++){
		Kthermal_diag[i] = Mthermal_diag[i] + alpha_thermal*dt_calor_sec*Kthermal2_diag[i];
		if (Kthermal_diag[i]>0 || Kthermal_diag[i]<=0);
		else estranho=1;
		pos_aux = Kthermal_posconec[i];
		for (iaux=0;iaux<pos_aux;iaux++){
			Kthermal[i][iaux] = Mthermal[i][iaux] + alpha_thermal*dt_calor_sec*Kthermal2[i][iaux];
		}			
	}
	printf("Kthermal_diag: %ld",estranho);
	
	FILE *F_K;
    
    F_K = fopen("F_K.txt","w");
    for (t=0;t<nodes_thermal;t++){   
        fprintf(F_K,"%ld | %g | %g |",Kthermal_posconec[t],Mthermal_diag[t],Kthermal2_diag[t]);
        for (i=0;i<Kthermal_posconec[t];i++){
            fprintf(F_K,"%10.2g ",Kthermal[t][i]);
        }
        fprintf(F_K,"\n");
    }
    fclose(F_K);
	
	printf("Volume: %g Esperado: %g\n",soma_volume,9*9*9*1.0E9);
	
}

void montaKeThermal(double **Ke, double **Me, long **Tet_thermal, double **xyz_thermal, long t)
{
	long i,j;
	long aux_l;
	
	double volume;
	double aa[4],bb[4],cc[4],dd[4];
	double x[4],y[4],z[4],aux;
	double v_prov;//0*500.0/1.0E6/seg_per_ano;
	double vv[3];
	
	
	x[0]=xyz_thermal[Tet_thermal[t][0]][0];
	y[0]=xyz_thermal[Tet_thermal[t][0]][1];
	z[0]=xyz_thermal[Tet_thermal[t][0]][2];

	x[1]=xyz_thermal[Tet_thermal[t][1]][0];
	y[1]=xyz_thermal[Tet_thermal[t][1]][1];
	z[1]=xyz_thermal[Tet_thermal[t][1]][2];
	
	x[2]=xyz_thermal[Tet_thermal[t][2]][0];
	y[2]=xyz_thermal[Tet_thermal[t][2]][1];
	z[2]=xyz_thermal[Tet_thermal[t][2]][2];
	
	x[3]=xyz_thermal[Tet_thermal[t][3]][0];
	y[3]=xyz_thermal[Tet_thermal[t][3]][1];
	z[3]=xyz_thermal[Tet_thermal[t][3]][2];
	
	v_prov = v_adv[Tet_thermal[t][0]][0];
	v_prov+= v_adv[Tet_thermal[t][1]][0];
	v_prov+= v_adv[Tet_thermal[t][2]][0];
	v_prov+= v_adv[Tet_thermal[t][3]][0];
	v_prov/=4.0;
	vv[0]=v_prov; 
	
	v_prov = v_adv[Tet_thermal[t][0]][1];
	v_prov+= v_adv[Tet_thermal[t][1]][1];
	v_prov+= v_adv[Tet_thermal[t][2]][1];
	v_prov+= v_adv[Tet_thermal[t][3]][1];
	v_prov/=4.0;
	vv[1]=v_prov;
	
	v_prov = v_adv[Tet_thermal[t][0]][2];
	v_prov+= v_adv[Tet_thermal[t][1]][2];
	v_prov+= v_adv[Tet_thermal[t][2]][2];
	v_prov+= v_adv[Tet_thermal[t][3]][2];
	v_prov/=4.0;
	vv[2]=v_prov;
	
	aa[0]= x[1]*y[2]*z[3]+x[3]*y[1]*z[2]+x[2]*y[3]*z[1]-x[3]*y[2]*z[1]-x[1]*y[3]*z[2]-x[2]*y[1]*z[3];
	bb[0]=      y[2]*z[3]+     y[1]*z[2]+     y[3]*z[1]-     y[2]*z[1]-     y[3]*z[2]-     y[1]*z[3];
	cc[0]= -(x[1]     *z[3]+x[3]     *z[2]+x[2]     *z[1]-x[3]     *z[1]-x[1]     *z[2]-x[2]     *z[3]);
	dd[0]= x[1]*y[2]     +x[3]*y[1]     +x[2]*y[3]     -x[3]*y[2]     -x[1]*y[3]     -x[2]*y[1]     ;		
	
	
	volume=(aa[0]  -bb[0]*x[0]  +cc[0]*y[0]  -dd[0]*z[0])/6;
	
	if (volume<0) {
		volume*=-1.0;
	}
	else{
		aux=x[0]; x[0]=x[1]; x[1]=aux;
		aux=y[0]; y[0]=y[1]; y[1]=aux;
		aux=z[0]; z[0]=z[1]; z[1]=aux;
		aux_l = Tet_thermal[t][0]; Tet_thermal[t][0]=Tet_thermal[t][1]; Tet_thermal[t][1]=aux_l;		
	}
	
	soma_volume+=volume;
	
	double x24 = x[2-1]-x[4-1],x34 = x[3-1]-x[4-1],x14 = x[1-1]-x[4-1];	
	double y24 = y[2-1]-y[4-1],y34 = y[3-1]-y[4-1],y14 = y[1-1]-y[4-1];
	double z24 = z[2-1]-z[4-1],z34 = z[3-1]-z[4-1],z14 = z[1-1]-z[4-1];
	
	bb[0]=y24*z34-y34*z24; bb[1]=y34*z14-y14*z34; bb[2]=y14*z24-y24*z14;
	cc[0]=z24*x34-z34*x24; cc[1]=z34*x14-z14*x34; cc[2]=z14*x24-z24*x14;
	dd[0]=x24*y34-x34*y24; dd[1]=x34*y14-x14*y34; dd[2]=x14*y24-x24*y14;
	
	for (i=0;i<3;i++){
		bb[i]*=1.0/6.0/volume;
		cc[i]*=1.0/6.0/volume;
		dd[i]*=1.0/6.0/volume;
	}
	
	bb[3]=-bb[0]-bb[1]-bb[2];
	cc[3]=-cc[0]-cc[1]-cc[2];
	dd[3]=-dd[0]-dd[1]-dd[2];
	
	
	for (i=0;i<4;i++){
		for (j=0;j<4;j++){
			Me[i][j]=volume/20./kappa;
		}		
	}
	for (i=0;i<4;i++) Me[i][i]*=2;
	
	
	double mod_vv = sqrt(vv[0]*vv[0]+vv[1]*vv[1]+vv[2]*vv[2]);
	double unit_vv[3];
	unit_vv[0]=vv[0]/mod_vv;
	unit_vv[1]=vv[1]/mod_vv;
	unit_vv[2]=vv[2]/mod_vv;
	
	double dh=0.0,min_dh=.0,max_dh=.0,aux_dh;
	
	if (mod_vv>0){
		for (i=1;i<4;i++){
			aux_dh =(x[0]-x[i])*unit_vv[0];
			aux_dh+=(y[0]-y[i])*unit_vv[1];
			aux_dh+=(z[0]-z[i])*unit_vv[2];
			if (min_dh>aux_dh) min_dh=aux_dh;
			if (max_dh<aux_dh) max_dh=aux_dh;
		}
		dh = max_dh-min_dh;
	}
	dh/=4.0;
	
	
	/*double dz = fabs(z14);
	if (fabs(z24)>dz) dz=fabs(z24);
	if (fabs(z34)>dz) dz=fabs(z34);
	dz/=2.0;*/
	
	double Peclet = mod_vv*dh/kappa;
	double alpha_P;
	if (Peclet<0.01) 
		alpha_P = Peclet/3-Peclet*Peclet*Peclet/45;
	else 
		alpha_P = (exp(Peclet)+exp(-Peclet))/(exp(Peclet)-exp(-Peclet))-1.0/Peclet;
	
	
	
	for (i=0;i<4;i++){
		for (j=0;j<4;j++){
			Ke[i][j]=  volume*(bb[j]*vv[0]+cc[j]*vv[1]+dd[j]*vv[2])/4/kappa;
			if (mod_vv>0) Ke[i][j]+= volume*(alpha_P*dh*(bb[i]*vv[0]+cc[i]*vv[1]+dd[i]*vv[2])/mod_vv)*(bb[j]*vv[0]+cc[j]*vv[1]+dd[j]*vv[2])/kappa;
			Ke[i][j]+= volume*(bb[i]*bb[j]+cc[i]*cc[j]+dd[i]*dd[j]);
		}		
	}
	
}