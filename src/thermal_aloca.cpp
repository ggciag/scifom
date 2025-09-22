/*
 *  thermal_aloca.cpp
 *  thermal1
 *
 *  Created by Victor Sacek on 12/03/10.
 *  Copyright 2010 IAG/USP. All rights reserved.
 *
 */
#include <stdio.h>
#include <math.h>

double **Aloc_matrix_real (long p, long n);
long **Aloc_matrix_long (long p, long n);
double *Aloc_vector_real (long n);
long *Aloc_vector_long (long n);

void lateral(double **xyz, long ntet, long **tetra);

extern long tetra_thermal;
extern long **Tet_thermal;
extern long nodes_thermal;
extern double **xyz_thermal;
extern double **xyz_thermal_fut;

extern long tri_flex;
extern long nodes_flex;

extern long **Tri_flex;
extern double **xy_flex;
extern double *h_top;
extern double *v_adv2D;

extern double h_bot;

extern double **v_adv;

extern double **Kthermal;
extern double *Kthermal_diag;
extern double **Mthermal;
extern double *Mthermal_diag;
extern double **Kthermal1;
extern double *Kthermal1_diag;
extern double **Kthermal2;
extern double *Kthermal2_diag;

extern double *fthermal;

extern long *cond_thermal;

extern double *T_vec;
extern double *T_vec_fut;
extern double *T_vec0;

extern double **TKe;
extern double **TMe;

extern long **Kthermal_conec;
extern long *Kthermal_posconec;

extern long thermal_max_conec_p;

extern double deltaz;

extern long layers;


extern long **tri_lat;
extern long ntri_lat;



void thermal_aloca()
{
	
	
	printf("e");
	
	nodes_thermal=nodes_flex*layers;
	tetra_thermal=tri_flex*(layers-1)*3;
	
	
	long i,j,k,l,m,nn,t,cont,verif;
	
	
	printf("d");
	
	h_top = Aloc_vector_real(nodes_flex);
	v_adv2D = Aloc_vector_real(nodes_flex);
	
	TKe = Aloc_matrix_real(4, 4);
	TMe = Aloc_matrix_real(4, 4);
	
	Tet_thermal=Aloc_matrix_long(tetra_thermal, 4);
	xyz_thermal=Aloc_matrix_real(nodes_thermal, 3);
	xyz_thermal_fut=Aloc_matrix_real(nodes_thermal, 3);
	
	v_adv = Aloc_matrix_real(nodes_thermal, 3);
	
	printf("c");
	
	Kthermal_conec = Aloc_matrix_long (nodes_thermal,thermal_max_conec_p);
	Kthermal_posconec = Aloc_vector_long (nodes_thermal);
	
	printf("b");
	
	Kthermal = Aloc_matrix_real (nodes_thermal,thermal_max_conec_p);
	Kthermal_diag = Aloc_vector_real(nodes_thermal);
	Mthermal = Aloc_matrix_real (nodes_thermal,thermal_max_conec_p);
	Mthermal_diag = Aloc_vector_real(nodes_thermal);
	Kthermal1 = Aloc_matrix_real (nodes_thermal,thermal_max_conec_p);
	Kthermal1_diag = Aloc_vector_real(nodes_thermal);
	Kthermal2 = Aloc_matrix_real (nodes_thermal,thermal_max_conec_p);
	Kthermal2_diag = Aloc_vector_real(nodes_thermal);
	
	fthermal = Aloc_vector_real(nodes_thermal);
	T_vec = Aloc_vector_real(nodes_thermal);
	T_vec_fut = Aloc_vector_real(nodes_thermal);
	T_vec0 = Aloc_vector_real(nodes_thermal);
	
	cond_thermal = Aloc_vector_long(nodes_thermal);
	
	
	printf("b");
	
	for (j=0;j<layers;j++){
		for (i=0;i<nodes_flex;i++){
			xyz_thermal[j*nodes_flex+i][0]=xy_flex[i][0];
			xyz_thermal[j*nodes_flex+i][1]=xy_flex[i][1];
			xyz_thermal_fut[j*nodes_flex+i][0]=xy_flex[i][0];
			xyz_thermal_fut[j*nodes_flex+i][1]=xy_flex[i][1];
			xyz_thermal[j*nodes_flex+i][2]=h_top[i] + j*(h_bot-h_top[i])/(layers-1);
			if (j==0 || j==layers-1) cond_thermal[j*nodes_flex+i]=1;
		}
	}
	
	printf("a %ld %ld",Tri_flex[0][0],Tri_flex[tri_flex-1][2]);
	
	long pc0,pc1,pc2,pc4,pc5,pc6;
	
	for (j=0,cont=0;j<layers-1;j++){
		for (t=0;t<tri_flex;t++){
			pc0 = j*nodes_flex+Tri_flex[t][0];
			pc1 = j*nodes_flex+Tri_flex[t][1];
			pc2 = j*nodes_flex+Tri_flex[t][2];
			
			pc4 = (j+1)*nodes_flex+Tri_flex[t][0];
			pc5 = (j+1)*nodes_flex+Tri_flex[t][1];
			pc6 = (j+1)*nodes_flex+Tri_flex[t][2];
			
			Tet_thermal[cont][0]=pc0; Tet_thermal[cont][1]=pc4;	Tet_thermal[cont][2]=pc5; Tet_thermal[cont][3]=pc6; cont++;
			Tet_thermal[cont][0]=pc0; Tet_thermal[cont][1]=pc1;	Tet_thermal[cont][2]=pc5; Tet_thermal[cont][3]=pc6; cont++;
			Tet_thermal[cont][0]=pc0; Tet_thermal[cont][1]=pc1;	Tet_thermal[cont][2]=pc2; Tet_thermal[cont][3]=pc6; cont++;
			
		}
	}
	
	FILE *F_tri;
	
	F_tri = fopen("F_tri.txt", "w");
	for (t=0;t<tetra_thermal;t++){
		fprintf(F_tri, "%ld %ld %ld %ld\n", Tet_thermal[t][0], Tet_thermal[t][1], Tet_thermal[t][2], Tet_thermal[t][3] );
		//fprintf(F_tri, "%ld\n",Tet_thermal[
	}
	fclose(F_tri);
	
	for (t=0;t<nodes_thermal;t++){
        Kthermal_posconec[t]=1;
        Kthermal_conec[t][0]=t;
    }
	
	for (t=0;t<tetra_thermal;t++){        
        for (i=0;i<4;i++){
            if (i==0) {k=Tet_thermal[t][0]; l=Tet_thermal[t][1]; m=Tet_thermal[t][2]; nn=Tet_thermal[t][3];}
            if (i==1) {k=Tet_thermal[t][1]; l=Tet_thermal[t][2]; m=Tet_thermal[t][3]; nn=Tet_thermal[t][0];}
            if (i==2) {k=Tet_thermal[t][2]; l=Tet_thermal[t][3]; m=Tet_thermal[t][0]; nn=Tet_thermal[t][1];}            
            if (i==3) {k=Tet_thermal[t][3]; l=Tet_thermal[t][0]; m=Tet_thermal[t][1]; nn=Tet_thermal[t][2];}
            for (verif=0,j=0;j<Kthermal_posconec[k];j++){
                if (Kthermal_conec[k][j]==l) verif=1;
            }
            if (verif==0){
                Kthermal_conec[k][Kthermal_posconec[k]]=l;
                Kthermal_posconec[k]++;
                if (Kthermal_posconec[k]>=thermal_max_conec_p){
                    printf("perigo!");
					
                }
            }
            for (verif=0,j=0;j<Kthermal_posconec[k];j++){
                if (Kthermal_conec[k][j]==m) verif=1;
            }
            if (verif==0){
                Kthermal_conec[k][Kthermal_posconec[k]]=m;
                Kthermal_posconec[k]++;
                if (Kthermal_posconec[k]>=thermal_max_conec_p){
                    printf("perigo!");
					
                }
            }
			for (verif=0,j=0;j<Kthermal_posconec[k];j++){
                if (Kthermal_conec[k][j]==nn) verif=1;
            }
            if (verif==0){
                Kthermal_conec[k][Kthermal_posconec[k]]=nn;
                Kthermal_posconec[k]++;
                if (Kthermal_posconec[k]>=thermal_max_conec_p){
                    printf("perigo!");
					
                }
            }
			
        }
    }
	
	FILE *F_conec;
    
    F_conec = fopen("Flex_conec.txt","w");
    for (t=0;t<nodes_thermal;t++){   
        fprintf(F_conec,"%ld | %ld |",Kthermal_posconec[t],cond_thermal[t]);
		fprintf(F_conec,"%10.2f %10.2f %10.2f|",xyz_thermal[t][0],xyz_thermal[t][1],xyz_thermal[t][2]);
        for (i=0;i<Kthermal_posconec[t];i++){
            fprintf(F_conec,"%ld ",Kthermal_conec[t][i]);
        }
        fprintf(F_conec,"\n");
    }
    fclose(F_conec);
	
	
}


