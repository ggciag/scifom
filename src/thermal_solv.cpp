/*
 *  thermal_solv.cpp
 *  thermal1
 *
 *  Created by Victor Sacek on 15/03/10.
 *  Copyright 2010 IAG/USP. All rights reserved.
 *
 */


#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void solv_diag_T(double *xx,double *bb, double *Kd);
void a_times_T(double *xx, double *vetor, double *Kd, double **Kf);
void zera_b_cond_T(double *b_vec);

double *Aloc_vector_real (long n);

extern long nodes_thermal;
extern long **Kthermal_conec;
extern long *Kthermal_posconec;

extern long *cond_thermal;

//extern long verif_first;

extern double tempo;

void thermal_solv(double **Kf, double *Kd, double *b_vec, double *x_vec)
{
    
    long i,j,GL,it=0,itmax=2000;
    double erro;
    double ak,akden,bk,bkden,bknum,bnrm;
    double *p,*r,*z;
    
    GL = nodes_thermal;
    
    p=Aloc_vector_real(GL);
    r=Aloc_vector_real(GL);
    z=Aloc_vector_real(GL);
    
    zera_b_cond_T(b_vec);   
    
    /*if (verif_first==0){
        solv_diag(x_vec,b_vec,Kd);   
        verif_first=1;
    }*/
	
	FILE *F_K;    
    F_K = fopen("F_vec1.txt","w");
    for (i=0;i<nodes_thermal;i++){
        fprintf(F_K,"%10.2f %10.2f\n",x_vec[i],b_vec[i]);
    }
    fclose(F_K);
	
    a_times_T(x_vec,r,Kd,Kf);
    for (j=0;j<GL;j++) {
		r[j]=b_vec[j]-r[j];
	}
	
	for (j=0,bnrm=0.0;j<GL;j++) bnrm += b_vec[j]*b_vec[j];
	bnrm=sqrt(bnrm);
	
	
	printf("bnrm = %f tempo = %f\n",bnrm,tempo);
	
	
	//if (bnrm>.100){
		
        solv_diag_T(z,r,Kd); 
        
        //printf("\n");
        while (it<=itmax){
            it++;  
    		for(j=0,bknum=0.0;j<GL;j++) bknum += z[j]*r[j];
    		if (it==1){
    			for(j=0;j<GL;j++){
    				p[j]=z[j];
    			}
    		}
    		else {
                bk=bknum/bkden;
                for (j=0;j<GL;j++) {
                    p[j]=bk*p[j]+z[j];
                }
            }
            bkden=bknum;
            a_times_T(p,z,Kd,Kf);
            for (akden=0.0,j=0;j<GL;j++) akden += z[j]*p[j];
            ak=bknum/akden;
            for (j=0;j<GL;j++){
				x_vec[j]+=ak*p[j];
				r[j]-=ak*z[j];			
            }
            solv_diag_T(z,r,Kd);
            
            for (j=0,erro=0.0;j<GL;j++) erro += r[j]*r[j];
        	
			if (bnrm>0) erro=sqrt(erro)/bnrm;
			else erro=sqrt(erro);
        	
        	if (erro<1.0E-12) itmax=1;
        	
        	
        }
    /*}
    else {
        for (j=0;j<GL;j++){
            x_vec[j]=0;
        } 
    }*/
    
    printf("\r%f %ld %f      ",tempo,it,erro);
	
    free(p);
    free(r);
    free(z);
    
    
    
    
    
    
    
}

void solv_diag_T(double *xx,double *bb, double *Kd)
{
    long i;
    
    for (i=0;i<nodes_thermal;i++) {        
		if (cond_thermal[i]==0){
			xx[i]=(Kd[i] != 0.0 ? bb[i]/Kd[i] : bb[i]); 
		}	
    }
} 

void a_times_T(double *xx, double *vetor, double *Kd, double **Kf)
{
    long i,pos_aux,cont,j_aux;
    double vetor_aux;
    
    for (i=0;i<nodes_thermal;i++){
        pos_aux=Kthermal_posconec[i];        
		if (cond_thermal[i]==0){
			vetor_aux=xx[i]*Kd[i];
			for (cont=0;cont<pos_aux;cont++){
				j_aux=Kthermal_conec[i][cont];
				vetor_aux+=Kf[i][cont]*xx[j_aux];				
			}            
			vetor[i]=vetor_aux;
		}                
    }
}

void zera_b_cond_T(double *b_vec)
{
    long i;
    for (i=0;i<nodes_thermal;i++) {	
		if (cond_thermal[i]==1){ 
			b_vec[i]=0;
		}	
    }
}
