#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void solv_diag(double *xx,double *bb, double *Kd);
void a_times(double *xx, double *vetor, double *Kd, double **Kf);
void zera_b_cond(double *bflex);

double *Aloc_vector_real (long n);

extern long nodes_flex;
extern long **Kconec;
extern long *Kposconec;
//extern double *wflex;
extern double *wflex_aux;
extern double *wflex_cumula;
extern long **cond_c;

extern long verif_first;

extern double tempo;

void flexura_solv(double **Kf, double *Kd, double *bflex, double *wflex)
{
    
    long i,j,GL,it=0,itmax=2000;
    double erro;
    double ak,akden,bk,bkden,bknum,bnrm;
    double *p,*r,*z;
    
    GL = 3*nodes_flex;
    
    p=Aloc_vector_real(GL);
    r=Aloc_vector_real(GL);
    z=Aloc_vector_real(GL);
    
    zera_b_cond(bflex);   
    
    if (verif_first==0){
        solv_diag(wflex,bflex,Kd);   
        verif_first=1;
    }
    a_times(wflex,r,Kd,Kf);
    for (j=0;j<GL;j++) {
		r[j]=bflex[j]-r[j];
	}
		
	for (j=0,bnrm=0.0;j<GL;j++) bnrm += bflex[j]*bflex[j];
	bnrm=sqrt(bnrm);
	
	
	printf("bnrm = %f tempo = %f \n",bnrm,tempo);
	
	
	if (bnrm>100){
	
        solv_diag(z,r,Kd); 
        
        printf("\n");
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
            a_times(p,z,Kd,Kf);
            for (akden=0.0,j=0;j<GL;j++) akden += z[j]*p[j];
            ak=bknum/akden;
            for (j=0;j<GL;j++){
                wflex[j]+=ak*p[j];
                r[j]-=ak*z[j];
            }
            solv_diag(z,r,Kd);
            
            for (j=0,erro=0.0;j<GL;j++) erro += r[j]*r[j];
        	erro=sqrt(erro)/bnrm;
        	
        	if (erro<1.0E-2) itmax=1;
        	printf("\r%ld %f      ",it,erro);
        	
        }
    }
    else {
        for (j=0;j<GL;j++){
            wflex[j]=0;
        } 
    }
         
    free(p);
    free(r);
    free(z);
    
    
    for (i=0;i<nodes_flex;i++){
        wflex_aux[i]=wflex[i*3];
        //wflex_cumula[i]+=wflex[i*3];
    }
    
    
    
    
    
}

void solv_diag(double *xx,double *bb, double *Kd)
{
    long i,c1,i_aux;
    
    for (i=0;i<nodes_flex;i++) {
        for (c1=0;c1<3;c1++){
            if (cond_c[i][c1]==0){
                i_aux = i*3+c1;
                xx[i_aux]=(Kd[i_aux] != 0.0 ? bb[i_aux]/Kd[i_aux] : bb[i_aux]); 
            }
        }
    }
} 

void a_times(double *xx, double *vetor, double *Kd, double **Kf)
{
    long i,c1,pos_aux,cont,i_aux,j_aux,contt;
    double vetor_aux;
    
    for (i=0;i<nodes_flex;i++){
        pos_aux=Kposconec[i];
        for (c1=0;c1<3;c1++){ 
            if (cond_c[i][c1]==0){
                i_aux=i*3+c1;
                vetor_aux=xx[i_aux]*Kd[i_aux];
                for (cont=0;cont<pos_aux;cont++){
                    j_aux=Kconec[i][cont]*3;
                    contt=cont*3;
                    vetor_aux+=Kf[i_aux][contt]*xx[j_aux];contt++;j_aux++;
                    vetor_aux+=Kf[i_aux][contt]*xx[j_aux];contt++;j_aux++;
                    vetor_aux+=Kf[i_aux][contt]*xx[j_aux];                    
                    /*for (c2=0;c2<3;c2++){                                               
                        vetor_aux+=Kf[i_aux][cont*3+c2]*xx[j_aux+c2]; 
                    }*/
                }            
                vetor[i_aux]=vetor_aux;
            }
        }        
    }
}

void zera_b_cond(double *bflex)
{
    long i,c1;
    for (i=0;i<nodes_flex;i++) {
        for (c1=0;c1<3;c1++){
            if (cond_c[i][c1]==1){ 
               bflex[i*3+c1]=0;
            }
        }
    }
}
