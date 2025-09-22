#include <stdlib.h>
#include <stdio.h>


extern long nodes;
extern long **conec;
extern long *pos_conec;

extern double *area_vor;
extern double **aresta_vor;
extern double **dist_vor;

extern long *ordem_fluvi;
extern long *direc_fluvi;
extern double *dist_fluvi;

extern double *Qr;
extern double *Qf;

extern double *Df;


extern double *h_topo;
extern double *h_bed;

extern double dt;

extern double nivel;

extern double vR;
extern double Kf;
extern double ls;
extern double lb;

void fluvi()
{
    long i,j,jj;
    double dist_aux;
    long cont_fluvi=nodes;
    double Qeqb;

    for (i=0;i<nodes;i++){  
        direc_fluvi[i]=i;
        Qr[i]=dt*area_vor[i]*vR;   
        Qf[i]=0;   
        Df[i]=0;  
        ordem_fluvi[i]=1;
        if (h_topo[i]>nivel){      
            jj=i; 
            dist_aux=100; //qualquer numero;
            for (j=0;j<pos_conec[i];j++){
                if (h_topo[jj]>h_topo[conec[i][j]]){
                    if (aresta_vor[i][j]>0){
                        jj=conec[i][j];
                        dist_aux=dist_vor[i][j];
                    }
                }
            }
            direc_fluvi[i]=jj;
            dist_fluvi[i]=dist_aux;
        }       
    }       
    
    for (i=0;i<nodes;i++){        
        j = direc_fluvi[i];
        if (j!=i){
            ordem_fluvi[j]=0;                
        }        
    }
    
    while (cont_fluvi>0){
        
        for (i=0;i<nodes;i++){
            if (ordem_fluvi[i]==1){
                j=direc_fluvi[i];                   
                if (j!=i){
                    if (h_topo[j]>nivel)
                        Qeqb = Kf*Qr[i]*(h_topo[i]-h_topo[j])/dist_fluvi[i];
                    else
                        Qeqb = Kf*Qr[i]*(h_topo[i]-nivel)/dist_fluvi[i];
                    if (Qeqb<=Qf[i]){
                        Df[i]=(Qf[i]-Qeqb)/area_vor[i];
                        Qf[j]+=Qeqb;
                    }
                    else {
                        if (h_bed[i]<h_topo[i]) {
                            Df[i]=((Qf[i]-Qeqb)/area_vor[i])*(dist_fluvi[i]/ls);
                            Qf[j]+=Qf[i]+(Qeqb-Qf[i])*(dist_fluvi[i]/ls);   							
                        }                                
                        else {
                            Df[i]=((Qf[i]-Qeqb)/area_vor[i])*(dist_fluvi[i]/lb);
                            Qf[j]+=Qf[i]+(Qeqb-Qf[i])*(dist_fluvi[i]/lb);
                        }  
                    }                  
                    Qr[j]+=Qr[i];
                }
                else {
                    Df[i]=Qf[i]/area_vor[i]; 
                }
                ordem_fluvi[i]=2;
                cont_fluvi--;
                
                //printf("%f\n",Df[i]);
            }
        }
        
        for (i=0;i<nodes;i++){
            if (ordem_fluvi[i]==0){
                ordem_fluvi[i]=1;
            }
        }
        
        for (i=0;i<nodes;i++){
            j = direc_fluvi[i];
            if (ordem_fluvi[i]!=2 && j!=i){
                ordem_fluvi[j]=0;
            }
        }
        
        
    }
                
                 
    
}
    
                                    
           
         
