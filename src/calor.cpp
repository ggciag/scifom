#include <stdio.h>

extern long Nx;
extern long Ny;
extern long Nz;

extern double *Temper;

extern long **Kthermal_conec;
extern long *Kthermal_posconec;

extern double **xyz_thermal;
extern double **v_adv;

extern long *cond_borda_Temper;
extern long *cond_borda_Temper2;

extern double *Dtemper;

extern long nodes_thermal;

extern double kappa;
extern double seg_per_ano;

extern double tempo;
extern double tempo_muda_falha;

void calor(double dt_calor){
     
    //double Tvar;
    long ti,tj,j;
    //double vv[3], nn[3];
    
    for (ti=0;ti<nodes_thermal;ti++){ 
        Dtemper[ti]=0;
    }
    for (ti=0;ti<nodes_thermal;ti++){ 
        for (j=0;j<Kthermal_posconec[ti];j++){
            tj = Kthermal_conec[ti][j];
            /*
            //Difusao
            Tvar = Temper[ti]-Temper[tj];
            if (Tvar>0){
                Dtemper[ti]-=kappa*Area_Temper[ti][j]*(Tvar/Dist_Temper[ti][j]);
                Dtemper[tj]+=kappa*Area_Temper[ti][j]*(Tvar/Dist_Temper[ti][j]);
            }
            
            
            //Advecção
            if (tempo<=3000000.0){
                nn[0]=(xyz_thermal[tj][0]-xyz_thermal[ti][0])/Dist_Temper[ti][j];
                nn[1]=(xyz_thermal[tj][1]-xyz_thermal[ti][1])/Dist_Temper[ti][j];
                nn[2]=(xyz_thermal[tj][2]-xyz_thermal[ti][2])/Dist_Temper[ti][j];
                
                vv[0]= (v_adv[ti][0]*Temper[ti]+v_adv[tj][0]*Temper[tj])/2;
                vv[1]= (v_adv[ti][1]*Temper[ti]+v_adv[tj][1]*Temper[tj])/2;            
                vv[2]= (v_adv[ti][2]*Temper[ti]+v_adv[tj][2]*Temper[tj])/2;
                
                Dtemper[ti]-=(nn[0]*vv[0]+nn[1]*vv[1]+nn[2]*vv[2])*Area_Temper[ti][j];
            }
			 */
            
        }
    }
    /*for (ti=0;ti<nodes_thermal;ti++){        
        Temper[ti]+=cond_borda_Temper[ti]*Dtemper[ti]*dt_calor*seg_per_ano/Volume_Temper[ti];
    }*/
    
    for (ti=0;ti<nodes_thermal;ti++){  
        if (cond_borda_Temper2[ti]>-1){
           Temper[ti]=Temper[cond_borda_Temper2[ti]];
        }
    }
    
    //long i,k;
    /*for (i=0; i<Nx; i++){
        for (j=0; j<Ny; j++){
            for (k=0; k<Nz;k++){
                printf("%d %d %d %f\n",i,j,k,Dtemper[i*Ny*Nz+j*Nz+k]/Volume_Temper[i*Ny*Nz+j*Nz+k]);
                
            }
        }
    }*/
    /*i=Nx/2;j=Ny/2;
    for (k=0; k<Nz;k++){
        printf("%d %d %d %f\n",i,j,k,Temper[i*Ny*Nz+j*Nz+k]);
        
    }*/
}
                
    
