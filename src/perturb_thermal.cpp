
#include <math.h>

extern long Nx;
extern long Ny;
extern long Nz;

extern double *Temper;
extern double **xyz_thermal;
extern double depth;

void perturb_thermal(double tempo){
    
    long i;
    
    for (i=0;i<Nx*Ny*Nz;i++){
        if (xyz_thermal[i][2]<-depth+5.0){
            if (xyz_thermal[i][0]<400000.0){
                if (tempo<0.5E6){
                   Temper[i]=1300+200*(tempo/0.5E6);
                }
                if (tempo>=0.5E6 && tempo<4.5E6){
                   Temper[i]=1500;
                }
                if (tempo>=4.5E6 && tempo<=5.0E6){
                   Temper[i]=1300+200*(5.0E6-tempo)/0.5E6;
                }
            }
        }
    }
    
    
    /*long i,j,k,ii,jj;
    ii=Nx/4;
    jj=Ny/2;
    k=Nz-1;
    long marg=10;
    double id,jd;
    
    for (i=ii-marg; i<=ii+marg;i++){
        for (j=jj-marg; j<=jj+marg;j++){
            id=double(i)-ii; jd=double(j)-jj;
            if (id*id+jd*jd<25){
                if (tempo<1E5) Temper[i*Ny*Nz + j*Nz + k]=1300.0+1000.0*tempo/1E5;
                if (tempo>=3E6) {
                    if (tempo<=3E6+1E5) Temper[i*Ny*Nz + j*Nz + k]=2300.0-1000.0*(tempo-3E6)/1E5;
                }
            }
        }
    }*/
}
