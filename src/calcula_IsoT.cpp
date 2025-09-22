
extern long Nx;
extern long Ny;
extern long Nz;

extern double depth;

extern double *Temper;

extern double *IsoT;

void calcula_IsoT()
{
    long cond;
    long i,j,k;
    double d1,d2,T2,T1;
    double Tr = 1200.0;
    for (i=0; i<Nx; i++){
        for (j=0; j<Ny; j++){
            for (k=0,cond=0; k<Nz && cond==0;k++){
                if (Temper[i*Ny*Nz + j*Nz + k]>Tr){
                    d1 = double(k-1)*depth/double(Nz-1);
                    d2 = double(k)*depth/double(Nz-1);
                    T1 = Temper[i*Ny*Nz + j*Nz + k-1];
                    T2 = Temper[i*Ny*Nz + j*Nz + k];
                    IsoT[i*Ny+j]=-((T2-T1)*(d2-d1)/(Tr-T1)+d1);
                    cond=1;
                }
            }  
            IsoT[i*Ny+j]=Temper[i*Ny*Nz + j*Nz + Nz-1];       
        }
    }
}
    
    
