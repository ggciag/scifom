#include <stdio.h>

extern double maxx;
extern double minx;

extern double pos_falha;
extern double pos_falha_ant;

extern double inclina;
double var_fault(double v);

extern double **pos_falha_vec;
extern long cont_falha;

extern long *side;

void desloca(long n, double **xy,double pos_falha,double despl,long tri,long **Tri)
{
    long i,t,cont;
    long k,l,m;
	if (pos_falha<maxx) maxx+=despl;
	if (pos_falha<minx) minx+=despl;

    for (i=0;i<n;i++){
        //if (xy[i][0]>pos_falha+var_fault(xy[i][1])){
		if (side[i]*pos_falha_vec[cont_falha][1]>0){
			xy[i][0]+=despl;
        }
        xy[i][2]=0;
    }
    for (t=0; t<tri; t++){
        Tri[t][3]=0;
        cont=0;
        k=Tri[t][0];
        l=Tri[t][1];
        m=Tri[t][2];
        //if (xy[k][0]>pos_falha+var_fault(xy[k][1])) cont++;
        //if (xy[l][0]>pos_falha+var_fault(xy[l][1])) cont++;
        //if (xy[m][0]>pos_falha+var_fault(xy[m][1])) cont++;
		
		if (side[k]*pos_falha_vec[cont_falha][1]>0) cont++;
		if (side[l]*pos_falha_vec[cont_falha][1]>0) cont++;
		if (side[m]*pos_falha_vec[cont_falha][1]>0) cont++;		
	
        if (cont>0 && cont<3){
            //Tri[t][3]=1;
            xy[Tri[t][0]][2]=1;
            xy[Tri[t][1]][2]=1;
            xy[Tri[t][2]][2]=1;
        }
    }
    for (t=0; t<tri; t++){
        if (xy[Tri[t][0]][2]==1 || xy[Tri[t][1]][2]==1 || xy[Tri[t][2]][2]==1){
            Tri[t][3]=1;
        }
    }
    for (t=0; t<tri; t++){
        if (Tri[t][3]==1){
            xy[Tri[t][0]][2]=1;
            xy[Tri[t][1]][2]=1;
            xy[Tri[t][2]][2]=1;
        }
    }
    
    
    
    
    cont=0;
    for (t=0; t<tri; t++){
        if (Tri[t][3]==1) cont++;
    }
    
    printf("num tri: %ld %ld\nnodes %ld",tri,cont,n);    
    
    
}
