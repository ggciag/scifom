#include <stdio.h>
#include <stdlib.h>

double *Aloc_vector_real (long n);
long *Aloc_vector_long (long n);
int incircle(double **xy, long i, long j, long k, long l);
double area(double **xy, long **Tri, long contador);

extern double Edge;

long gera_malha(long n, double **xy, long **Tri)
{
    long i,j,k,l;
    long ii,t,verif;
    long *tri_prov, contt, contt2,tt;
    long *n_prov,contn;
    
    Tri[0][0]=0; Tri[0][1]=1;Tri[0][2]=3;
    Tri[1][0]=1; Tri[1][1]=2;Tri[1][2]=3;
    
    tri_prov = Aloc_vector_long(400);
    n_prov = Aloc_vector_long(400);
    
    printf("Creating triangles\n");
    t=2;
    for (ii=4;ii<n;ii++){
        
        
        
        printf("\r%ld",ii);
        contn=0;
        contt=0;
        
        for (tt=0;tt<t;tt++){
            if (incircle(xy,Tri[tt][0],Tri[tt][1],Tri[tt][2],ii)==0){
                tri_prov[contt]=tt;
                contt++;
                for (i=0;i<3;i++){
                    verif=0;
                    for (j=0;j<contn;j++){
                        if (Tri[tt][i]==n_prov[j]){
                            verif=1;
                        }
                    } 
                    if (verif==0){
                        n_prov[contn]=Tri[tt][i];
                        contn++;
                    }
                }
            }
        }
        //printf("teste");
        //printf(" %d %d %d",contn,contt,t);            
        
        n_prov[contn]=ii;
        contn++;
        contt2=0;
        for (i=0;i<contn-2;i++){
            for (j=i+1;j<contn-1;j++){
                //for (k=j+1;k<contn;k++){
                k=contn-1;      
                verif=1;
                for (l=0;l<contn&&verif==1;l++){
                    if (i!=l && j!=l && k!=l){
                        verif=incircle(xy,n_prov[i],n_prov[j],n_prov[k],n_prov[l]);
                    }
                }
                if (verif==1){
                    //printf("+");
                    if (contt2<contt){
                        Tri[tri_prov[contt2]][0]=n_prov[i];
                        Tri[tri_prov[contt2]][1]=n_prov[j];
                        Tri[tri_prov[contt2]][2]=n_prov[k];
                        contt2++;
                    }
                    else {
                        Tri[t][0]=n_prov[i];
                        Tri[t][1]=n_prov[j];
                        Tri[t][2]=n_prov[k];
                        t++;
                    } 
                }
                //}
            }
        }
       
        //printf(" teste2");
    }
    
    free(tri_prov);
    free(n_prov);
    
    double area_parc;
    
    for (i=0;i<t;i++){
        area_parc=area(xy, Tri, i);           
        if (area_parc<0.0001*Edge*Edge) {
           printf("Menor %ld\n\n",i);
           if (i<t-1){
               Tri[i][0]=Tri[t-1][0];
               Tri[i][1]=Tri[t-1][1];
               Tri[i][2]=Tri[t-1][2];                          
           }
           t--;
        }
        
    }
    
    return(t);
} 

double area(double **xy, long **Tri, long contador)
{
    double delta; 
    double xi,xj,xk,yi,yj,yk;  
    
    xi=xy[Tri[contador][0]][0]; yi=xy[Tri[contador][0]][1];
	xj=xy[Tri[contador][1]][0]; yj=xy[Tri[contador][1]][1];
	xk=xy[Tri[contador][2]][0]; yk=xy[Tri[contador][2]][1];
    
    delta=0.5*((xj*yk-yj*xk)+xi*(yj-yk)+yi*(xk-xj));   
    if (delta<0) delta*=-1.0;
    return(delta);
}

long importa_malha(long **Tri){
    
	FILE *f_malha;

    long n_tri;
	
	f_malha = fopen("malha_externa.txt","r");

    if (f_malha==NULL){
        printf("File malha_externa.txt does not exist\n");
        exit(-1);
    }

    fscanf(f_malha,"%ld",&n_tri);
	
	for (long i=0;i<n_tri;i++){
		fscanf(f_malha, "%ld %ld %ld",&Tri[i][0],&Tri[i][1],&Tri[i][2]);
	}
	
	fclose(f_malha);

    return (n_tri);
}