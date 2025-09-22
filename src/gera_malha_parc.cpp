
#include <stdio.h>
#include <stdlib.h>


extern double maxx;
extern double minx;
extern double maxy;
extern double miny;

extern double Edge;

double *Aloc_vector_real (long n);
long *Aloc_vector_long (long n);
int incircle(double **xy, long i, long j, long k, long l);
double area(double **xy, long **Tri, long contador);

long gera_malha_parc(long n, long n_old, double **xy, long **Tri,long tri)
{
    long i,j,k,l;
    long ii,t,verif;
    long *tri_prov, contt, contt2,tt;
    long *n_prov,contn;
        
        
        
        
    tri_prov = Aloc_vector_long(400);
    n_prov = Aloc_vector_long(400);
    
    printf("Creating triangles\n");
    t=tri;
    for (ii=n_old;ii<n;ii++){
        
        
        
        printf("\r%ld %ld: %ld",n,n_old,ii);
        contn=0;
        contt=0;
        
        //printf(" 1 ");
        
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
                
        
        //printf(" 2 %d | ",contn);
        
        n_prov[contn]=ii;
        contn++;
        contt2=0;
        for (i=0;i<contn-2;i++){
            for (j=i+1;j<contn-1;j++){
                
                k=contn-1;      
                verif=1;
                for (l=0;l<contn&&verif==1;l++){
                    if (i!=l && j!=l && k!=l){
                        verif=incircle(xy,n_prov[i],n_prov[j],n_prov[k],n_prov[l]);
                    }
                }
                if (verif==1){
                              
                    if (contt2<contt){
                        Tri[tri_prov[contt2]][0]=n_prov[i];
                        Tri[tri_prov[contt2]][1]=n_prov[j];
                        Tri[tri_prov[contt2]][2]=n_prov[k];
                        //printf("a%d|",tri_prov[contt2]);
                        contt2++;
                    }
                    else {
                        Tri[t][0]=n_prov[i];
                        Tri[t][1]=n_prov[j];
                        Tri[t][2]=n_prov[k];
                        //printf("n%d|",t);
                        t++;
                    } 
                }
            }
        }
        /*printf(" 3 %7.0f %7.0f | %7.0f %7.0f | %7.0f %7.0f",
                 xy[Tri[t-1][0]][0],xy[Tri[t-1][0]][1],
                 xy[Tri[t-1][1]][0],xy[Tri[t-1][1]][1],
                 xy[Tri[t-1][2]][0],xy[Tri[t-1][2]][1]);*/
    }
    
    free(tri_prov);
    free(n_prov);
    
    
    
    double area_parc;
    
    for (i=0;i<t;i++){
        area_parc=area(xy, Tri, i);
		if (Tri[i][0]==5165 || Tri[i][1]==5165 || Tri[i][2]==5165){
			printf("\n\nAREA%lf\n\n", area_parc);
			printf("%ld %lf %lf\n",Tri[i][0],xy[Tri[i][0]][0],xy[Tri[i][0]][1]);
			printf("%ld %lf %lf\n",Tri[i][1],xy[Tri[i][1]][0],xy[Tri[i][1]][1]);
			printf("%ld %lf %lf\n",Tri[i][2],xy[Tri[i][2]][0],xy[Tri[i][2]][1]);
		}
        if (area_parc<0.0001*Edge*Edge) {
           printf("Menor %ld %lf\n\n",i,area_parc);
           if (i<t-1){
               Tri[i][0]=Tri[t-1][0];
               Tri[i][1]=Tri[t-1][1];
               Tri[i][2]=Tri[t-1][2];
			   i--;
           }
           t--;
        }
    }
	
    
    return(t);
} 
