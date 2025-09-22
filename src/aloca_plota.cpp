#include <stdio.h>
#include <stdlib.h>
#include <stdio.h>

#define NSTACK 50
#define M 7
#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;
#define SWAPint(a,b) inttemp=(a);(a)=(b);(b)=inttemp;

double *Aloc_vector_real (long n);
long *Aloc_vector_long (long n);
void sort(unsigned long n, double arr[], long *Ox);

long *aloca_plota(long tri,double **xy,long **Tri)
{
    long t,i;
    double *dist_plota;
    long *ordem_plota;
    printf("Aloca_plota: ");
    dist_plota = Aloc_vector_real (tri+1);  
    ordem_plota = Aloc_vector_long (tri+1);
    
    
    for (t=0;t<tri;t++){
        ordem_plota[t+1]=t;
        for (i=0;i<3;i++){
            dist_plota[t+1]+=xy[Tri[t][i]][0]*xy[Tri[t][i]][0];
            dist_plota[t+1]+=xy[Tri[t][i]][1]*xy[Tri[t][i]][1];
        }
    }
    sort(tri,dist_plota,ordem_plota);
    
    for (t=0;t<tri;t++){
        dist_plota[t]=dist_plota[t+1];
        ordem_plota[t]=ordem_plota[t+1];
    }   
    
    printf("Fim\n");
    free(dist_plota);
    return(ordem_plota);
    //printf("testePlota");
        
}
void sort(unsigned long n, double arr[], long *Ox)
{
	unsigned long i,ir=n,j,k,l=1;
	long *istack;
	int jstack=0;
	double a,temp;
	long xa,inttemp;

	istack=Aloc_vector_long(NSTACK+1);
	for (;;) {
		if (ir-l < M) {
			for (j=l+1;j<=ir;j++) {
				a=arr[j];
				xa=Ox[j]; 
				for (i=j-1;i>=l;i--) {
					if (arr[i] <= a) break;
					arr[i+1]=arr[i];
					Ox[i+1]=Ox[i]; 
				}
				arr[i+1]=a;
				Ox[i+1]=xa; 
			}
			if (jstack == 0) break;
			ir=istack[jstack--];
			l=istack[jstack--];
		} 
		else {
			k=(l+ir) >> 1;
			SWAP(arr[k],arr[l+1])
			SWAPint(Ox[k],Ox[l+1])
			if (arr[l] > arr[ir]) {
				SWAP(arr[l],arr[ir])
				SWAPint(Ox[l],Ox[ir])
			}
			if (arr[l+1] > arr[ir]) {
				SWAP(arr[l+1],arr[ir])
				SWAPint(Ox[l+1],Ox[ir])
			}
			if (arr[l] > arr[l+1]) {
				SWAP(arr[l],arr[l+1])
				SWAPint(Ox[l],Ox[l+1])
			}
			i=l+1;
			j=ir;
			a=arr[l+1];
			xa=Ox[l+1]; 
			for (;;) {
				do i++; while (arr[i] < a);
				do j--; while (arr[j] > a);
				if (j < i) break;
				SWAP(arr[i],arr[j]);
				SWAPint(Ox[i],Ox[j]);
			}
			arr[l+1]=arr[j];
			arr[j]=a;
			Ox[l+1]=Ox[j]; 
			Ox[j]=xa;      			
			jstack += 2;
			if (jstack > NSTACK) {
				printf("NSTACK too small in sort.");
				exit(-1);
			}
			if (ir-i+1 >= j-l) {
				istack[jstack]=ir;
				istack[jstack-1]=i;
				ir=j-1;
			}
			else {
				istack[jstack]=j-1;
				istack[jstack-1]=l;
				l=i;
			}
		}
	}
	free(istack);
}
