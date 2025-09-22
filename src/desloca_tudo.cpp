#include <stdio.h>
#include <stdlib.h>

#include <math.h>

extern double maxx;
extern double minx;
extern double maxy;
extern double miny;

extern double Edge;


int intriangle(double **pontos, long i, long j, long k, double xx, double yy);


extern double *area_vor;
extern double **aresta_vor;
extern double **dist_vor;

extern long **conec;
extern long *pos_conec;
extern long **conec_Tri;
extern long *pos_conec_Tri;

extern long max_conec_p;

extern long tri_flex;
extern long **Tri_flex;
extern long nodes_flex;
extern double **xy_flex;

extern double **peso;
extern long **peso_pos;

extern double *h_topo;
extern double *h_bed;

extern double dx_flex;
extern double dy_flex;

extern double minx_flex;
extern double miny_flex;

extern double pos_falha;
extern double pos_falha_ant;

extern long Nx;
extern long Ny;
extern long Nz;


void desloca_tudo(long n, double **xy, long tri, long **Tri)
{
    
    long ii,i;//,j,k,l,m,t;
    //long verif;

    
    
    //Pesos
    //printf("pesos\n");
    
    double a1,a2,a3,b1,b2,b3,c1,c2,c3;
    double x,y;
    double xi,xj,xk,yi,yj,yk,delta;
    
    long i_flex,j_flex;
    
    long in_tri;
    
    for (i=0;i<n;i++){
        x = xy[i][0];
        y = xy[i][1];
        if (y>maxy-Edge*0.0001) y-=0.0001*Edge;
        if (y<miny+Edge*0.0001) y+=0.0001*Edge;
        
        i_flex = long(floor((x-minx_flex)/dx_flex));
        j_flex = long(floor((y-miny_flex)/dy_flex));
		
        ii=2*(i_flex*(Ny-1)+j_flex);
        in_tri=0;
        if (intriangle(xy_flex,Tri_flex[ii][0],Tri_flex[ii][1],Tri_flex[ii][2],x,y)==1){
            ii=2*(i_flex*(Ny-1)+j_flex)+1;
            if (intriangle(xy_flex,Tri_flex[ii][0],Tri_flex[ii][1],Tri_flex[ii][2],
						   x,y)==0){
            }
            else {
                /*printf("Epa!\n");
				 printf("%d %d ",i_flex,j_flex); 
				 printf("| %d %d| ",ii,tri_flex);
				 printf("%f %f | \n%f %f\n",x,y,xy_flex[Tri_flex[ii][0]][0],xy_flex[Tri_flex[ii][0]][1]);
				 printf("%f %f\n\n",xy_flex[Tri_flex[ii][0]][0]+dx_flex,xy_flex[Tri_flex[ii][0]][1]+dy_flex);*/
                in_tri=1;
            }                                                       
        }
        
        
        peso_pos[i][0]=Tri_flex[ii][0];
        peso_pos[i][1]=Tri_flex[ii][1];
        peso_pos[i][2]=Tri_flex[ii][2];
        
        xi = xy_flex[Tri_flex[ii][0]][0]; 
        yi = xy_flex[Tri_flex[ii][0]][1];
	    xj = xy_flex[Tri_flex[ii][1]][0]; 
        yj = xy_flex[Tri_flex[ii][1]][1];
	    xk = xy_flex[Tri_flex[ii][2]][0]; 
        yk = xy_flex[Tri_flex[ii][2]][1];
        
        delta=0.5*((xj*yk-yj*xk)+xi*(yj-yk)+yi*(xk-xj));
        
        a1=xj*yk-xk*yj; b1=yj-yk; c1=xk-xj;
		a2=xk*yi-xi*yk; b2=yk-yi; c2=xi-xk;
		a3=xi*yj-xj*yi; b3=yi-yj; c3=xj-xi;
		
		peso[i][0] = (a1+b1*x+c1*y)/2/delta;
		peso[i][1] = (a2+b2*x+c2*y)/2/delta;
		peso[i][2] = (a3+b3*x+c3*y)/2/delta;
		
		if (in_tri==1) {
            printf("%f %f %f\n\n",peso[i][0],peso[i][1],peso[i][2]);
            
        }
		
		
        
    }
    
} 


