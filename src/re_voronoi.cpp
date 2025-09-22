#include <stdio.h>
#include <stdlib.h>

#include <math.h>

extern double maxx;
extern double minx;
extern double maxy;
extern double miny;

extern double Edge;

long *Aloc_vector_long (long n);
long **Aloc_matrix_long (long p, long n);

double *Aloc_vector_real (long n);

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


void re_voronoi(long n, double **xy, long tri, long **Tri)
{
    
    long i,j,k,l,m,t;
    long verif;
    
    for (i=0;i<n;i++){
        if (xy[i][2]==1){
             pos_conec[i]=0;  
             pos_conec_Tri[i]=0;          
        }
    }   
    
    for (t=0;t<tri;t++){        
        for (i=0;i<3;i++){
            if (i==0) {k=Tri[t][0]; l=Tri[t][1]; m=Tri[t][2];}
            if (i==1) {k=Tri[t][1]; l=Tri[t][2]; m=Tri[t][0];}
            if (i==2) {k=Tri[t][2]; l=Tri[t][0]; m=Tri[t][1];}              
            if (xy[k][2]==1){         
                
                for (verif=0,j=0;j<pos_conec_Tri[k];j++){
                    if (conec_Tri[k][j]==t) verif=1;
                }
                if (verif==0){
                    conec_Tri[k][pos_conec_Tri[k]]=t;
                    pos_conec_Tri[k]++;
                }           
                            
                
                for (verif=0,j=0;j<pos_conec[k];j++){
                    if (conec[k][j]==l) verif=1;
                }
                if (verif==0){
                    conec[k][pos_conec[k]]=l;
                    pos_conec[k]++;
                    if (pos_conec[k]>=max_conec_p){
                        printf("perigo!");
                                             
                    }
                }
                for (verif=0,j=0;j<pos_conec[k];j++){
                    if (conec[k][j]==m) verif=1;
                }
                if (verif==0){
                    conec[k][pos_conec[k]]=m;
                    pos_conec[k]++;
                    if (pos_conec[k]>=max_conec_p){
                        printf("perigo!");
                                             
                    }
                }
            }
        }
    } 
    
    /////////////Voronoi
    //printf("voronoi ");
    
    double *Ax,*Ay,*B,alpha,beta;
    double x1,y1,x2,y2,minba,maxba;
    double area_old_aux;
    double sed_thick_aux;
    
    long maxj,ii,jj; 
    
    Ax = Aloc_vector_real(max_conec_p);
    Ay = Aloc_vector_real(max_conec_p);
    B = Aloc_vector_real(max_conec_p);  
    
    for (i=0;i<n;i++){
        if (xy[i][2]==1){
            x1=xy[i][0]; y1=xy[i][1];
            area_old_aux = area_vor[i];
            area_vor[i]=0;
            
            for (j=0;j<pos_conec[i];j++){
                x2=xy[conec[i][j]][0]; y2=xy[conec[i][j]][1];
                Ax[j]=x2-x1;
                Ay[j]=y2-y1;
                B[j]=(x2*x2-x1*x1+y2*y2-y1*y1)/2;
                if (Ax[j]*x1+Ay[j]*y1>B[j]){
                    Ax[j]*=-1; Ay[j]*=-1; B[j]*=-1;            
                }
            }
            maxj=j;
            /*if (x1==minx){Ax[maxj]=-1; Ay[maxj]= 0; B[maxj]=-minx; maxj++;}
            if (x1==maxx){Ax[maxj]= 1; Ay[maxj]= 0; B[maxj]= maxx; maxj++;}
            if (y1==miny){Ax[maxj]= 0; Ay[maxj]=-1; B[maxj]=-miny; maxj++;}
            if (y1==maxx){Ax[maxj]= 0; Ay[maxj]= 1; B[maxj]= maxy; maxj++;}*/
            
            if (x1<minx+Edge*0.01){Ax[maxj]=-1; Ay[maxj]= 0; B[maxj]=-minx; maxj++;}
            if (x1>maxx-Edge*0.01){Ax[maxj]= 1; Ay[maxj]= 0; B[maxj]= maxx; maxj++;}
            if (y1<miny+Edge*0.01){Ax[maxj]= 0; Ay[maxj]=-1; B[maxj]=-miny; maxj++;}
            if (y1>maxy-Edge*0.01){Ax[maxj]= 0; Ay[maxj]= 1; B[maxj]= maxy; maxj++;}
            
            
            for (j=0;j<maxj;j++){
                minba = 1.0E50; maxba = -1.0E50;
                if (Ax[j]!=0){
                    for (jj=0;jj<maxj;jj++){
                        if (jj!=j){
                            alpha=Ay[jj]-Ay[j]*Ax[jj]/Ax[j];
                            beta =B[jj] -B[j] *Ax[jj]/Ax[j];
                            if (alpha>0){
                                if (beta/alpha<minba) minba = beta/alpha;
                            }
                            if (alpha<0){
                                if (beta/alpha>maxba) maxba = beta/alpha;
                            }
                        }
                    }
                    if (minba-maxba>0) aresta_vor[i][j]=minba-maxba;
                    else aresta_vor[i][j]=0;    
                    area_vor[i]+=(B[j]/fabs(Ax[j]))*aresta_vor[i][j]/2;  
                    aresta_vor[i][j]*=sqrt(Ax[j]*Ax[j]+Ay[j]*Ay[j])/fabs(Ax[j]);     
                }
                else {
                    for (jj=0;jj<maxj;jj++){
                        if (jj!=j){
                            alpha=Ax[jj]-Ax[j]*Ay[jj]/Ay[j];
                            beta =B[jj] -B[j] *Ay[jj]/Ay[j];
                            if (alpha>0){
                                if (beta/alpha<minba) minba = beta/alpha;
                            }
                            if (alpha<0){
                                if (beta/alpha>maxba) maxba = beta/alpha;
                            }
                        }
                    }
                    if (minba-maxba>0) aresta_vor[i][j]=minba-maxba;
                    else aresta_vor[i][j]=0;    
                    area_vor[i]+=(B[j]/fabs(Ay[j]))*aresta_vor[i][j]/2; 
                    aresta_vor[i][j]*=sqrt(Ax[j]*Ax[j]+Ay[j]*Ay[j])/fabs(Ay[j]);            
                }
            }
            
            if (pos_falha-pos_falha_ant!=0){
                sed_thick_aux = h_topo[i]-h_bed[i];            
                h_topo[i]=h_bed[i]+sed_thick_aux*area_old_aux/area_vor[i]; 
            }            
            
        }                          
    }
    
    
    for (i=0;i<n;i++){
        if (xy[i][2]==1){
            x1=xy[i][0]; y1=xy[i][1];
            for (j=0;j<pos_conec[i];j++){
                x2=xy[conec[i][j]][0]; y2=xy[conec[i][j]][1];
                dist_vor[i][j]=sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
            }
        }
    }
    
    
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
    
    /*FILE *f_vor;    
    double soma_area=0;
    long verif_zero;
    f_vor = fopen("voronoi.txt","w");    
    for (i=0;i<n;i++){
        fprintf(f_vor,"%f %f |",xy[i][0],xy[i][1]);
        fprintf(f_vor,"%f | ",area_vor[i]);
        for (j=0;j<pos_conec[i];j++){
            fprintf(f_vor,"%d ",conec[i][j]);
        }
        fprintf(f_vor," | ");
        for (j=0;j<pos_conec[i];j++){
            fprintf(f_vor,"%f ",aresta_vor[i][j]);
            if (aresta_vor[i][j]==0) verif_zero=1;
        }
        if (xy[i][0]==minx || xy[i][0]==maxx || xy[i][1]==miny || xy[i][1]==maxy){
            fprintf(f_vor,"*\n");
        }
        else fprintf(f_vor,"\n");
            
        soma_area+=area_vor[i];
    }    
    fprintf(f_vor,"\n%f",soma_area);
    fclose(f_vor);*/
      
    
    free(Ax);
    free(Ay);
    free(B);
    
} 
