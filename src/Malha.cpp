///////////////// em i=0, topo nao muda


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <time.h>

double **Aloc_matrix_real (long p, long n);
long **Aloc_matrix_long (long p, long n);
double *Aloc_vector_real (long n);
long *Aloc_vector_long (long n);
int incircle(double **xy, long i, long j, long k, long l);
int intriangle(double **pontos, long i, long j, long k, double xx, double yy);
double det (double p1, double p2, double p3, double p4);

void malha_regular(double minx,double maxx,double miny,double maxy,long Nx, long Ny);

long gera_malha(long n, double **xy, long **Tri);
long importa_malha(long **Tri);

void thermal_aloca();
void thermal_modiftopo();
void thermal_monta_K();
void thermal_createf();
void thermal_firstK1();
void thermal_inic();

void stream();

extern double maxx;
extern double minx;
extern double maxy;
extern double miny;

extern long n_lat;
extern long n_latx;

extern double Edge;

extern double *h_isost;

extern long **Tri;
extern long tri;
extern long tri_max_aloca;
extern long **conec;
extern long *pos_conec;
extern long **conec_Tri;
extern long *pos_conec_Tri;
extern double **xy;
extern long nodes;
extern long nodes_max_aloca;
extern double *area_vor;
extern double **aresta_vor;
extern double **dist_vor;

extern long max_conec_p;

extern long tri_flex;
extern long **Tri_flex;
extern long nodes_flex;
extern double **xy_flex;

extern double **peso;
extern long **peso_pos;

extern long Nx;
extern long Ny;
extern long Nz;

extern double depth;

extern long *cond_topo_modif;

extern long *in_tri; // Numero do triangulo para o respectivo ponto novo

extern double tempo;

extern int ext_mesh;

void malha()
{
    
    long i,j,jj,maxj,k,l,m,ii;
    double x,y;    
    double area_tri,aresta,aresta_aux;
    double *Ax,*Ay,*B,alpha,beta;
    
    printf("Digite o numero de pontos da lateral\n");
    long n,t;
    //scanf("%d",&n_lat);
    
    n = n_lat*n_latx;
    nodes_max_aloca = n*3;
    if (n<4) {
        printf("Numero invalido de pontos");
        exit(-1);
    }
    nodes=n;
    tri=2*n;
    tri_max_aloca=tri*3;    
    
    
    area_tri = (maxx-minx)*(maxy-miny)/tri;
    //aresta = sqrt(4*1.73*area_tri/3);
    aresta = (maxy-miny)/(n_lat-1);

    printf("\nmaxx before: %f\n",maxx);
	
	minx = 100000.0;
	maxx = minx;//aresta*(n_latx-1)+minx;
    printf("\nmaxx after: %f\n",maxx);

    Edge = aresta;
    
    long verif;
    printf("Creating nodes\n");
	
	in_tri = Aloc_vector_long(n_lat);
	
    xy = Aloc_matrix_real(nodes_max_aloca,3);
    conec = Aloc_matrix_long(nodes_max_aloca,max_conec_p);
    pos_conec = Aloc_vector_long(nodes_max_aloca);
    conec_Tri = Aloc_matrix_long(nodes_max_aloca,max_conec_p);
    pos_conec_Tri = Aloc_vector_long(nodes_max_aloca);
    
    cond_topo_modif = Aloc_vector_long (nodes_max_aloca);
    
    dist_vor = Aloc_matrix_real(nodes_max_aloca,max_conec_p);
    
    Ax = Aloc_vector_real(max_conec_p);
    Ay = Aloc_vector_real(max_conec_p);
    B = Aloc_vector_real(max_conec_p);   
    
    area_vor = Aloc_vector_real(nodes_max_aloca);
    aresta_vor = Aloc_matrix_real(nodes_max_aloca,max_conec_p);
    
    Tri= Aloc_matrix_long(tri_max_aloca,4);
	
	
	
	FILE *f_pontos;
	
	f_pontos = fopen("pontos.txt","r");

    minx = 1.0E24;
    maxx = -1.0E24;
	
	for (i=0;i<n;i++){
		fscanf(f_pontos, "%lf %lf",&xy[i][0],&xy[i][1]);
		if (xy[i][0]!=minx){///////////// em i=0, topo nao muda
		cond_topo_modif[i]=1; 
		}
        if (xy[i][0]<minx) minx = xy[i][0];
        if (xy[i][0]>maxx) maxx = xy[i][0];
	}
	
	fclose(f_pontos);
	
	
    aresta_aux=aresta*aresta*0.7*0.7;
    
    printf("Aresta min: %f\n",sqrt(aresta_aux));
    
   
    srand(1);
	
	
	
	
    printf("%ld %ld\n",i,nodes);
    
    if (ext_mesh==0) {
        tri = gera_malha(n, xy, Tri);
    }
    else {
        tri = importa_malha(Tri);
    }
       
    long aux1;    
    double xi,xj,xk,yi,yj,yk,delta;
    
    
    for (t=0;t<tri;t++){  
        xi=xy[Tri[t][0]][0]; yi=xy[Tri[t][0]][1]; 
    	xj=xy[Tri[t][1]][0]; yj=xy[Tri[t][1]][1]; 
    	xk=xy[Tri[t][2]][0]; yk=xy[Tri[t][2]][1]; 
        delta=0.5*((xj*yk-yj*xk)+xi*(yj-yk)+yi*(xk-xj));  
        if (fabs(delta)<0.01*area_tri) {
            printf("\n\nZero! %ld %ld %ld\n",Tri[t][0],Tri[t][1],Tri[t][2]);
            printf("%f %f\n%f %f\n%f %f\n\n",xi,yi,xj,yj,xk,yk);
            if (t!=tri-1){
                Tri[t][0]=Tri[tri-1][0];
                Tri[t][1]=Tri[tri-1][1];
                Tri[t][2]=Tri[tri-1][2];
                t--; 
            }
            tri--;            
        }
    	if (delta<0) {
    	    aux1=Tri[t][0]; Tri[t][0]=Tri[t][2]; Tri[t][2]=aux1;    		 
        }
    }  
       
	
    printf("Numero de triangulos %ld\n",tri);
    printf("Numero de nos %ld\n",n);
    printf("t/n = %f", double(t)/double(n));
    
    for (t=0;t<tri;t++){        
        for (i=0;i<3;i++){
            if (i==0) {k=Tri[t][0]; l=Tri[t][1]; m=Tri[t][2];}
            if (i==1) {k=Tri[t][1]; l=Tri[t][2]; m=Tri[t][0];}
            if (i==2) {k=Tri[t][2]; l=Tri[t][0]; m=Tri[t][1];}   
            
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
    
    FILE *f_conec;    
    f_conec = fopen("conec.txt","w");    
    for (i=0;i<n;i++){
        for (j=0;j<pos_conec[i];j++){
            fprintf(f_conec,"%ld ",conec[i][j]);
        }
        fprintf(f_conec,"\n");
    }    
    fclose(f_conec);
    
    double x1,x2,y1,y2,minba,maxba;

    printf("\nLimits: %f %f %f %f\n",minx,maxx,miny,maxy);
    
    for (i=0;i<n;i++){
        x1=xy[i][0]; y1=xy[i][1];
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
        if (x1==minx){Ax[maxj]=-1; Ay[maxj]= 0; B[maxj]=-minx; maxj++;}
        if (x1==maxx){Ax[maxj]= 1; Ay[maxj]= 0; B[maxj]= maxx; maxj++;}
        if (y1==miny){Ax[maxj]= 0; Ay[maxj]=-1; B[maxj]=-miny; maxj++;}
        if (y1==maxy){Ax[maxj]= 0; Ay[maxj]= 1; B[maxj]= maxy; maxj++;}
        
        
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
        
                     
    }
    
    
    for (i=0;i<n;i++){
        x1=xy[i][0]; y1=xy[i][1];
        for (j=0;j<pos_conec[i];j++){
            x2=xy[conec[i][j]][0]; y2=xy[conec[i][j]][1];
            dist_vor[i][j]=sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
        }
    }
            
    
    long verif_zero=0;
    
    FILE *f_vor;    
    double soma_area=0;
    f_vor = fopen("voronoi0.txt","w");    
    for (i=0;i<n;i++){
        fprintf(f_vor,"%f %f |",xy[i][0],xy[i][1]);
        fprintf(f_vor,"%f | ",area_vor[i]);
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
    fclose(f_vor);
    printf("Verif_zero: %ld\n",verif_zero);
    
    printf("\nArea total calculada:%g",soma_area);///1000000.0);
    printf("\nArea total real:     %g\n",(maxx-minx)*(maxy-miny));///1000000.0);
    
    
    
    double bordaL=5.0;    
    malha_regular(minx-bordaL,maxx+bordaL,miny,maxy,Nx,Ny);   
	thermal_aloca();
	thermal_firstK1();
	thermal_inic();
	
	
	
    stream();
	printf("\n");
    
    //malha_regular(minx,maxx*2,miny,maxy,Nx,Ny);    
    //aloca_calor(minx,maxx*2,miny,maxy,depth,Nx,Ny,Nz);
    
    peso = Aloc_matrix_real(nodes_max_aloca,3);
    peso_pos = Aloc_matrix_long(nodes_max_aloca,3);
    
    double a1,a2,a3,b1,b2,b3,c1,c2,c3;
          
          
    printf("Calcula peso: ");  
    for (i=0;i<nodes;i++){
        x = xy[i][0];
        y = xy[i][1];
        if (y>maxy-Edge*0.0001) y-=0.0001*Edge;
        if (y<miny+Edge*0.0001) y+=0.0001*Edge;
        for (verif=0,ii=0;ii<tri_flex && verif==0;ii++){
            if (intriangle(xy_flex,Tri_flex[ii][0],Tri_flex[ii][1],Tri_flex[ii][2],
                           x,y)==0){
                                                  
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
        		
        		verif=1;
            }            
        }
        if (verif==0){
            printf("out: %f %f %ld\n", x,y,i);
                       
        }
    }
    printf("Fim\n");
    
    
}



int incircle(double **xy, long i, long j, long k, long l)
{
	double Ax,Ay,Bx,By,Cx,Cy,Dx,Dy,aux;
	double M[3][3];
	double Det1,Det2;
	int verif;
	Ax = xy[i][0]; Ay = xy[i][1];
	Bx = xy[j][0]; By = xy[j][1];
	Cx = xy[k][0]; Cy = xy[k][1];
	Dx = xy[l][0]; Dy = xy[l][1];

	Det1 = Ax*By + Ay*Cx + Bx*Cy - By*Cx - Ax*Cy - Ay*Bx;

	if(Det1<0){
		aux=Bx; Bx=Cx; Cx=aux;
		aux=By; By=Cy; Cy=aux;
	}


	M[0][0]=Ax-Dx;
	M[1][0]=Bx-Dx;
	M[2][0]=Cx-Dx;

	M[0][1]=Ay-Dy;
	M[1][1]=By-Dy;
	M[2][1]=Cy-Dy;

	M[0][2]=(Ax-Dx)*(Ax-Dx)+(Ay-Dy)*(Ay-Dy);
	M[1][2]=(Bx-Dx)*(Bx-Dx)+(By-Dy)*(By-Dy);
	M[2][2]=(Cx-Dx)*(Cx-Dx)+(Cy-Dy)*(Cy-Dy);

	Det2 =  M[0][0]*M[1][1]*M[2][2];
	Det2 += M[0][1]*M[1][2]*M[2][0];
	Det2 += M[0][2]*M[1][0]*M[2][1];

	Det2 += -M[0][2]*M[1][1]*M[2][0];
	Det2 += -M[0][0]*M[1][2]*M[2][1];
	Det2 += -M[0][1]*M[1][0]*M[2][2];

	//printf("Det: %g \n",Det2);


	if(Det2<0){
		verif=1;
	}
	else{
		verif=0;
	}

	return(verif);
}

int intriangle(double **pontos, long i, long j, long k, double xx, double yy)
{
	double v[2],v0[2],v1[2],v2[2];
	double a,b;
		
	v[0]  = xx;				                 v[1]  = yy;
	v0[0] = pontos[i][0];				     v0[1] = pontos[i][1];
	v1[0] = pontos[j][0]-pontos[i][0];		 v1[1] = pontos[j][1]-pontos[i][1];
	v2[0] = pontos[k][0]-pontos[i][0];		 v2[1] = pontos[k][1]-pontos[i][1];
	a =  (det(v[0],v2[0],v[1],v2[1])-det(v0[0],v2[0],v0[1],v2[1]))/det(v1[0],v2[0],v1[1],v2[1]);
	b = -(det(v[0],v1[0],v[1],v1[1])-det(v0[0],v1[0],v0[1],v1[1]))/det(v1[0],v2[0],v1[1],v2[1]);
	
    
	if (a>=0 && b>=0 && a+b<=1)	return(0);
	else return(1);		

}
double det (double p1, double p2, double p3, double p4)
{
	return(p1*p4-p2*p3);
}
