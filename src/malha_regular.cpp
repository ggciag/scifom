///!!!!!! Borda livre para x_min!!!!!!


#include <stdio.h>
#include <stdlib.h>

double **Aloc_matrix_real (long p, long n);
long **Aloc_matrix_long (long p, long n);
double *Aloc_vector_real (long n);
long *Aloc_vector_long (long n);

extern double TeConstante;
extern double TeConstante2;

extern double RHOM;

extern long tri_flex;
extern long **Tri_flex;
extern long nodes_flex;
extern double **xy_flex;

extern long **cond_c;

extern double **Kflex;
extern double **Kflex_c;
extern long **Kconec;
extern long *Kposconec;

extern double *Kdiag;
extern double *Kdiag_c;

extern double **Ke;
extern double *Te;
extern double *props;

//extern double *qflex;
//extern double *qflex_a;
extern double *bflex;
extern double *uflex;
extern double *wflex;
extern double *wflex_aux;
extern double *wflex_cumula;
extern double *wflex_fault;

extern double *IsoT;

extern long max_conec_p;

extern double area_ele_flex;

extern double dx_flex;
extern double dy_flex;

extern double minx_flex;
extern double miny_flex;

extern double *temper_q;
extern double *temper_q_a;

extern double *moho_flex;
extern double *moho_flex_aux;

extern double *Temper35;


void malha_regular(double minx,double maxx,double miny,double maxy,long Nx, long Ny)
{
	printf("%f %f %f %f\n",minx,maxx,miny,maxy);
    if (maxx<=minx) {
        printf("Erro nas dimensoes x da malha de elementos para o problema flexural\n");

    }
    if (maxy<=miny) {
        printf("Erro nas dimensoes y da malha de elementos para o problema flexural\n");

    }
    if (Nx<=1) {
        printf("Nx nao valido\n");

    }
    if (Ny<=1) {
        printf("Ny nao valido\n");

    }

    minx_flex=minx;
    miny_flex=miny;

    long i,j,k,l,m,t,tri_aux;
    long verif;

	double dx,dy;
	dx = (maxx-minx)/(Nx-1);
	dy = (maxy-miny)/(Ny-1);

	dx_flex = dx;
	dy_flex = dy;

	FILE *saida;
	saida = fopen("malha_flex.txt","w");
	fprintf(saida,"%ld\n",Nx*Ny);

	tri_flex = (Nx-1)*(Ny-1)*2;
	Tri_flex = Aloc_matrix_long(tri_flex,3);
    nodes_flex = Nx*Ny;
    xy_flex = Aloc_matrix_real(nodes_flex,2);

	temper_q = Aloc_vector_real(nodes_flex);
	temper_q_a = Aloc_vector_real(nodes_flex);

	moho_flex = Aloc_vector_real(nodes_flex);
	moho_flex_aux = Aloc_vector_real(nodes_flex);


	Temper35 = Aloc_vector_real(nodes_flex);

    area_ele_flex = (maxx-minx)*(maxy-miny)/tri_flex;

    cond_c = Aloc_matrix_long(nodes_flex,3);

	for (i=0;i<Nx;i++){
		for (j=0;j<Ny;j++){
			fprintf (saida,"%ld %f %f ", i*Ny+j,dx*i+minx,dy*j+miny);  // 0 ate (Nx-1)+Nx*(Ny-1) = Nx*Ny-1
			xy_flex[i*Ny+j][0]=dx*i+minx;
			xy_flex[i*Ny+j][1]=dy*j+miny;
            fprintf(saida,"0 ");
            cond_c[i*Ny+j][0]=0;
            if (j==0 || j==Ny-1){
                fprintf(saida,"1 ");
                cond_c[i*Ny+j][1]=1;
            }
            else {
                fprintf(saida,"0 ");
                cond_c[i*Ny+j][1]=0;
            }
			if (i==Nx-1){//////////////////Borda livre para i==0
				fprintf(saida,"1\n");
				cond_c[i*Ny+j][2]=1;
            }
			else {
				fprintf(saida,"0\n");
				cond_c[i*Ny+j][2]=0;
            }
		}
	}

	for (i=0;i<Nx;i++){
		for (j=0;j<Ny;j++){
			fprintf (saida,"%ld ", i*Ny+j);
			fprintf(saida,"0\n");
		}
	}

	fprintf(saida,"%ld\n",(Nx-1)*(Ny-1)*2);

    tri_aux=0;
	for (i=0;i<Nx-1;i++){
		for (j=0;j<Ny-1;j++){
            Tri_flex[tri_aux][0]=i*Ny+j;
            Tri_flex[tri_aux][1]=(i+1)*Ny+j;
            Tri_flex[tri_aux][2]=(i+1)*Ny+(j+1);
            tri_aux++;
			fprintf(saida,"%ld %ld %ld 0 %f\n",i*Ny+j,(i+1)*Ny+j,(i+1)*Ny+(j+1),TeConstante);
			Tri_flex[tri_aux][0]=i*Ny+j;
            Tri_flex[tri_aux][1]=(i+1)*Ny+(j+1);
            Tri_flex[tri_aux][2]=i*Ny+(j+1);
            tri_aux++;
			fprintf(saida,"%ld %ld %ld 0 %f\n",i*Ny+j,(i+1)*Ny+(j+1),i*Ny+(j+1),TeConstante);
		}
	}


	props = Aloc_vector_real(7);
	fprintf(saida,"1\n1E+11 0.25 3300 3300 9.8");

	props[0]=1E+11;
	props[1]=0.25;
	props[2]=RHOM;
	props[3]=RHOM;
	props[4]=9.8;


	fclose(saida);

	Kconec = Aloc_matrix_long (nodes_flex,max_conec_p);
	Kposconec = Aloc_vector_long (nodes_flex);
	Kflex = Aloc_matrix_real (3*nodes_flex,3*max_conec_p);
	Kflex_c= Aloc_matrix_real (3*nodes_flex,3*max_conec_p);
	Ke = Aloc_matrix_real(9,9);

	Kdiag = Aloc_vector_real (3*nodes_flex);
	Kdiag_c= Aloc_vector_real (3*nodes_flex);
	//qflex = Aloc_vector_real (3*nodes_flex);
    //qflex_a = Aloc_vector_real (3*nodes_flex);
	bflex = Aloc_vector_real (3*nodes_flex);
	uflex = Aloc_vector_real (3*nodes_flex);
	wflex = Aloc_vector_real (3*nodes_flex);
	wflex_fault = Aloc_vector_real (3*nodes_flex);
	wflex_aux = Aloc_vector_real (nodes_flex);
	wflex_cumula = Aloc_vector_real (nodes_flex);


    IsoT = Aloc_vector_real (nodes_flex);

	Te = Aloc_vector_real(tri_flex);
	for (t=0;t<tri_flex;t++){
        if ((xy_flex[Tri_flex[t][0]][1]+
             xy_flex[Tri_flex[t][1]][1]+
             xy_flex[Tri_flex[t][2]][1])/3>(maxy+miny)/2)
             Te[t]=TeConstante;
        else
             Te[t]=TeConstante2;
    }

    for (t=0;t<nodes_flex;t++){
        Kposconec[t]=1;
        Kconec[t][0]=t;
    }
	for (t=0;t<tri_flex;t++){
        for (i=0;i<3;i++){
            if (i==0) {k=Tri_flex[t][0]; l=Tri_flex[t][1]; m=Tri_flex[t][2];}
            if (i==1) {k=Tri_flex[t][1]; l=Tri_flex[t][2]; m=Tri_flex[t][0];}
            if (i==2) {k=Tri_flex[t][2]; l=Tri_flex[t][0]; m=Tri_flex[t][1];}
            for (verif=0,j=0;j<Kposconec[k];j++){
                if (Kconec[k][j]==l) verif=1;
            }
            if (verif==0){
                Kconec[k][Kposconec[k]]=l;
                Kposconec[k]++;
                if (Kposconec[k]>=max_conec_p){
                    printf("perigo!");

                }
            }
            for (verif=0,j=0;j<Kposconec[k];j++){
                if (Kconec[k][j]==m) verif=1;
            }
            if (verif==0){
                Kconec[k][Kposconec[k]]=m;
                Kposconec[k]++;
                if (Kposconec[k]>=max_conec_p){
                    printf("perigo!");

                }
            }
        }
    }

    FILE *F_conec;

    F_conec = fopen("Flex_conec.txt","w");
    for (t=0;t<nodes_flex;t++){
        fprintf(F_conec,"%ld | ",Kposconec[t]);
        for (i=0;i<Kposconec[t];i++){
            fprintf(F_conec,"%ld ",Kconec[t][i]);
        }
        fprintf(F_conec,"\n");
    }
    fclose(F_conec);


}
