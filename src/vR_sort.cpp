/*
 *  vR_sort.cpp
 *  Orogenia1.17_Te_variavel_3_m0.5
 *
 *  Created by Victor Sacek on 08/04/13.
 *  Copyright 2013 IAG/USP. All rights reserved.
 *
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#define NSTACK 50
#define M 7
#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;
#define SWAPint(a,b) inttemp=(a);(a)=(b);(b)=inttemp;

double **Aloc_matrix_real (long p, long n);
long **Aloc_matrix_long (long p, long n);
double *Aloc_vector_real (long n);
long *Aloc_vector_long (long n);

extern double *vR_aux;
extern long *vR_ordem;

extern double windx;
extern double windy;

extern long nodes;
extern long nodes_max_aloca;

extern double **xy;

extern double *vR_map;
extern double *vR_flow;

extern long **conec;
extern long *pos_conec;

extern double **aresta_vor;

extern double **vR_maps;
extern double *h_vR_external;
extern int nvR_maps;
extern int vR_external_flag;

void sort(unsigned long n, double arr[], long *Ox);

void vR_sort()
{
	long i,ii,j,jj;
	
	double wind_aux = sqrt(windx*windx+windy*windy);
	
	windx = windx/wind_aux;
	windy = windy/wind_aux;
	
	vR_aux = Aloc_vector_real(nodes_max_aloca);
	vR_ordem = Aloc_vector_long(nodes_max_aloca);
	
	vR_map = Aloc_vector_real(nodes_max_aloca);
	vR_flow = Aloc_vector_real(nodes_max_aloca);
	
	for (i=0;i<nodes;i++){
		vR_aux[i+1] = -xy[i][0]*windx - xy[i][1]*windy;
		vR_ordem[i+1] = i;
	}
	
	printf("vR_aux: %.1lf %.1lf %.1lf %.1lf %.1lf %.1lf\n",vR_aux[0],vR_aux[1],vR_aux[2],vR_aux[3],vR_aux[4],vR_aux[5]);
	
	sort(nodes, vR_aux, vR_ordem);
	
	for (i=0;i<nodes;i++){
		vR_aux[i] = -xy[i][0]*windx - xy[i][1]*windy;
		vR_ordem[i]=vR_ordem[i+1];
	}
	
	
	printf("vR_aux: %.1lf %.1lf %.1lf %.1lf %.1lf %.1lf %.1lf %.1lf\n",vR_aux[0],vR_aux[1],vR_aux[2],vR_aux[3],vR_aux[4],vR_aux[5],vR_aux[6],vR_aux[7]);
	printf("vR_ordem: %ld %ld %ld %ld %ld %ld %ld %ld\n",vR_ordem[0],vR_ordem[1],vR_ordem[2],vR_ordem[3],vR_ordem[4],vR_ordem[5],vR_ordem[6],vR_ordem[7]);
	for (ii=nodes-1;ii>=0;ii--){
		i=vR_ordem[ii];
		jj=i; 
		for (j=0;j<pos_conec[i];j++){
			if (vR_aux[jj]>vR_aux[conec[i][j]]){
				if (aresta_vor[i][j]>0){
					jj=conec[i][j];
				}
			}
		}
		vR_map[jj]=0;
	}
	
	long soma=0;
	
	printf("%ld %ld\n",soma,nodes);
	
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

void read_vR_external()
{
	FILE *f;
	nvR_maps=0;
	f = fopen("vR_external.txt","r");
	if (f!=NULL){
		
		fscanf(f,"%d",&nvR_maps);
		if (nvR_maps>0){
			vR_external_flag=1;
			printf("\n\nExternal vR on\n\n");
			h_vR_external = Aloc_vector_real(nvR_maps);
			vR_maps = Aloc_matrix_real(nvR_maps,nodes_max_aloca);
			
			for (int t=0;t<nvR_maps;t++){
				fscanf(f,"%lf",&h_vR_external[t]);
				printf("%lf\n",h_vR_external[t]);
			}
			for (int i=0; i<nodes; i++){
				for (int t=0;t<nvR_maps;t++){
					fscanf(f,"%lf",&vR_maps[t][i]);
				}
			}
		}
		else {
			vR_external_flag=-1;
			fscanf(f,"%lf %lf",&windx,&windy);
		}
		fclose(f);
	}
	else {
		vR_external_flag=0;
		printf("\n\nExternal vR off\n\n");
	}
}

extern double *tempo_sea;
extern double *h_sea;
extern long n_sea_levels;

void read_sea_level()
{
	FILE *f;
	n_sea_levels=0;
	f = fopen("sea_level_external.txt","r");
	if (f!=NULL){
		
		fscanf(f,"%ld",&n_sea_levels);
		if (n_sea_levels>0){
			printf("\n\nExternal sea level on\n\n");
			h_sea = Aloc_vector_real(n_sea_levels);
			tempo_sea = Aloc_vector_real(n_sea_levels);
			
			for (int t=0;t<n_sea_levels;t++){
				fscanf(f,"%lf %lf",&tempo_sea[t],&h_sea[t]);
				//printf("%lf\n",h_vR_external[t]);
			}
		}
		fclose(f);
	}
	else {
		printf("\n\nExternal sea level off\n\n");
	}
}