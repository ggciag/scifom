//MODIFICADO  vR modificado aqui!!!! 

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>

extern long **Tri;
extern double **xy;
extern long tri;


extern long nodes;
extern long **conec;
extern long *pos_conec;

extern double *area_vor;
extern double **aresta_vor;
extern double **dist_vor;

extern long *ordem_fluvi;
extern long *direc_fluvi;
extern double *dist_fluvi;

extern double *Qr;
extern double *Qr_prov;
extern double *Qf;

extern double *Df;

extern long *basins;
extern long **global_basin;
extern double **h_global_basin;
extern long **global_index;
extern long *n_global;
extern long *global_flag;

extern double *h_topo;
extern double *h_bed;

extern long *lagos;

extern double dt;

extern double nivel;

extern double vR;
extern double Kf;
extern double ls;
extern double lb;

extern double *h_min;
extern long *pos_min;

extern long *cond_topo_modif;

extern long *stack_fluvi;

extern long n_sub_dt;

extern double *h_topo_prov;

extern double *Lf_vec;

extern double *vR_map;

extern long *mar;

extern double *lsr;
extern double *depth_lsr;
extern int nsr;
extern double **h_sr;
extern double *lsr_map;

extern int lith_flag;


void fluvi_min()
{
    long i,j,jj;
    double dist_aux;
    long cont_fluvi=nodes;
    double Qeqb;
    long verif,cont_min;
	
	long n_min;
	
	
	/*for (i=0;i<nodes;i++){
		if (mar[i]==0){
			if (h_topo[i]<nivel){
				nivel=h_topo[i]-1.0;
			}
		}
	}*/
	
	for (i=0;i<nodes;i++){
		if (h_topo[i]<-3000.0)
			mar[i]=1;
		else
			mar[i]=0;
	}
	
	long var=1;
	while (var>0){
		var=0;
		for (i=0;i<nodes;i++){
			if (mar[i]==0){
				if (h_topo[i]<nivel){
					for (j=0;j<pos_conec[i];j++){
						if (mar[conec[i][j]]==1){
							mar[i]=1;
							var++;
						}
					}
				}
			}
		}		
	}
	
	

    for (i=0;i<nodes;i++){          
        direc_fluvi[i]=i;
		Qr[i]=dt*vR_map[i]/n_sub_dt; /// sem area!!!!! já está no vRmap!!!

        Qf[i]=0;   
        Df[i]=0;  
        ordem_fluvi[i]=1;
        //if (h_topo[i]>nivel && cond_topo_modif[i]==1){      
		if (mar[i]==0 && cond_topo_modif[i]==1){      
			jj=i; 
            dist_aux=100; //qualquer numero;
            for (j=0;j<pos_conec[i];j++){
                if (h_topo[jj]>h_topo[conec[i][j]]){
                    if (aresta_vor[i][j]>0){
                        jj=conec[i][j];
                        dist_aux=dist_vor[i][j];
                    }
                }
            }
            direc_fluvi[i]=jj;
            dist_fluvi[i]=dist_aux;
        }
    }

	

    for (cont_min=2,i=0;i<nodes;i++){  
        //if (h_topo[i]<nivel || cond_topo_modif[i]==0){
        if (mar[i]==1 || cond_topo_modif[i]==0){
			basins[i]=1;  
        }
        else {
            if (direc_fluvi[i]==i){
                basins[i]=cont_min;
                h_min[cont_min]=1000000.0;
                pos_min[cont_min]=i;                
                cont_min++;
				
            }
            else basins[i]=0;
        }
    }
	
	n_min=cont_min;
	
    verif=0;
    while (verif==0){
        verif=1;
        for (i=0;i<nodes;i++){
            if (basins[i]==0){
                j = direc_fluvi[i];
                if (basins[j]!=0){
                    basins[i]=basins[j];
					
                }
                else verif=0;
            }
        }
    }
	
	
    
    /*for (i=0;i<nodes;i++){
        if (basins[i]!=1){
            cont_min=basins[i];
            for (j=0;j<pos_conec[i];j++){
                jj=conec[i][j];
                if (basins[jj]==1){
                    if (h_topo[jj]<h_min[cont_min]){
                        h_min[cont_min]=h_topo[jj];
                        direc_fluvi[pos_min[cont_min]]=jj;
                    }
                }
            }
        }
    }*/
	
	long cont_min2;
	
	long aux_n;
	
	long cont,c1,c2;
	
	double h_max;
	
	
	for (c1=2;c1<n_min;c1++) n_global[c1]=0;
	
	for (i=0;i<nodes;i++){
        if (basins[i]!=1){
            cont_min=basins[i];
            for (j=0;j<pos_conec[i];j++){
                jj=conec[i][j];
                if (basins[jj]!=cont_min){
					
					if (h_topo[jj]>h_topo[i]) h_max=h_topo[jj];
					else h_max=h_topo[i];
						
					
					cont_min2=basins[jj];
					
					c1=cont_min;
					c2=cont_min2;
					
					aux_n=-1;
					for (cont=0;cont<n_global[c1];cont++){
						if (global_basin[c1][cont]==c2)
							aux_n=cont;
					}
					if (aux_n==-1){
						aux_n=n_global[c1];
						global_basin[c1][aux_n]=c2;
						h_global_basin[c1][aux_n]=h_max;
						global_index[c1][aux_n]=jj;//						<----- j
						n_global[c1]++;	
					}
					else if(h_max<h_global_basin[c1][aux_n]){
						h_global_basin[c1][aux_n]=h_max;
						global_index[c1][aux_n]=jj;//						<----- j
					}
					
					c1=cont_min2;
					c2=cont_min;
					
					aux_n=-1;
					for (cont=0;cont<n_global[c1];cont++){
						if (global_basin[c1][cont]==c2)
							aux_n=cont;
					}
					if (aux_n==-1){
						aux_n=n_global[c1];
						global_basin[c1][aux_n]=c2;
						h_global_basin[c1][aux_n]=h_max;
						global_index[c1][aux_n]=i;//						<----- i
						n_global[c1]++;									 
					}
					else if(h_max<h_global_basin[c1][aux_n]){
						h_global_basin[c1][aux_n]=h_max;
						global_index[c1][aux_n]=i;//						<----- i
					}
                }
            }
        }
    }
	
	/*for (c1=2;c1<n_min;c1++){
		printf("%ld %ld: ",c1,n_global[c1]);
		for (cont=0;cont<n_global[c1];cont++){
			printf("%ld ",global_basin[c1][cont]);
		}
		printf("\n");
	}
	printf("\n");*/
	
	
	long flood_basin;
	long flood_index;
	double h_flood;
	
	for (c1=2;c1<n_min;c1++) global_flag[c1]=1;
	
	for (i=2;i<n_min;i++){
		h_flood=10000000.0;
		for (c1=2;c1<n_min;c1++){
			if (global_flag[c1]==1){
				for (cont=0;cont<n_global[c1];cont++){
					if (global_basin[c1][cont]==1){
						if (h_flood>h_global_basin[c1][cont]){
							h_flood=h_global_basin[c1][cont];
							flood_basin=c1;
							flood_index=global_index[c1][cont];
						}
					}
				}
			}
		}
		h_min[flood_basin]=h_flood;
		direc_fluvi[pos_min[flood_basin]]=flood_index;
		global_flag[flood_basin]=0;

		for (c1=2;c1<n_min;c1++){
			if (global_flag[c1]==1){
				for (cont=0;cont<n_global[c1];cont++){
					if (global_basin[c1][cont]==flood_basin){
						global_basin[c1][cont]=1;
					}
				}
			}
		}
	
		
		
		//h_min[cont_min]=h_topo[jj];
		//direc_fluvi[pos_min[cont_min]]=jj;
	}
	
	
	for (i=0;i<nodes;i++) lagos[i]=1;
	
	for (i=0;i<nodes;i++){
		if (basins[i]!=1){
			if (h_topo[i]<=h_min[basins[i]]){
				lagos[i]=0;
			}
		}
	}
    
    
    for (i=0;i<nodes;i++){        
        j = direc_fluvi[i];
        if (j!=i){
            ordem_fluvi[j]=0;                
        }        
    }
    
    long cont_stack=0;
    
    while (cont_fluvi>0){
        
        for (i=0;i<nodes;i++){
            if (ordem_fluvi[i]==1){
				stack_fluvi[cont_stack]=i;
				cont_stack++;
                j=direc_fluvi[i];                   
                
                if (j!=i){
                    Qr[j]+=Qr[i];
                }
                ordem_fluvi[i]=2;
                cont_fluvi--;
                
                //printf("%f\n",Df[i]);
            }
        }
        
        for (i=0;i<nodes;i++){
            if (ordem_fluvi[i]==0){
                ordem_fluvi[i]=1;
            }
        }
        
        for (i=0;i<nodes;i++){
            j = direc_fluvi[i];
            if (ordem_fluvi[i]!=2 && j!=i){
                ordem_fluvi[j]=0;
            }
        }
        
        
    }
	
	
	
	/*
	///m=0.5
	for (i=0;i<nodes;i++){
		Qr_prov[i]=sqrt(Qr[i]);
	}
	*/
	///m=1.0
	for (i=0;i<nodes;i++){
		Qr_prov[i]=Qr[i];
	}
	
	
	for (i=0;i<nodes;i++)
		h_topo_prov[i]=h_topo[i];

	for (i=0;i<nodes;i++){
		lsr_map[i]=1.0;
	}
	if (lith_flag==1){
		for (i=0;i<nodes;i++){
			for (int lit=0; lit<nsr; lit++){
				if (h_topo[i]-depth_lsr[lit] > h_sr[i][lit]) h_sr[i][lit] = h_topo[i]-depth_lsr[lit];
				if (h_sr[i][lit]>h_topo[i]) h_sr[i][lit]=h_topo[i];
			}
		}
		for (i=0;i<nodes;i++){
			for (int lit=0;lit<nsr;lit++){
				if (h_topo[i]==h_sr[i][lit]) lsr_map[i]=lsr[lit];
			}
		}
	}

	for (i=0;i<nodes;i++){
		lsr_map[i]*=ls;
		Qr_prov[i]*=Kf/dist_fluvi[i];
	}
	double hi,hj;
	
	for (cont=0;cont<n_sub_dt;cont++){
		
		for (i=0;i<nodes;i++){
			Qf[i]=0;
		}
		
		
		for (cont_stack=0;cont_stack<nodes;cont_stack++){
			i=stack_fluvi[cont_stack];
			j=direc_fluvi[i];
			
			if (j!=i){
				hi=h_topo_prov[i];
				hj=h_topo_prov[j];
				if (hi>hj){      
					if (hj>nivel){
						//Qeqb = Kf*Qr_prov[i]*(hi-hj)/dist_fluvi[i];
						Qeqb = Qr_prov[i]*(hi-hj);
					}
					else{
						//Qeqb = Kf*Qr_prov[i]*(hi-nivel)/dist_fluvi[i];
						Qeqb = Qr_prov[i]*(hi-nivel);
					}
					
					if (hi<nivel){
						//Qeqb = Kf*Qr_prov[i]*(hi-hj)/dist_fluvi[i];
						Qeqb = Qr_prov[i]*(hi-hj);
					}
						
						
					if (Qeqb<=Qf[i]){
						Df[i]=(Qf[i]-Qeqb)/area_vor[i];
						Qf[j]+=Qeqb;
					}
					else {
						if (h_bed[i]<hi)
							lb=lsr_map[i];
						else
							lb=Lf_vec[i];
						//if (hi>4500.0)
						//	lb = 10000.0;
						Df[i]=((Qf[i]-Qeqb)/area_vor[i])*(dist_fluvi[i]/lb);
						Qf[j]+=Qf[i]+(Qeqb-Qf[i])*(dist_fluvi[i]/lb);
					}               
				}    
				else {
					Df[i]=Qf[i]/area_vor[i]; 
				}     
			}
			else {
				Df[i]=Qf[i]/area_vor[i]; 
			}
			
		}
		
		for (i=0;i<nodes;i++){
			h_topo_prov[i]+=cond_topo_modif[i]*Df[i];
		}
	
	}
	
	
	for (i=0;i<nodes;i++){
		Df[i]=h_topo_prov[i]-h_topo[i];
	}

	for (i=0;i<nodes;i++){
		lsr_map[i]/=ls;
	}
	
                 
    
}
    
                                    
           
         
