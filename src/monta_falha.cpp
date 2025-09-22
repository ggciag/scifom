#include <math.h>
#include <stdio.h>

extern double **pos_falha_vec;
extern long cont_falha;
extern long num_falha;
extern double pos_falha;
extern double pos_falha_ant;
extern double tempo;

extern double *h_foot;
extern double *h_bed;
extern double H_brit;
extern double **xy;

extern double axis_stream;

extern long nodes;

extern double **falha_plot; //posição da falha em y = 0
extern long **falha_plot_pos;

double var_fault(double v);

extern long *side;

void monta_falha()
{
    long aux;
	double escala_brit=30;
    if (cont_falha<num_falha){
        if (tempo >=pos_falha_vec[cont_falha][2] && pos_falha_vec[cont_falha][4]==0){
            pos_falha_vec[cont_falha][4]=1;      
            pos_falha =     pos_falha_vec[cont_falha][0];
            pos_falha_ant = pos_falha_vec[cont_falha][0];
            double dip =    pos_falha_vec[cont_falha][1];   
			
			//axis_stream = pos_falha;
			
            
            long i;
            if (dip<0) {
                for (i=0;i<nodes;i++){      
                    if ((xy[i][0]-pos_falha-var_fault(xy[i][1]))/(H_brit*escala_brit/fabs(dip))>-1.0){
                        h_foot[i]=dip*(pos_falha+var_fault(xy[i][1])-xy[i][0]);
						/*if (xy[i][0]-pos_falha-var_fault(xy[i][1])>0 && h_bed[i]>h_foot[i]){
							h_foot[i]=h_bed[i]+0.01;
						}
						if (xy[i][0]-pos_falha-var_fault(xy[i][1])<0 && h_bed[i]<h_foot[i]){
							h_foot[i]=h_bed[i]-0.01;
						}*/
                        if (xy[i][1]==0 && h_foot[i]<h_bed[i] && h_foot[i]>-H_brit+5000.0){
                            aux = falha_plot_pos[cont_falha][0]+1;
                            falha_plot[cont_falha][aux]=h_foot[i];
                            falha_plot_pos[cont_falha][aux]=i;
                            falha_plot_pos[cont_falha][0]++;
                        }
                    }
                    else
                       h_foot[i]=-H_brit*escala_brit; 
                    if (h_foot[i]>40000.0) h_foot[i]=40000;
                }
            }
            else {                
                for (i=0;i<nodes;i++){ 
                    if ((xy[i][0]-pos_falha-var_fault(xy[i][1]))/(H_brit*escala_brit/fabs(dip))<1.0){
                        h_foot[i]=dip*(pos_falha+var_fault(xy[i][1])-xy[i][0]);
						/*if (xy[i][0]-pos_falha-var_fault(xy[i][1])<0 && h_bed[i]>h_foot[i]){
							h_foot[i]=h_bed[i]+0.01;
						}
						if (xy[i][0]-pos_falha-var_fault(xy[i][1])>0 && h_bed[i]<h_foot[i]){
							h_foot[i]=h_bed[i]-0.01;
						}*/
                        if (xy[i][1]==0 && h_foot[i]<h_bed[i] && h_foot[i]>-H_brit+5000.0){
                            aux = falha_plot_pos[cont_falha][0]+1;
                            falha_plot[cont_falha][aux]=h_foot[i];
                            falha_plot_pos[cont_falha][aux]=i;
                            falha_plot_pos[cont_falha][0]++;
                        }
                    }
                    else
                        h_foot[i]=-H_brit*escala_brit;
                    if (h_foot[i]>40000.0) h_foot[i]=40000.0;
                } 
					
				
                printf("falha: %ld %ld\n", falha_plot_pos[cont_falha][0],cont_falha);
            }
			for (i=0;i<nodes;i++){
				if (h_bed[i]>h_foot[i]) side[i]=1;
				else					side[i]=-1;
			}
        }
        if (tempo > pos_falha_vec[cont_falha][3]) cont_falha++;
    }          
}
                
                
                
                
