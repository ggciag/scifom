/*
 *  fault_analit.cpp
 *  F5_8_3
 *
 *  Created by Victor Sacek on 15/09/10.
 *  Copyright 2010 IAG/USP. All rights reserved.
 *
 */
#include <math.h>

extern double tempo_despl;
extern double **pos_falha_vec;
extern long cont_falha;
extern double despl;
extern double tgdip;
extern double TeConstante;
extern double *props;

extern double RHOC;

extern long *side;

extern double tempo;

extern double V_meio;

extern double H_brit;



double fault_analit(double x, long i){
	
	double w_max = tgdip*despl/2;
	
	double fator=2.0;
	
	double Te_analit = TeConstante*(1.0-pow(fator*tgdip*V_meio*2*(tempo-pos_falha_vec[cont_falha][2])/H_brit,2.0));
	
	if (Te_analit<10.0) 
		Te_analit=10.0;
	
	//if (tempo>1.0E6) Te_analit = 2000.0;
	
	double E     = props[0];
	//double v     = props[1];
	double g     = props[4];
	
	double D = E*Te_analit*Te_analit*Te_analit/12.0/(1-0.25*0.25);
	
	double alpha = pow(D/RHOC/g,0.25);
	
	if (x<0) x*=-1.0;
	
	if (pos_falha_vec[cont_falha][1]>0){
		if (side[i]*pos_falha_vec[cont_falha][1]<0)
			return (w_max*exp(-0.701*x/alpha)*cos(0.701*x/alpha));
		else 
			return (w_max*2.0-w_max*exp(-0.701*x/alpha)*cos(0.701*x/alpha));
	}
	else {
		if (side[i]*pos_falha_vec[cont_falha][1]>0)
			return (w_max*exp(-0.701*x/alpha)*cos(0.701*x/alpha));
		else 
			return (w_max*2.0-w_max*exp(-0.701*x/alpha)*cos(0.701*x/alpha));
	}

}