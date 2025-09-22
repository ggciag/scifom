/*
 *  var_fault_stream.cpp
 *  F5_8_2
 *
 *  Created by Victor Sacek on 13/09/10.
 *  Copyright 2010 IAG/USP. All rights reserved.
 *
 */


#include <math.h>

extern double **pos_falha_vec;
extern long cont_falha;

double var_fault_stream(double v)
{
	return(0.0);
	//return(-v*35.0/200.0-v*0.2+20000.0);//teste2
	//return(-v*35.0/200.0-v*0.2+20000.0+(-v*30.0/200.0+30000.0));
	
}
