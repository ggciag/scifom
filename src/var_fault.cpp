#include <math.h>

double var_fault_stream(double v);

extern double **pos_falha_vec;
extern long cont_falha;

double var_fault(double v)
{
    return(0);
	
	
	if (cont_falha==0) {
		if (v<102500.0) return(7000.0);
		else {
			return(-20000.0);
		}
	}
	if (cont_falha==1) {
		if (v<102500.0) return(5000.0);
		else {
			return(-10000.0);
		}
	}
	if (cont_falha==2) return(0.0);
	if (cont_falha==3) {
		if (v<102500.0) return(10000.0);
		else {
			return(-35000.0);
		}
	}
	
	if (cont_falha==4) {
		if (v<102500.0) return(40000.0);
		else {
			return(15000.0);
		}
	}

	
	


	
	//return(0.0);

	
	/*double p1=-30000.0,p2=30000.0;
    double l1=50000.0,l2=150000.0;
	
	if (cont_falha==1) {
		l1 = 20000.0;
		l2 = 120000.0;
	}
	
	
	if (v>l2)
		return(p2);
	else {
	   if (v>l1){
		  double v_aux = p1+(p2-p1)*(v-l1)/(l2-l1);
		  return(v_aux);
	   }
	   else return(p1); 
	}*/
}
