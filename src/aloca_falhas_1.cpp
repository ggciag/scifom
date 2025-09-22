double **Aloc_matrix_real (long p, long n);
long **Aloc_matrix_long (long p, long n);

extern double maxx;
extern double minx;

extern double tgdip;

extern long num_falha;
extern double **pos_falha_vec;

extern double **falha_plot; //posi‹o da falha em y = 0
extern long **falha_plot_pos;

void aloca_falhas()
{
    num_falha = 4;
    
    pos_falha_vec = Aloc_matrix_real(num_falha,5);
    
    falha_plot = Aloc_matrix_real(num_falha,100);
    falha_plot_pos = Aloc_matrix_long(num_falha,100);        
    
	//double space_shift=160000.0;
	double space_shift=70000.0;
	
    /*pos_falha_vec[0][0]= space_shift+minx+105000.0;
	 pos_falha_vec[1][0]= space_shift+minx+80000.0;
	 pos_falha_vec[2][0]= space_shift+minx+60000.0;
	 pos_falha_vec[3][0]= space_shift+minx+45000.0; 
	 pos_falha_vec[4][0]= space_shift+minx+30000.0;*/
    

	if (num_falha>0) pos_falha_vec[0][0]= space_shift+minx+60000.0;
	if (num_falha>1) pos_falha_vec[1][0]= space_shift+minx+80000.0;
	if (num_falha>2) pos_falha_vec[2][0]= space_shift+minx+100000.0;
	if (num_falha>3) pos_falha_vec[3][0]= space_shift+minx-200000000.0;
	if (num_falha>4) pos_falha_vec[4][0]= space_shift+minx+60000.0;
	/*if (num_falha>5) pos_falha_vec[5][0]= space_shift+minx+60000.0;
	if (num_falha>6) pos_falha_vec[6][0]= space_shift+minx+60000.0;
	if (num_falha>7) pos_falha_vec[7][0]= space_shift+minx+60000.0;
	if (num_falha>8) pos_falha_vec[8][0]= space_shift+minx+60000.0;
	if (num_falha>9) pos_falha_vec[9][0]= space_shift+minx+60000.0;*/
	
    if (num_falha>0) pos_falha_vec[0][1]= -tgdip;
    if (num_falha>1) pos_falha_vec[1][1]= -tgdip;
    if (num_falha>2) pos_falha_vec[2][1]= -tgdip;
	if (num_falha>3) pos_falha_vec[3][1]= tgdip;
    if (num_falha>4) pos_falha_vec[4][1]= -tgdip;
	/*if (num_falha>5) pos_falha_vec[5][1]= tgdip;
	if (num_falha>6) pos_falha_vec[6][1]= -tgdip;
	if (num_falha>7) pos_falha_vec[7][1]= tgdip;
	if (num_falha>8) pos_falha_vec[8][1]= -tgdip;
	if (num_falha>9) pos_falha_vec[9][1]= tgdip;*/
    
    double time_shift=1000.0E6, scale_time=10.0; 
    
    if (num_falha>0) pos_falha_vec[0][2]=    0000.0/scale_time+time_shift;
    if (num_falha>1) pos_falha_vec[1][2]=  4000000.0/scale_time+time_shift;
    if (num_falha>2) pos_falha_vec[2][2]=  8000000.0/scale_time+time_shift;
    if (num_falha>3) pos_falha_vec[3][2]= 12000000.0/scale_time+time_shift;
    if (num_falha>4) pos_falha_vec[4][2]= 16000000.0/scale_time+time_shift;
    /*if (num_falha>5) pos_falha_vec[5][2]= 3000000.0/scale_time+time_shift;
	if (num_falha>6) pos_falha_vec[6][2]= 4000000.0/scale_time+time_shift;
	if (num_falha>7) pos_falha_vec[7][2]= 5000000.0/scale_time+time_shift;
	if (num_falha>8) pos_falha_vec[8][2]= 6000000.0/scale_time+time_shift;
	if (num_falha>9) pos_falha_vec[9][2]= 7000000.0/scale_time+time_shift;*/
   
    
    if (num_falha>0) pos_falha_vec[0][3]=   4000000.0/scale_time+time_shift;
    if (num_falha>1) pos_falha_vec[1][3]=   8000000.0/scale_time+time_shift;
    if (num_falha>2) pos_falha_vec[2][3]=  12000000.0/scale_time+time_shift;
    if (num_falha>3) pos_falha_vec[3][3]=  6000000000.0/scale_time+time_shift;
    if (num_falha>4) pos_falha_vec[4][3]=  80000000.0/scale_time+time_shift;
	/*if (num_falha>5) pos_falha_vec[5][3]=  4000000.0/scale_time+time_shift;
	if (num_falha>6) pos_falha_vec[6][3]=  5000000.0/scale_time+time_shift;
	if (num_falha>7) pos_falha_vec[7][3]=  6000000.0/scale_time+time_shift;
	if (num_falha>8) pos_falha_vec[8][3]=  7000000.0/scale_time+time_shift;
	if (num_falha>9) pos_falha_vec[9][3]=  8000000.0/scale_time+time_shift;*/
    
}



