double maxx;//600000;
double minx;//200000;
double maxy;//=400000.0;//400000;
double miny;//=0.0;

double Edge;

long n_lat;
long n_latx;

double tempo_max;//=50.1E6;

long *side;

double *moho_flex;
double *moho_flex_aux;


long **Tri;
long tri;
long tri_max_aloca;
long **conec_Tri;
long *pos_conec_Tri;
long **conec;
long *pos_conec;
double **xy;
long nodes;
long nodes_old;
long nodes_max_aloca;
double *area_vor;
double *inv_area;
double **aresta_vor;
double **dist_vor;

double *h_topo;
double *h_bed;
double *h_flex_cumulat;
double *h_topo_din_cumulat;
double *hu_flex_cumulat;
double *water;


double *tempo_sea;
double *h_sea;
long n_sea_levels;

double *h_topo_prov;
// variaveis adicionadas para uplift variavel

double **uplift_map;
double *tempos_uplift_min;
double *tempos_uplift_max;
long numero_uplift;
double *uplift_factor;

// variaveis adicionadas para a topografia dinamica
double **topo_din_map;
double *tempo_topo_din;
long numero_intervalos;
double fac_topo_din;
//tempo que representa o inicio da simula��o (Ma)
double tempo_inicial;

int lith_flag;

//
double *Lf_vec;

long *lagos;

long *mar;

double *moho;
double *load_Temper;

long *cond_topo_modif;
//double *h_isost;
//double *h_original;


long n_falhas;
double **falhas_param;


double *h_q;
double *h_q_a;
double *h_w;

double *h_u;

double *Ds;
double *Df;

long *basins;
double *h_min;
long *pos_min;

long **global_basin;
double **h_global_basin;
long *n_global;
long **global_index;
long *global_flag;

long *ordem_fluvi;
long *direc_fluvi;
double *dist_fluvi;

long *stack_fluvi;

double *Qr;
double *Qr_prov;
double *Qf;

double vR;//=1;
double time_ofchangevR;
double vR2;

double vRandes;//=1;
double time_ofchangevRandes;
double vR2andes;

double **vR_maps;
double *h_vR_external;
int nvR_maps;
int vR_external_flag=0;

double Terigida;
double Teoffshore;

double *vR_aux;
long *vR_ordem;

double *vR_map;
double *vR_flow;

double windx=-1.0;
double windy=0.0;

double K_d=0.0;// = 0.3; /////alterado
double K_m;// = 200.0;
double Kf;// = 0.03;
double ls;// = 10000.0; ////alterado

double lb;// = 100000.0; ////alterado
double lb2;

double *lsr;
double *depth_lsr;
int nsr;
double **h_sr;
double *lsr_map;

double nivel;
double dt;
long n_sub_dt = 1000;

double tempo = 0;

double uplift_scale;

double time_ofchangeu;

double uplift_scale2;

long max_conec_p=40;

double topo_seno;


FILE *f_h_topo;
FILE *f_h_sed;



//////////////////////////////////////////
//Flexura

double TeConstante=5000.0;
double TeConstante2=5000.0;

double Telitho;

long Te_uc_on=1;

long tri_flex;
long **Tri_flex;
long nodes_flex;
double **xy_flex;

long **cond_c;


double **peso;
long **peso_pos;

double **Kflex;

double **Kflex_c;
long **Kconec;
long *Kposconec;

double **Ke;

double *Te;
double *Te_map;

double *props;
double *Kdiag;
double *Kdiag_c;

//double *qflex;
//double *qflex_a;
double *bflex;
double *wflex;
double *wflex_aux;
double *wflex_cumula;
double *uflex;
double *wflex_fault;

double *Temper35;

double *h_temper35;

double area_ele_flex;


double *IsoT;

long verif_first=0;

///////////////////////////

double Sed_esq=0;
double Sed_dir=0;
double Sed_in =0;




//////////////////////////

long tetra_thermal;
long **Tet_thermal;
long nodes_thermal;
double **xyz_thermal;
double **xyz_thermal_fut;

double **Kthermal;
double *Kthermal_diag;
double **Mthermal;
double *Mthermal_diag;
double **Kthermal1;
double *Kthermal1_diag;
double **Kthermal2;
double *Kthermal2_diag;

double *T_vec;
double *T_vec_fut;
double *T_vec0;


double *fthermal;

long *cond_thermal;

double **TKe;
double **TMe;

double *Temper;

long thermal_max_conec_p=40;

long **Kthermal_conec;
long *Kthermal_posconec;
long *pos_temper;

double *Dtemper;
long *cond_borda_Temper;
long *cond_borda_Temper2;

double *temper_q;
double *temper_q_a;

double **v_adv;

double *h_top;
double *v_adv2D;

double depth = 150000.0;

double h_bot = -depth;

double alpha_thermal=0.5;
double comp_alpha_thermal = 1.0 - alpha_thermal;

double H_C = 35000.0;
double H_brit = 12000.0;

//////////////////////////

double RHOW = 1030;
double RHOS = 2700;
double RHOC = 2700;
double RHOM = 3300;
double alpha_exp_thermo = 3.28E-5;
double kappa = 1.0E-6;
double seg_per_ano = 365.0*24.0*3600.0;


double dt_calor=20000.0;
double dt_calor_sec=dt_calor*seg_per_ano;


//////////////////////////

double soma_volume;

long Nx = 60;
long Ny = 42;
long Nz = 12;

long layers = Nz;

double dx_flex;
double dy_flex;

double minx_flex;
double miny_flex;

//////////////////////////

double **pos_falha_vec;
long cont_falha=0;
long num_falha;
double pos_falha;
double pos_falha_ant=pos_falha;
double inclina=0.5;

double **falha_plot; //posi��o da falha em y = 0
long **falha_plot_pos;

double tgdip = 1;//1.73;//12.0/20.0;

//double tempo_muda_falha=500000;

double V_meio=0.01;
double tempo_despl = 2000;
//double despl=20;
double despl=2*V_meio*tempo_despl;


double axis_stream=300000.0;

double tempo_max_stream=40.0E6;

////////////////////////////


double *h_foot;


double *h_crust_sup;

double *moho_aux;


long **tri_p;
long *cond_p;
long **aresta_p;

long ntri_lat;

long **tri_lat;

int ext_mesh=0;



