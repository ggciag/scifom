#include <stdio.h>
#include <stdlib.h>
#include <string.h>

extern double maxx;
extern double minx;
extern double maxy;
extern double miny;


extern long n_lat;
extern long n_latx;

extern double axis_stream;
extern double Terigida;
extern double Teoffshore;


extern double vR;
extern double time_ofchangevR;
extern double vR2;

extern double vRandes;
extern double time_ofchangevRandes;
extern double vR2andes;

extern double Kf;

extern double K_m;
extern double K_d;

extern double ls;
extern double lb;
extern double lb2;

extern double uplift_scale;

extern double time_ofchangeu;
extern double uplift_scale2;
extern double tempo_max;
extern double dt;

extern long Nx;
extern long Ny;

extern long n_sub_dt;

void reader(){
    FILE *entra_var;
	entra_var = fopen("param_OrogSedFlex_1.1.txt", "r");

    int nline=0;
	int size = 1024;
	void *nread;
	char *tkn_w, *tkn_v;
	char line[size];

    while (!feof(entra_var)){

        // Increment line number nline
        nline += 1;

        // Read each line from the parameters f_parameters, trim undesireble characters
        // and split line in thk_w[] and thk_v[].
        nread = fgets(line, size, entra_var);
        if ((nread == NULL) || (line[0] == '\n') || (line[0] == '#')) continue;
        tkn_w 	= strtok(line, " \t=\n");
        if (tkn_w[0] == '#') continue;
        tkn_v 	= strtok(NULL, " \t=#\n");

        if (strcmp(tkn_w, "maxy") == 0) {maxy = atof(tkn_v);}
        //fscanf(entra_var,"%lf",&maxy);
        else if (strcmp(tkn_w, "miny") == 0) {miny = atof(tkn_v);}
        //fscanf(entra_var,"%lf",&miny);

        else if (strcmp(tkn_w, "ny") == 0) {n_lat = atoi(tkn_v);}
        //fscanf(entra_var,"%ld",&n_lat);
        else if (strcmp(tkn_w, "nx") == 0) {n_latx = atoi(tkn_v);}
        //fscanf(entra_var,"%ld",&n_latx);

        else if (strcmp(tkn_w, "ny_flexural") == 0) {Ny = atoi(tkn_v);}
        
        else if (strcmp(tkn_w, "nx_flexural") == 0) {Nx= atoi(tkn_v);}
        

        else if (strcmp(tkn_w, "axis_stream") == 0) {axis_stream = atof(tkn_v);}
        //fscanf(entra_var, "%lf",&axis_stream);
        else if (strcmp(tkn_w, "Te_rigida") == 0) {Terigida = atof(tkn_v);}
        //fscanf(entra_var, "%lf",&Terigida);
        else if (strcmp(tkn_w, "Te_offshore") == 0) {Teoffshore = atof(tkn_v);}
        //fscanf(entra_var, "%lf",&Teoffshore);
        //TeConstante2=TeConstante;
        //fscanf(entra_var, "%lf",&Telitho);
        //Telitho = TeConstante;


        else if (strcmp(tkn_w, "vR") == 0) {vR = atof(tkn_v);}
        //fscanf(entra_var,"%lf",&vR);
        else if (strcmp(tkn_w, "time_ofchangevR") == 0) {time_ofchangevR = atof(tkn_v);}
        //fscanf(entra_var,"%lf",&time_ofchangevR);
        else if (strcmp(tkn_w, "vR2") == 0) {vR2 = atof(tkn_v);}
        //fscanf(entra_var,"%lf",&vR2);

        else if (strcmp(tkn_w, "vRandes") == 0) {vRandes = atof(tkn_v);}
        //fscanf(entra_var,"%lf",&vRandes);
        else if (strcmp(tkn_w, "time_ofchangevRandes") == 0) {time_ofchangevRandes = atof(tkn_v);}
        //fscanf(entra_var,"%lf",&time_ofchangevRandes);
        else if (strcmp(tkn_w, "vR2andes") == 0) {vR2andes = atof(tkn_v);}
        //fscanf(entra_var,"%lf",&vR2andes);

        else if (strcmp(tkn_w, "Kf") == 0) {Kf = atof(tkn_v);}
        //fscanf(entra_var, "%lf",&Kf);
        else if (strcmp(tkn_w, "Kd") == 0) {K_d = atof(tkn_v);}
        //fscanf(entra_var, "%lf",&K_d);
        else if (strcmp(tkn_w, "Km") == 0) {K_m = atof(tkn_v);}
        //fscanf(entra_var, "%lf",&K_m);

        else if (strcmp(tkn_w, "ls") == 0) {ls = atof(tkn_v);}
        //fscanf(entra_var, "%lf",&ls);
        else if (strcmp(tkn_w, "lb") == 0) {lb = atof(tkn_v);}
        //fscanf(entra_var, "%lf",&lb);
        else if (strcmp(tkn_w, "lb2") == 0) {lb2 = atof(tkn_v);}
        //fscanf(entra_var, "%lf",&lb2);

        else if (strcmp(tkn_w, "uplift_scale") == 0) {uplift_scale = atof(tkn_v);}
        //fscanf(entra_var, "%lf",&uplift_scale);

        else if (strcmp(tkn_w, "time_ofchangeu") == 0) {time_ofchangeu = atof(tkn_v);}
        //fscanf(entra_var, "%lf",&time_ofchangeu);
        else if (strcmp(tkn_w, "uplift_scale2") == 0) {uplift_scale2 = atof(tkn_v);}
        //fscanf(entra_var, "%lf",&uplift_scale2);
        else if (strcmp(tkn_w, "tempo_max") == 0) {tempo_max = atof(tkn_v);}
        //fscanf(entra_var, "%lf",&tempo_max);
        else if (strcmp(tkn_w, "dt") == 0) {dt = atof(tkn_v);}
        //fscanf(entra_var, "%lf",&dt);

        else if (strcmp(tkn_w, "n_sub_dt") == 0) {n_sub_dt = atoi(tkn_v);}

        // Else
        else
        {
            fprintf(stderr, "Error. Unrecognized keyword <%s> on line <%d> in the parameter file.\n", tkn_w, nline);
            exit(1);
        }
    }

    printf("Parameters reading: complete\n");

    fclose(entra_var);

	
}