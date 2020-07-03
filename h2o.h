#include"quaternions.h"

void ReadParamsH2O(char *confFile, float *param)
{
	FILE *fp;   
	fp = fopen(confFile,"r");
	if (fp == NULL) {
		fprintf(stderr, "Can't open input file %s\n", confFile);
		exit(1);
	}

	char buff[256];
	int linenum = 0;
	while ((fgets (buff, sizeof buff, fp)) != NULL)
	{
		linenum++;
		if (buff[0] == '\n' || buff[0] == '#')
			continue;

		char tmp1[256], tmp2[256], tmp3[256] ;
		if (sscanf(buff, "%s %s %s",tmp1, tmp2, tmp3) != 3)
		{
			fprintf(stderr, "Syntax error, line %d, config file\n", linenum);
			exit(1);
		}
		else if (strcmp(tmp1,"h2o_w") == 0) {param[0]= atof(tmp3);}
		else if (strcmp(tmp1,"h2o_vx") == 0) {param[1]= atof(tmp3);}
		else if (strcmp(tmp1,"h2o_vy") == 0) {param[2]= atof(tmp3);}
		else if (strcmp(tmp1,"h2o_vz") == 0) {param[3]= atof(tmp3);}
		else if (strcmp(tmp1,"h2o_m_o") == 0) {param[4]= atof(tmp3);}
		else if (strcmp(tmp1,"h2o_m_h1") == 0) {param[5]= atof(tmp3);}
		else if (strcmp(tmp1,"h2o_m_h2") == 0) {param[6]= atof(tmp3);}
		else if (strcmp(tmp1,"h2o_d_oh") == 0) {param[7]= atof(tmp3);}
		else if (strcmp(tmp1,"h2o_a_oh") == 0) {param[8]= atof(tmp3);}

		else if (strcmp(tmp1,"nbMaille") == 0) {continue;}
		else if (strcmp(tmp1,"T") == 0) {continue;}
		else if (strcmp(tmp1,"dt") == 0) {continue;}
		else if (strcmp(tmp1,"nstep") == 0)	{continue;}
		else if (strcmp(tmp1,"a") == 0)	{continue;}
		else if (strcmp(tmp1,"m_gr") == 0)	{continue;}
		else if (strcmp(tmp1,"rc") == 0)	{continue;}
		else if (strcmp(tmp1,"epsilon") == 0)	{continue;}
		else if (strcmp(tmp1,"sigma") == 0)	{continue;}

		else
		{
			fprintf(stderr,"Unknown parameter: %s, line %d",tmp1, linenum);
			exit(1);
		}
	}
	fclose(fp);
}

void InitH2O(char *confFile, H2O *h2o, System *sys)
{
	float *param = alloc_1d_float_array(9);
	ReadParamsH2O(confFile, param);

	// Initialisation des array
	h2o->q     = malloc(sizeof(double)*4); // quaternion
	h2o->r = alloc_2darray(4, 3); // position du centre de masse et des 3 atomes
	h2o->v = alloc_2darray(4, 3); // vitesse du centre de masse 
	h2o->f = alloc_2darray(4, 3); // force sur le centre de masse et les 3 atomes
	h2o->torque = malloc(sizeof(double)*3); // torque
	h2o->pq     = malloc(sizeof(double)*4); // quaternion de momentum
	h2o->fq     = malloc(sizeof(double)*4); // quaternion de force
	zero_2darray(h2o->r, 4, 3);
	zero_2darray(h2o->v, 4, 3);
	zero_2darray(h2o->f, 4, 3);
	zero_1darray(h2o->torque, 3);
	zero_1darray(h2o->fq, 4);
	zero_1darray(h2o->pq, 4);
	// Valeurs par défaut des paramètres
	double qw  = 1;
	double qvx = 0;
	double qvy = 0;
	double qvz = 0;
	h2o->w    = malloc(sizeof(double)*3); // vitesse angulaire
	h2o->w[0] = 0.0;
	h2o->w[1] = 0.0;
	h2o->w[2] = 0.0;
	h2o->m_o  = 16; // masse des atomes
	h2o->m_h1 = 1;
	h2o->m_h2 = 1;
	h2o->d_oh = 0.9576; // distance OH et angle de h2O
	h2o->a_oh = 1.824;
	h2o->r[0][0] = sys->r[sys->id_sub][0]; // position du centre de masse de H2O
	h2o->r[0][1] = sys->r[sys->id_sub][1];
	h2o->r[0][2] = sys->r[sys->id_sub][2];
//	h2o->r[0][0] = 10;
//	h2o->r[0][1] = 10;
//	h2o->r[0][2] = 10;
//	h2o->v[0][0] = 0;
//	h2o->v[0][1] = 0;
//	h2o->v[0][2] = 0;
	//////////////////////////////////////
	qw  = param[0];
	qvx = param[1];
	qvy = param[2];
	qvz = param[3];
	h2o->m_o  = param[4];
	h2o->m_h1 = param[5];
	h2o->m_h2 = param[6];
	h2o->m_h2o = h2o->m_o + h2o->m_h1 + h2o->m_h2;
	h2o->d_oh = param[7];
	h2o->a_oh = param[8];
	double vecQ[3] = {qvx, qvy, qvz}; // partie vectorielle du quaternion
	InitQuat(h2o->q, qw, vecQ);
	
	// Axe principaux d'inertie de H2O
	double a[3] = {0,1,0};
	double b[3] = {1,0,0};
	double c[3] = {0,0,1};
	h2o->I = alloc_2darray(3,3); // tenseur d'inertie de H2O
	zero_2darray(h2o->I, 3, 3);
	h2o->I[0][0] = Inertie(h2o, b);
	h2o->I[1][1] = Inertie(h2o, a);
	h2o->I[2][2] = Inertie(h2o, c);

	PosH2O(h2o); // place les atomes de h2o dans le SFF
	free(param);
//	InitAngVels(sys, h2o);
}


void InitAngVels(System *sys, H2O *h2o)
{
	double k = 8.32049*pow(10, -7); // cste de Boltzmann en uma.A^2.fs^-2
	h2o->w[0] = 2*PI*sqrt(k*sys->T/h2o->I[1][1]);
	h2o->w[2] = 2*PI*sqrt(k*sys->T/h2o->I[2][2]);
}

double Inertie(H2O *h2o, double *axe)
/* Calcul le moment d'inertie de h2o selon un axe donné
 * Cet axe se trouve dans le MFF.
 */
{
    // === Position des atomes de H2O dans le MFF ===
	double dcm = ((h2o->m_h2 + h2o->m_h2) * h2o->d_oh *cos(h2o->a_oh/2))/h2o->m_h2o;
    double o[3]  = {-dcm, 0., 0.};  // centre de masse a l'origine du MFF
    double h1[3] = {o[0] + h2o->d_oh*cos(h2o->a_oh/2), o[1] +  h2o->d_oh*sin(h2o->a_oh/2), o[2]};
    double h2[3] = {h1[0], -h1[1], o[2]};

	double i = h2o->m_o*(pow(o[0],2) + pow(o[1],2) + pow(o[2],2) - pow(DotProd(o, axe),2))
               +  h2o->m_h1*(pow(h1[0],2) + pow(h1[1],2) + pow(h1[2],2) - pow(DotProd(h1, axe),2))
               +  h2o->m_h2*(pow(h2[0],2) + pow(h2[1],2) + pow(h2[2],2) - pow(DotProd(h2, axe),2));

	return i;
}

void PosH2O(H2O *h2o)
/* Place les atomes de H2O dans le SFF
 */
{
	double dcm = ((h2o->m_h2 + h2o->m_h2) * h2o->d_oh *cos(h2o->a_oh/2))/h2o->m_h2o;
    double o[3]  = {-dcm, 0., 0.};  // centre de masse a l'origine du MFF
	double h1[3] = {o[0] + h2o->d_oh*cos(h2o->a_oh/2), o[1] +  h2o->d_oh*sin(h2o->a_oh/2), o[2]};
	double h2[3] = {h1[0], -h1[1], o[2]};
	MFF2SFF(h2o->q, o);
	MFF2SFF(h2o->q, h1);
	MFF2SFF(h2o->q, h2);
	h2o->r[1][0] = o[0]  + h2o->r[0][0];
	h2o->r[1][1] = o[1]  + h2o->r[0][1];
	h2o->r[1][2] = o[2]  + h2o->r[0][2];
	h2o->r[2][0] = h1[0] + h2o->r[0][0];
	h2o->r[2][1] = h1[1] + h2o->r[0][1];
	h2o->r[2][2] = h1[2] + h2o->r[0][2];
	h2o->r[3][0] = h2[0] + h2o->r[0][0];
	h2o->r[3][1] = h2[1] + h2o->r[0][1];
	h2o->r[3][2] = h2[2] + h2o->r[0][2];
	
}
