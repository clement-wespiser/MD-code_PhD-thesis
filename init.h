/* Routine d'initialisation du système
 * Retourne un pointeur vers un array de positions, vitesses, etc
 */

#include <string.h>
#include <math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include "array.h"

void InitSystem(char *confFile, System *sys)
{
	float *param = alloc_1d_float_array(9); // array contenant param de maille, nbMaille, T.
	ReadParams(confFile, param);
	int nbMaille = param[1];
	sys->T = param[2];
    float m_gr = param[3], paramMaille = param[0];
	
	sys->epsilon = param[8];
	sys->sigma = param[7];
	sys->rc = param[6];
	sys->sizeBox = nbMaille*paramMaille;  // taille d'un côté de la cellule de simulation
	sys->natoms  =  pow(nbMaille,3)*4; // 4 atomes par maille + les 2 hydrogènes
//	sys->id_sub = 0;
//	sys->natoms = 2;

	sys->type = alloc_1d_int_array(sys->natoms);
	for(int i=0; i<sys->natoms; i++)
		sys->type[i] = 1;

	sys->m_gr  = param[3];
	sys->dt    = param[4];
	sys->nstep = param[5];
	sys->a     = param[0]; // paramètres de maille

	sys->r = alloc_2darray(sys->natoms, 3); // array des positions
	InitCoords(nbMaille, sys);

//---------------- Comment this line to NOT put water in the system
//---------------- Adapt output.h for snapshot dumping
	sys->type[sys->id_sub] = 2;
//-----------------------------------------------------------------

//	sys->r[1][0] = 14;
//	sys->r[1][1] = 10;
//	sys->r[1][2] = 10;
//	sys->r[0][0] = 18;
//	sys->r[0][1] = 10;
//	sys->r[0][2] = 10;
//	sys->r[2][0] = 6;
//	sys->r[2][1] = 10;
//	sys->r[2][2] = 10;

	sys->v = alloc_2darray(sys->natoms, 3); // array des vitesses 
	sys->f = alloc_2darray(sys->natoms, 3); // array des forces
	sys->neighList = malloc(sizeof(int*)*sys->natoms);

	free(param);
}

void ReadParams(char *confFile, float *param)
{
	FILE *fp;   
	fp = fopen(confFile,"r");
	if (fp == NULL) {
		fprintf(stderr, "Can't open input file %s\n", confFile);
		exit(1);
	}

	// Valeurs par défaut des paramètres
	float nbMaille = 4;
	float T = 10;
	double  dt = 1;
	int nstep = 500 ; 
	float a = 5.25 ;   // parametre de maille
	float m_gr = 40.0 ;   // masse des atomes de gaz rare
	float rc = 10;  // cutoff radius
	double epsilon = 87.37; // unité : cm^-1
	double sigma = 3.345;
	//////////////////////////////////////

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
		else if (strcmp(tmp1,"nbMaille") == 0) {nbMaille = atof(tmp3);}
		else if (strcmp(tmp1,"T") == 0) {T = atof(tmp3);}
		else if (strcmp(tmp1,"dt") == 0) {dt = atof(tmp3);}
		else if (strcmp(tmp1,"nstep") == 0)	{nstep = atoi(tmp3);}
		else if (strcmp(tmp1,"a") == 0)	{a = atof(tmp3);}
		else if (strcmp(tmp1,"m_gr") == 0)	{m_gr = atof(tmp3);}
		else if (strcmp(tmp1,"rc") == 0)	{rc = atof(tmp3);}
		else if (strcmp(tmp1,"epsilon") == 0)	{epsilon = atof(tmp3);}
		else if (strcmp(tmp1,"sigma") == 0)	{sigma= atof(tmp3);}
		else if (strcmp(tmp1,"h2o_w") == 0) {continue;}
		else if (strcmp(tmp1,"h2o_vx") == 0) {continue;}
		else if (strcmp(tmp1,"h2o_vy") == 0) {continue;}
		else if (strcmp(tmp1,"h2o_vz") == 0) {continue;}
		else if (strcmp(tmp1,"h2o_m_o") == 0) {continue;}
		else if (strcmp(tmp1,"h2o_m_h1") == 0) {continue;}
		else if (strcmp(tmp1,"h2o_m_h2") == 0) {continue;}
		else if (strcmp(tmp1,"h2o_d_oh") == 0) {continue;}
		else if (strcmp(tmp1,"h2o_a_oh") == 0) {continue;}
		else
		{
			fprintf(stderr,"Unknown parameter: %s, line %d",tmp1, linenum);
			exit(1);
		}
	}
	fclose(fp);

	*param = a;
	*(param+1) = nbMaille;
	*(param+2) = T;
	*(param+3) = m_gr;
	*(param+4) = dt;
	printf("dt : %f\n", dt);
	*(param+5) = nstep;
	*(param+6) = rc;
	*(param+7) = sigma;
	*(param+8) = epsilon;

}

void InitCoords(int nbMaille, System *sys)
/* Maille CFC */
{
	int id_atom=0;
	float a = sys->a;

	int mailleSub = 0; // id de la maille a substituer 
	if (nbMaille % 2 != 0)
		mailleSub = (nbMaille-1)/2;
	else
		mailleSub = nbMaille/2;

	for(int i=0; i<nbMaille; i++)
	{
		for(int j=0; j<nbMaille; j++)
		{
			for(int k=0; k<nbMaille; k++)
			{
				if(i == mailleSub && j == mailleSub && k ==mailleSub)
				   	sys->id_sub = id_atom+1;

				sys->r[id_atom][0] = 0+i*a;
				sys->r[id_atom][1] = 0+j*a;
				sys->r[id_atom][2] = 0+k*a;
				sys->r[id_atom+1][0] = a/2+i*a;
				sys->r[id_atom+1][1] = 0+j*a;
				sys->r[id_atom+1][2] = a/2+k*a;
				sys->r[id_atom+2][0] = a/2+i*a;
				sys->r[id_atom+2][1] = a/2+j*a;
				sys->r[id_atom+2][2] = 0+k*a;
				sys->r[id_atom+3][0] = 0+i*a;
				sys->r[id_atom+3][1] = a/2+j*a;
				sys->r[id_atom+3][2] = a/2+k*a;
				id_atom = id_atom+4;
			}
		}
	}
}

void InitVels(System *sys, H2O *h2o)
/* Distribution gaussienne des composantes 
 * de la vitesse pour chaque atome.
 */
{
	float T_eng = sys->T*8.32049*pow(10,-7) ;   // température en unité d'énergie (uma.A^2.fs^-2)
	gsl_rng *rng = gsl_rng_alloc(gsl_rng_taus);  // gaussian random number generator
	double v_moy[3] = {0,0,0};
	double sigma = 0; // écart-type de la distribution gaussienne pour les vitesses  

	for(int i = 0; i< sys->natoms; i++)
	{
		if(sys->type[i] == 2)
		{
			sigma = sqrt(T_eng/h2o->m_h2o); 
			h2o->v[0][0] = gsl_ran_gaussian(rng, sigma); // composantes x des vitesses
			h2o->v[0][1] = gsl_ran_gaussian(rng, sigma); // composantes y des vitesses
			h2o->v[0][2] = gsl_ran_gaussian(rng, sigma); // composantes z des vitesses
			v_moy[0] += h2o->v[0][0];
			v_moy[1] += h2o->v[0][1];
			v_moy[2] += h2o->v[0][2];
		}
		else
		{
			sigma = sqrt(T_eng/sys->m_gr);
			sys->v[i][0] = gsl_ran_gaussian(rng, sigma); // composantes x des vitesses
			sys->v[i][1] = gsl_ran_gaussian(rng, sigma); // composantes y des vitesses
			sys->v[i][2] = gsl_ran_gaussian(rng, sigma); // composantes z des vitesses
			v_moy[0] += sys->v[i][0];
			v_moy[1] += sys->v[i][1];
			v_moy[2] += sys->v[i][2];
		}
	}

	v_moy[0] /= sys->natoms;
	v_moy[1] /= sys->natoms;
	v_moy[2] /= sys->natoms;

	for(int i = 0; i< sys->natoms; i++)
	{
		if(sys->type[i] == 2)
		{
			h2o->v[0][0] -= v_moy[0];
			h2o->v[0][1] -= v_moy[1];
			h2o->v[0][2] -= v_moy[2];
		}
		else
		{
			sys->v[i][0] -= v_moy[0];
			sys->v[i][1] -= v_moy[1];
			sys->v[i][2] -= v_moy[2];
		}
	}
	gsl_rng_free(rng); 
}

void FreeSystem(System *sys, H2O *h2o)
{
	free(h2o->q);
	free(h2o->r);
	free(h2o->v);
	free(h2o->w);
	free(h2o->torque);
	free2dArray(h2o->I, 3);

	free(sys->type);
	free2dArray(sys->r, 3);
	free2dArray(sys->f, 3);
	free2dArray(sys->v, 3);

	for(int i=0; i<sys->natoms; i++)
		free(sys->neighList[i]);
	free(sys->neighList);
}

void ReadResFile(char *resFile, char *confFile, System *sys, H2O *h2o)
{
	/* Initialisation gaz rare */
	float *param = alloc_1d_float_array(9); // array contenant param de maille, nbMaille, T.
	ReadParams(confFile, param);
	int nbMaille = param[1];
	sys->T = param[2];
    float m_gr = param[3], paramMaille = param[0];
	
	sys->epsilon = param[8];
	sys->sigma = param[7];
	sys->rc = param[6];
	sys->sizeBox = nbMaille*paramMaille;  // taille d'un côté de la cellule de simulation
	sys->natoms  =  pow(nbMaille,3)*4; // 4 atomes par maille + les 2 hydrogènes

	sys->type = alloc_1d_int_array(sys->natoms);
	sys->m_gr = param[3];
	sys->dt = param[4];
	sys->nstep = param[5];
	sys->a = param[0]; // paramètres de maille

	sys->r = alloc_2darray(sys->natoms, 3); // array des positions
	sys->v = alloc_2darray(sys->natoms, 3); // array des vitesses 
	sys->f = alloc_2darray(sys->natoms, 3); // array des forces
	sys->neighList = malloc(sizeof(int*)*sys->natoms);

	free(param);
	/*-------------------------*/

	/* Initialisation H2O*/
	float *param_h2o = alloc_1d_float_array(9); // array contenant param de maille, nbMaille, T.
	ReadParamsH2O(confFile, param_h2o);
	h2o->r = alloc_2darray(4,3);
	h2o->v = alloc_2darray(3,3);
	h2o->w = malloc(sizeof(double)*3); // vitesse angulaire
	h2o->q = malloc(sizeof(double)*4);
	h2o->f = alloc_2darray(4,3);
	h2o->torque = malloc(sizeof(double)*3);
	// Axe principaux d'inertie de H2O
	double a[3] = {0,1,0};
	double b[3] = {1,0,0};
	double c[3] = {0,0,1};
	h2o->I = alloc_2darray(3,3); // tenseur d'inertie de H2O
	zero_2darray(h2o->I, 3, 3);
	h2o->I[0][0] = Inertie(h2o, b);
	h2o->I[1][1] = Inertie(h2o, a);
	h2o->I[2][2] = Inertie(h2o, c);
	h2o->d_oh = param[7];
	h2o->a_oh = param[8];
	h2o->m_o  = param_h2o[4];
	h2o->m_h1 = param_h2o[5];
	h2o->m_h2 = param_h2o[6];
	h2o->m_h2o = h2o->m_o + h2o->m_h1 + h2o->m_h2;

	// Lecture du fichier restart
	char buff[256];
	int linenum = 0;
	FILE *fp = fopen(resFile, "r");
	while ((fgets (buff, sizeof buff, fp)) != NULL)
	{
		if (buff[0] == '\n' || buff[0] == '#')
			continue;

		char tmp1[256], tmp2[256],tmp3[256],tmp4[256],tmp5[256],tmp6[256],tmp7[256], tmp8[256], tmp9[256], tmp10[256], tmp11[256], tmp12[256], tmp13[256], tmp14[256]; // type, position et vitesses
		sscanf(buff, "%s %s %s %s %s %s %s %s %s %s %s %s %s %s",tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8,  tmp9, tmp10,tmp11, tmp12, tmp13, tmp14);

		sys->type[linenum] = atoi(tmp1);
		if(sys->type[linenum] == 1)
		{
			sys->r[linenum][0] = atof(tmp2);
			sys->r[linenum][1] = atof(tmp3);
			sys->r[linenum][2] = atof(tmp4);
			sys->v[linenum][0] = atof(tmp5);
			sys->v[linenum][1] = atof(tmp6);
			sys->v[linenum][2] = atof(tmp7);
		}
		else if(sys->type[linenum] == 2)
		{
			h2o->r[0][0] = atof(tmp2);
			h2o->r[0][1] = atof(tmp3);
			h2o->r[0][2] = atof(tmp4);
			h2o->v[0][0] = atof(tmp5);
			h2o->v[0][1] = atof(tmp6);
			h2o->v[0][2] = atof(tmp7);
			h2o->q->w    = atof(tmp8);
			h2o->q->vx   = atof(tmp9);
			h2o->q->vy   = atof(tmp10);
			h2o->q->vz   = atof(tmp11);
			h2o->w[0]    = atof(tmp12);
			h2o->w[1]    = atof(tmp13);
			h2o->w[2]    = atof(tmp14);
		}

		linenum++;
	}

	fclose(fp);
	free(param_h2o);
	PosH2O(h2o);
}
