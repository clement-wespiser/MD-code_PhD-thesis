typedef struct System System;
typedef struct H2O H2O;
typedef struct Element Element;
typedef struct chainList chainList;
typedef struct Quaternion Quaternion;

struct System
{
	double sizeBox, rc, m_gr, dt, a, T;
	double **r, **v, **f, epsilon, sigma, epot_lj, epot_maka, epot, ekin, etot, pTot[3], gr_ekin;
	int natoms, nstep, *type, id_sub;  // type d'atomes et id de l'atome a substituer
	chainList **neighList;
};

struct H2O
{
	double **r, **v, **f; // position, vitesse et force sur O, H1 et H2
	double *torque; // 3 vecteurs du gradient
	double *w, **I, ekin, ekin_trans, ekin_rot, TRot; // vitesse angulaire, tenseur d'inertie, ... 
	Quaternion *q; // quaternion d'orientation
	float m_o, m_h1, m_h2, d_oh, a_oh, m_h2o; // paramètre géométrique des masse des atomes
	int nbNeigh ; // nb de voisins de H2O;
	double *fq;  // quaternion de force
	double *pq;  // quaternion de momentum
};

struct Element
{
	int nb;
	Element *next;
};

struct chainList
{
	Element *first;
	int nbEle;
};

struct Quaternion
{
	double w, vx, vy, vz; 
};

// init.h
void InitSystem(char *confFile, System *sys);
void ReadParams(char * confFile, float *param);
void InitCoords(int nbMaille, System *sys);
void InitVels(System *sys, H2O *h2o);
void ReadResFile(char *resFile, char *configFile, System *sys, H2O *h2o);
void FreeSystem(System *sys, H2O *h2o);

// h2o.h
void InitH2O(char *confFile, H2O *h2o, System *sys);
double Inertie(H2O *h2o, double *axe);
void InitAngVels(System *sys, H2O *h2o);
void PosH2O(H2O *h2o);
void ReadParamsH2O(char *confFile, float *param);

// neighbors.h
void MakeNeighList(System *sys, H2O *h2o);
double CalcDistSqr(double *coordRel);
void CoordRel(double *coordRel, double *ri, double *rj);
void CheckPBCNeighList(double *coordRel, double sizeBox);

// calculator.h
void ComputeForces(System *sys, H2O *h2o);
void EvalProps(System *sys, H2O *h2o);

// integrator.h
void VelocityVerlet(System *sys, H2O *h2o);
static double *Torque(H2O *h2o);
int Sign(double x);
double **H2OSpheCoordNeigh(System *sys, H2O *h2o);
double *AngMomH2O(double **spheCoord_t, double **spheCoord_dt, H2O *h2o, System *sys);
void Nosquish(System *sys, H2O *h2o);
double **BuildSMatrix(H2O *h2o);

// array.h
double **alloc_2darray(int nrows, int ncolumns);
double *alloc_1d_double_array(int ncolumns);
int **alloc_2d_int_array(int nrows, int ncolumns);
float *alloc_1d_float_array(int ncolumns);
void zero_2darray(double **array, int nrows, int ncolumns);
void zero_1darray(double *array, int ncolumns);
void free2dArray(double **a, int m);
int *alloc_1d_int_array(int ncolumns);
chainList *InitChainList();
void InsertEleChainList(chainList *list, int ele);
void DelEleChainList(chainList *list);

// quaternion.h
void RotateVecQuat(double *r, Quaternion *q);
void InitQuat(Quaternion *q, double w, double *v);
void InitQuatAngVel(Quaternion *q, double *v, int dt);
void MultQuat(Quaternion *q3, Quaternion *p, Quaternion *q);
void NormalizeQuat(Quaternion *q);
double *CrossProd(double *v1, double *v2);
double DotProd(double *v1, double *v2);
void Normalize(Quaternion *q);
void SFF2MFF(Quaternion *q, double *v); 
void MFF2SFF(Quaternion *q, double *v); 
void Transpose3x3Matrix(double **m2, double **m1);
void Transpose4x4Matrix(double **m2, double **m1);

// output.h
void WriteData(char *dataFile, System *sys, H2O *h2o);
void WriteQuat(FILE *fp, H2O *h2o, int step);
void WriteRestart(char *resFile, System *sys, H2O *h2o);
//void WriteSnapshot(char *trajFile, char *mode, System *sys, H2O *h2o, int step);
void WriteSnapshot(FILE *fp, System *sys, H2O *h2o, int step);
void WriteEng(FILE *fp, System *sys, H2O *h2o, int step);
void WriteAngVel(FILE *fp, H2O *h2o, int step);


#include "init.h"
#include "h2o.h"
#include "neighbors.h"
#include "calculator.h"
#include "integrator2.h"
#include "output.h"

