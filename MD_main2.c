#define PI 3.14159265358979323846
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "prototypes.h"

int main (int argc, char **argv)
{
	argv[1] = "config_test";
	argv[2] = "traj_test.lammpstrj" ;
	argv[3] = "log_test_dt2_p.txt";
	argv[4] = "log_w_test.txt" ;
	argv[7] = "log_v-mff_test.txt" ;
	argv[5] = "velDist_test.txt";
	argv[6] = "resFile_test.txt";
	argv[8] = "log_q_test.txt" ;

	int restart = 0;

	System sys;  // structure System contient la taille de la cellule, les positions, vitesses, etc..
	H2O h2o;

	if(restart == 0)
	{
    	InitSystem(argv[1], &sys);  // initialise positions des atomes
		InitH2O(argv[1], &h2o, &sys);
		InitVels(&sys, &h2o); // initialise vitesse des atomes
	}
	else	
		ReadResFile(argv[6], argv[1], &sys, &h2o);

	MakeNeighList(&sys, &h2o); // Liste ds voisins

	/* Ouverture des fichiers d'output */
	FILE *snapshot_fp = OpenFile(argv[2]);
	FILE *eng_fp      = OpenFile(argv[3]);
	FILE *angVel_fp   = OpenFile(argv[4]);
	FILE *velMFF_fp   = OpenFile(argv[7]);
	FILE *quat_fp     = OpenFile(argv[8]);
	/*---------------------------------*/

	ComputeForces(&sys, &h2o); // Calcul énergie potentielle et forces
	EvalProps(&sys, &h2o); // Évalue les proprités (kinEng,etc )

	WriteEng(eng_fp, &sys, &h2o, 0);      // écriture des énergies
//	WriteSnapshot(snapshot_fp, &sys, &h2o, 0); // écritre trajectoire
//	WriteAngVel(angVel_fp, &h2o, 0);         // écriture vitesses angulaires
//	WriteVelMFF(velMFF_fp, &h2o, 0);       // écriture vitesse MFF
//	WriteQuat(quat_fp, &h2o, 0);         // ecriture quaternion d'orientation
//	WriteData("data", &sys, &h2o);
//

	for(int i=0; i<sys.nstep; i++)
	{
		fprintf(stdout,"step %d ", i);
		Integrate(&sys, &h2o);
		EvalProps(&sys, &h2o);
//		if(i%1000 == 0)
//			EvalVelDist(argv[5], &sys, &h2o, i, 0, 0.003, 0.00001);
//		WriteSnapshot(snapshot_fp, &sys, &h2o, i+1); // écritre trajectoire
		WriteEng(eng_fp, &sys, &h2o, i+1);
//		WriteAngVel(angVel_fp, &h2o, i+1); 
//		WriteVelMFF(velMFF_fp, &h2o, i+1); 
//		WriteQuat(quat_fp, &h2o, i+1);
	}
// Écriture des types, positions et vitesses pour chaque atomes
// dans un fichier, pour repartir une simulation à partir de ce point
//	WriteRestart(argv[6], &sys, &h2o);

	/* Closing output files */	
	CloseFile(snapshot_fp);
	CloseFile(eng_fp);
	CloseFile(angVel_fp);
	CloseFile(velMFF_fp);
	CloseFile(quat_fp);
	/*----------------------*/	


//	FreeSystem(&sys, &h2o);

	return 0;
}
