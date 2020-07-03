#include<math.h>

/* Liste de voisin  */

void MakeNeighList(System *sys, H2O *h2o)
{
	float rl = sys->rc + 0.5;  // verlet radius : rl = rc + skin
	h2o->nbNeigh = 0;

	/* Initialise les liste des chaque atome */
	for(int i=0; i< sys->natoms; i++)
	{
		sys->neighList[i] = malloc(sizeof(chainList*));
		sys->neighList[i] = InitChainList();
		DelEleChainList(sys->neighList[i]); // enlève le premier element de la liste, initialisé a 0.
	}

	for(int i=0; i< sys->natoms; i++)
	{
		if(sys->type[i] == 2)
		{
			sys->r[i][0] = h2o->r[0][0];
			sys->r[i][1] = h2o->r[0][1];
			sys->r[i][2] = h2o->r[0][2];
		}
		for(int j=i+1; j< sys->natoms; j++)
		{
			if(sys->type[j] == 2)
			{
				sys->r[j][0] = h2o->r[0][0];
				sys->r[j][1] = h2o->r[0][1];
				sys->r[j][2] = h2o->r[0][2];
			}
			double *coordRel = malloc(sizeof(double)*3);
			CoordRel(coordRel, sys->r[i], sys->r[j]);  // coordonnées relative de l'atome j par rapport a l'atome i.
			CheckPBCNeighList(coordRel, sys->sizeBox); // inclue les PBC dans le calcule des coordonnées relatives.
//			printf("i : %f %f %f\n", sys->r[i][0], sys->r[i][1],sys->r[i][2]);
//			printf("j : %f %f %f\n", sys->r[j][0], sys->r[j][1],sys->r[j][2]);
//			printf("i : %d j : %d coordrel : %f %f %f\n", i, j, coordRel[0], coordRel[1], coordRel[2]);
			double distSqr = CalcDistSqr(coordRel); // distance au carré
			if(distSqr <= pow(rl,2))
			{
				InsertEleChainList(sys->neighList[i], j);
				InsertEleChainList(sys->neighList[j], i);
				if(i == sys->id_sub || j == sys->id_sub) // si on a la molécule d'eau
					h2o->nbNeigh++;  // incémente le n bde voins de H2O
			}
			free(coordRel);
		}
	}
}

void CheckPBCNeighList(double *coordRel, double sizeBox)
/* Si la distance entre atomes i et j est plus grande que la moitié
 * de la boite, déplace la coordonnée relative pour placer l'atome j 
 * dans une cellule adjacente.
 */
{
	if(fabs(*coordRel) > sizeBox/2)
	{
		if(*coordRel > 0)
			*coordRel -= sizeBox;
		else if(*coordRel < 0)
			*coordRel += sizeBox;
	}
	if(fabs(*(coordRel+1)) > sizeBox/2)
	{
		if(*(coordRel+1) > 0)
			*(coordRel+1) -= sizeBox;
		else if(*(coordRel+1) < 0)
			*(coordRel+1) += sizeBox;
	}
	if(fabs(*(coordRel+2)) > sizeBox/2)
	{
		if(*(coordRel+2) > 0)
			*(coordRel+2) -= sizeBox;
		else if(*(coordRel+2) < 0)
			*(coordRel+2) += sizeBox;
	}
}

void CoordRel(double *coordRel, double *ri, double *rj)
/* Retourne les coordonnées de la particule j
 * par rapport a la particule i
 */
{
	*coordRel = rj[0] - ri[0];
	*(coordRel+1) = rj[1] - ri[1];
	*(coordRel+2) = rj[2] - ri[2];
}



double CalcDistSqr(double *coordRel)
{
	double distSqr = pow(coordRel[0],2) + pow(coordRel[1],2) + pow(coordRel[2],2)  ;    // distance au carré 
	return distSqr;
}







