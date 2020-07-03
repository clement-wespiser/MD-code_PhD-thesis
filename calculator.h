#include"potential.h"

/* Routine de calcul de force et énegies. */

void ComputeForces(System *sys, H2O *h2o)
{
	zero_2darray(h2o->f, 4, 3);
	zero_1darray(h2o->torque, 3);
	zero_2darray(sys->f, sys->natoms, 3);
	double epot_maka = 0;
	double epot_lj   = 0;

	for(int i=0; i<sys->natoms; i++)
	{
		Element *j = sys->neighList[i]->first;
		while(j != NULL)
		{
			if(sys->type[i]==2 && sys->type[j->nb]==1) // H2O-Ar interaction
			{
				double *spheCoord = SpheCoord(h2o, sys->r[j->nb]);
				epot_maka += Potential(spheCoord[0], spheCoord[1], spheCoord[2])/1.196266e-6;

				double **f = GradSphe(h2o, sys->r[j->nb]);
				MFF2SFF(h2o->q, f[1]);
				for(int k=0; k<3; k++) // loop over xyz component 
				{
					h2o->f[0][k] += f[0][0]*f[1][k];
				}
				
				double t[3] = {0,0,0};  // instantanous torque
				t[0] = f[0][1] * f[3][0] / spheCoord[0] / (2*PI);
				t[1] = f[0][1] * f[3][1] / spheCoord[0] / (2*PI);
				t[2] = f[0][2] / (spheCoord[0]*sin(spheCoord[1]) * 2*PI);

				for(int k=0; k<3; k++) // loop over xyz
					h2o->torque[k] += t[k];

				free(spheCoord);
				free2dArray(f, 4);

			}

			else if(sys->type[i]==1 && sys->type[j->nb]==2) // Ar-H2O interaction
			{
				double *spheCoord = SpheCoord(h2o, sys->r[i]);
				double r = spheCoord[0]; double theta = spheCoord[1]; double phi = spheCoord[2];
				free(spheCoord);

				double **f = GradSphe(h2o, sys->r[i]);
				for(int k=1; k<4; k++) // loop over 3 unit vectors of gradient
					MFF2SFF(h2o->q, f[k]);

				for(int k=0; k<3; k++) // loop over xyz
					sys->f[i][k] -= f[0][0]*f[1][k] + f[0][1]*f[2][k]/r/(2*PI) + f[0][2]*f[3][k]/(r*sin(theta)*2*PI);
				
				free2dArray(f, 4);
			}

			else if(sys->type[i]==1 && sys->type[j->nb]==1)
			{
				double *coordRel = malloc(sizeof(double)*3);
				CoordRel(coordRel, sys->r[i], sys->r[j->nb]);
				CheckPBCNeighList(coordRel, sys->sizeBox);

				double r2 = pow(coordRel[0],2) + pow(coordRel[1],2) +  pow(coordRel[2],2);
				double term6 = pow(pow(sys->sigma,2)/r2, 3);
				epot_lj += term6*(term6 -1);
				double epsilon = sys->epsilon*1.196266e-6;
				double f = 48*epsilon/r2 *term6*(term6 - 0.5);
				for(int k=0; k<3; k++) // loop over xyz
					sys->f[i][k] -= f*coordRel[k];

				free(coordRel);
			}
			j = j->next;
		}
		free(j);
	}
	sys->epot_lj = epot_lj*2*sys->epsilon; // factorize the 4 facor of LJ potential and divide by 2 because counted ij and ji
	sys->epot_maka = epot_maka;
	sys->epot = sys->epot_lj + sys->epot_maka;
}

void EvalProps(System *sys, H2O *h2o)
{
	sys->ekin       = 0;
	sys->gr_ekin    = 0;
	h2o->ekin       = 0;
	h2o->ekin_trans = 0;
	h2o->ekin_rot   = 0;
	double sum_mv2  = 0;  // sum of the product m*v^2

	for(int i=0; i<3; i++) // loop over xyz
		sys->pTot[i] = 0;

	for(int i=0; i<sys->natoms; i++)
	{
		double vSqr = 0; // squared velocity
		if(sys->type[i] == 1)
		{
			vSqr    = pow(sys->v[i][0],2) + pow(sys->v[i][1],2) + pow(sys->v[i][3],2);
			sum_mv2 += vSqr*sys->m_gr;
			sys->gr_ekin += vSqr;
			for(int j=0; j<3; j++)
				sys->pTot[j] += sys->m_gr*sys->v[i][j];
		}
		else if(sys->type[i] == 2)
		{
			vSqr = pow(h2o->v[0][0],2) + pow(h2o->v[0][1],2)+ pow(h2o->v[0][2],2);
			h2o->ekin_trans = (0.5*h2o->m_h2o*vSqr) / 1.196266e-6;
			sum_mv2 += vSqr*h2o->m_h2o;
			for(int j=0; j<3; j++)
			{
			//	sys->pTot[j] += h2o->m_h2o*h2o->v[0][j];
				sys->pTot[j] += h2o->m_o*h2o->v[1][j] + h2o->m_h1*h2o->v[2][j] + h2o->m_h2*h2o->v[3][j] ;
			}
		}
	}
	// angular velocity of water computed in integrator
	h2o->ekin_rot   = (0.5*(h2o->I[0][0]*pow(h2o->w[0],2) + h2o->I[1][1]*pow(h2o->w[1],2) + h2o->I[2][2]*pow(h2o->w[2],2))) / 1.196266e-6;
	h2o->ekin = h2o->ekin_rot + h2o->ekin_trans;

	sys->gr_ekin *= 0.5*sys->m_gr /1.196266e-6;
	sys->ekin = h2o->ekin + sys->gr_ekin;
	sys->etot = sys->ekin + sys->epot;

//	sys->T = sum_mv2/8.32048e-7/(3*sys->natoms-3);  // Boltzmann cst in uma.A^2.fs^-2.K^⁻1
	sys->T = 3*sys->ekin/0.69/(3*sys->natoms-3);  // Boltzmann cst in cm^-1.K^⁻1
}


void EvalVelDist(char *logFile, System *sys, H2O *h2o, int step, double vMin, double vMax, double dv)
/* Évalue la distribution de vitesse */
{
	char *mode; // mode d'ouverture du fichier
	if(step == 0)
		mode = "w";
	else
		mode = "a";

	FILE *fp = fopen(logFile, mode);
	fprintf(fp, "step%d\n", step);

	int nbins = (vMax-vMin) / dv;
	double *nv = malloc(sizeof(double)*nbins); // array pour contenir le nb de particule ayant une vitesse entre v et v+dv

	double v = 0; // borne inférieur de l'intervalle
	double sum=0; // aire totale sous la courbe de distribution, pour normalisation
	for(int n=0; n<nbins; n++)
	{
		v = n*dv;
		nv[n] = 0;
		for(int i=0; i<sys->natoms; i++)
		{
			double normV = 0;
			if(sys->type[i] == 1)
				normV = sqrt(pow(sys->v[i][0],2) + pow(sys->v[i][1],2) + pow(sys->v[i][2],2));
			else
				normV = sqrt(pow(h2o->v[0][0],2) + pow(h2o->v[0][1],2) + pow(h2o->v[0][2],2));

			if(normV>v && normV<=(v+dv))
				nv[n]++;
		}
		sum += nv[n]*dv;
	}

	for(int n=0; n<nbins; n++) // normalisation et écriture
	{
		nv[n] /= sum;
		fprintf(fp, "    %f %f\n", (n+0.5)*dv, nv[n]);
	}

	fprintf(fp, "\n");
	fclose(fp);
	free(nv);
}





