void WriteData(char *dataFile, System *sys, H2O *h2o)
{
	FILE *fp = fopen(dataFile, "w");
	fprintf(fp, "LAMMPS Description\n\n");
	fprintf(fp, "    %d atoms\n", sys->natoms+2);
	fprintf(fp, "    2 bonds\n");
	fprintf(fp, "    1 angles\n");
	fprintf(fp, "    0 dihedrals\n");
	fprintf(fp, "    0 impropers\n\n");
	fprintf(fp, "    3 atom types\n");
	fprintf(fp, "    1 bond types\n");
	fprintf(fp, "    1 angle types\n");
	fprintf(fp, "    0 dihedral types\n");
	fprintf(fp, "    0 improper types\n\n");
	fprintf(fp, "    0.0 36.75 xlo xhi\n");
	fprintf(fp, "    0.0 36.75 ylo yhi\n");
	fprintf(fp, "    0.0 36.75 zlo zhi\n\n");
	fprintf(fp, "Masses\n\n");
	fprintf(fp, "    1 39.948\n");
	fprintf(fp, "    2 16\n");
	fprintf(fp, "    3 1.01\n\n");
	fprintf(fp, "Pair Coeffs\n\n");
	fprintf(fp, "    1 0.238 3.405\n");
	fprintf(fp, "    2 0.16 3.12\n");
	fprintf(fp, "    3 0.0 0.0\n\n");
	fprintf(fp, "Bond Coeffs\n\n");
	fprintf(fp, "    1 600 0.9576\n\n");
	fprintf(fp, "Angle Coeffs\n\n");
	fprintf(fp, "    1 75 104.51\n\n");
	fprintf(fp, "Atoms\n\n");

	int id_o = 0;

	for(int i=0; i<sys->natoms; i++)
	{
		if(sys->type[i] == 1)
			fprintf(fp, "%d 1 1 0.0 %f %f %f\n",i+1 , sys->r[i][0]+sys->sizeBox/2, sys->r[i][1]+sys->sizeBox/2, sys->r[i][2]+sys->sizeBox/2);
		else if (sys->type[i] == 2)
		{
			double dcm = (h2o->m_h1+h2o->m_h2)*h2o->d_oh*cos(h2o->a_oh/2) /h2o->m_h2o ;
			double o[3]  = {-dcm, 0., 0.};  // centre de masse a l'origine du MFF
			double h1[3] = {o[0] + h2o->d_oh*cos(h2o->a_oh/2), o[1] +  h2o->d_oh*sin(h2o->a_oh/2), o[2]};
			double h2[3] = {h1[0], -h1[1], o[2]};
			RotateVecQuat(o, h2o->q);
			RotateVecQuat(h1, h2o->q);
			RotateVecQuat(h2, h2o->q);
			o[0] += h2o->r[0][0]+sys->sizeBox/2;
			o[1] += h2o->r[0][1]+sys->sizeBox/2;
			o[2] += h2o->r[0][2]+sys->sizeBox/2;
			h1[0] += h2o->r[0][0]+sys->sizeBox/2;
			h1[1] += h2o->r[0][1]+sys->sizeBox/2;
			h1[2] += h2o->r[0][2]+sys->sizeBox/2;
			h2[0] += h2o->r[0][0]+sys->sizeBox/2;
			h2[1] += h2o->r[0][1]+sys->sizeBox/2;
			h2[2] += h2o->r[0][2]+sys->sizeBox/2;

			fprintf(fp, "%d 2 2 0.0 %f %f %f\n", i+1, o[0], o[1], o[2]);
			id_o = i+1;
			fprintf(fp, "%d 2 3 0.0 %f %f %f\n", sys->natoms+1, h1[0], h1[1], h1[2]);
			fprintf(fp, "%d 2 3 0.0 %f %f %f\n", sys->natoms+2, h2[0], h2[1], h2[2]);
		}
	}

	fprintf(fp, "\nBonds\n\n1 1 %d %d \n2 1 %d %d\n\n", id_o, sys->natoms+1, id_o, sys->natoms+2);
	fprintf(fp, "Angles\n\n1 1 %d %d %d\n\n", sys->natoms+1, id_o, sys->natoms+2);
	fclose(fp);
}

FILE *OpenFile(char *file)
{
	FILE *fp = fopen(file, "w");
	return fp;
}

void CloseFile(FILE *fp)
{
	fclose(fp);
}

void WriteSnapshot(FILE *fp, System *sys, H2O *h2o, int step)
/* Écrit un snapshot en mode écriture ("w") ou append ("a") */
{
//	FILE *fp = fopen(trajFile, mode);
	fprintf(fp,"ITEM: TIMESTEP\n%d\n", step);
	fprintf(fp,"ITEM: NUMBER OF ATOMS\n%d\n", sys->natoms+2);
//	fprintf(fp,"ITEM: NUMBER OF ATOMS\n%d\n", sys->natoms);
	fprintf(fp,"ITEM: BOX BOUNDS pp pp pp\n0 %f\n0 %f\n0 %f\n", sys->sizeBox, sys->sizeBox, sys->sizeBox);
	fprintf(fp,"ITEM: ATOMS id element xu yu zu vx vy vz\n");
	for(int i=0; i< sys->natoms; i++)
	{
		if(sys->type[i] == 1)
			fprintf(fp,"%d Ar %f %f %f %f %f %f\n", i+1, sys->r[i][0], sys->r[i][1], sys->r[i][2], sys->v[i][0], sys->v[i][1], sys->v[i][2]);
		if(sys->type[i] == 2)
			fprintf(fp,"%d O %f %f %f %f %f %f\n", i+1, h2o->r[1][0], h2o->r[1][1], h2o->r[1][2], h2o->v[1][0], h2o->v[1][1], h2o->v[1][2]);
	}
	fprintf(fp,"%d H %f %f %f %f %f %f\n", sys->natoms+1, h2o->r[2][0], h2o->r[2][1], h2o->r[2][2], h2o->v[2][0], h2o->v[2][1], h2o->v[2][2]);
	fprintf(fp,"%d H %f %f %f %f %f %f\n", sys->natoms+2, h2o->r[3][0], h2o->r[3][1], h2o->r[3][2], h2o->v[3][0], h2o->v[3][1], h2o->v[3][2]);
}

void WriteEng(FILE *fp, System *sys, H2O *h2o, int step)
{
//	FILE *fp = fopen(logFile, mode);
	if(step == 0)
		fprintf(fp, "step, T, etot, epot, ekin, epot_lj, epot_maka, ekin_gr, h2o_ekin, h2o_ekinTrans, h2o_ekinRot, pTot_x, pTot_y, pTot_z\n");
	fprintf(fp,"%d %f %f %f %f %f %f %f %f %f %f %f %f %f\n", step, sys->T, sys->etot, sys->epot, sys->ekin, sys->epot_lj, sys->epot_maka, sys->gr_ekin, h2o->ekin, h2o->ekin_trans, h2o->ekin_rot, sys->pTot[0], sys->pTot[1], sys->pTot[2]);
//	fclose(fp);
}

void WriteAngVel(FILE *fp, H2O *h2o, int step)
{
//	FILE *fp = fopen(logFile, mode);
	fprintf(fp,"%d %f %f %f\n", step, h2o->w[0], h2o->w[1], h2o->w[2]);
//	fclose(fp);
}

void WriteVelMFF(FILE *fp, H2O *h2o, int step)
{
//	FILE *fp = fopen(logFile, mode);
	double v[3] = {h2o->v[0][0], h2o->v[0][1], h2o->v[0][2]};
	SFF2MFF(h2o->q, v);

	fprintf(fp,"%d %f %f %f\n", step, v[0], v[1], v[2]);
//	fclose(fp);
}

void WriteQuat(FILE  *fp, H2O *h2o, int step)
{
//	FILE *fp = fopen(logFile, mode);
	double q[4] = {h2o->q->w, h2o->q->vx, h2o->q->vy, h2o->q->vz};

	fprintf(fp,"%d %f %f %f %f\n", step, q[0], q[1], q[2], q[3]);
//	fclose(fp);
}


void WriteRestart(char *resFile, System *sys, H2O *h2o)
{
	FILE *fp  = fopen(resFile, "w");
	for(int i=0; i< sys->natoms; i++)
	{
		if(sys->type[i] == 1)
			fprintf(fp, "%d %f %f %f %f %f %f\n", sys->type[i],
						   	sys->r[i][0], sys->r[i][1], sys->r[i][2],
						   	sys->v[i][0], sys->v[i][1], sys->v[i][2]);
		if(sys->type[i] == 2)
			fprintf(fp, "%d %f %f %f %f %f %f %f %f %f %f %f %f %f\n", sys->type[i],
						   	h2o->r[0][0], h2o->r[0][1], h2o->r[0][2],
						   	h2o->v[0][0], h2o->v[0][1], h2o->v[0][2],
						   	h2o->q->w, h2o->q->vx, h2o->q->vy, h2o->q->vz,
						   	h2o->w[0], h2o->w[1], h2o->w[2]);
	}
	fclose(fp);
}


