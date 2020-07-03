/* Intégrateur de velocity-Verlet et NOSQUISH */

void Integrate(System *sys, H2O *h2o);
void NoSquish_0_dt2(System *sys, H2O *h2o, double **s);
void NoSquish_dt2_dt(System *sys, H2O *h2o, double **s);
void Create_fquat_com(System *sys, H2O *h2o, double **s);
void NoSquish_free_rotor(System *sys, H2O *h2o);
void NoSquish_rotate(System *sys, H2O *h2o, int k, double small_dt);
double *NoSquish_permute(int order, double *quat);
void VelocityVerlet_0_dt2(System *sys, H2O *h2o);
void VelocityVerlet_dt2_dt(System *sys, H2O *h2o);

void Integrate(System *sys, H2O *h2o)
{
	double **s = BuildSMatrix(h2o);
	NoSquish_0_dt2(sys, h2o, s);
	VelocityVerlet_0_dt2(sys, h2o);
	ComputeForces(sys, h2o);
	s = BuildSMatrix(h2o);
	NoSquish_dt2_dt(sys, h2o, s);
	VelocityVerlet_dt2_dt(sys, h2o);
	NormalizeQuat(h2o->q);
	free2dArray(s,4);
}

void NoSquish_0_dt2(System *sys, H2O *h2o, double **s)
{
	Create_fquat_com(sys, h2o, s);
	for(int i=0; i<4; i++) // loop over quaternion component
		h2o->pq[i] += h2o->fq[i] * sys->dt * 0.5; // Half-kick angular momentum

	for(int i=0; i<3; i++) // loop over x, y, z
	{
		h2o->v[0][i] += 0.5*sys->dt*h2o->f[0][i]/h2o->m_h2o;  // Update c-o-m velocities
		h2o->r[0][i] += h2o->v[0][i]*sys->dt;  // Update c-o-m positions 
	}
	NoSquish_free_rotor(sys, h2o); // update orientation

//	 Update H2O orientation and compute H2O atoms velocities 
	double pos_o[3]  = {h2o->r[1][0],h2o->r[1][1],h2o->r[1][2]};
	double pos_h1[3] = {h2o->r[2][0],h2o->r[2][1],h2o->r[2][2]};
	double pos_h2[3] = {h2o->r[3][0],h2o->r[3][1],h2o->r[3][2]};
	PosH2O(h2o); // H2O atoms in SFF

	double dr_o[3]  = {h2o->r[1][0] - pos_o[0] , h2o->r[1][1] - pos_o[1] , h2o->r[1][2] - pos_o[2]};
	double dr_h1[3] = {h2o->r[2][0] - pos_h1[0], h2o->r[2][1] - pos_h1[1], h2o->r[2][2] - pos_h1[2]};
	double dr_h2[3] = {h2o->r[3][0] - pos_h2[0], h2o->r[3][1] - pos_h2[1], h2o->r[3][2] - pos_h2[2]};
	for(int i=0; i<3; i++) // loop over xyz, update velocities O, H1 et H2
	{
		h2o->v[1][i] = dr_o[i]  / sys->dt;
		h2o->v[2][i] = dr_h1[i] / sys->dt;
		h2o->v[3][i] = dr_h2[i] / sys->dt;
	}

}

void Create_fquat_com(System *sys, H2O *h2o, double **s)
{
	for(int i=0; i<4; i++)
	{
		h2o->fq[i] = 2*s[i][1]*h2o->torque[0] + 2*s[i][2]*h2o->torque[1] + 2*s[i][3]*h2o->torque[2]; 
	}
}

void NoSquish_dt2_dt(System *sys, H2O *h2o, double **s)
{
	Create_fquat_com(sys, h2o, s);
	for(int i=0; i<4; i++) // loop over quaternion component
		h2o->pq[i] += h2o->fq[i] * sys->dt * 0.5; // Half-kick angular momentum

	for(int i=0; i<3; i++) // loop over x, y, z
		h2o->v[0][i] += 0.5*sys->dt*h2o->f[0][i]/h2o->m_h2o;  // Update c-o-m velocities



	/* Angular velocity computation */
	double w_q[4] = {0,0,0,0}; // angular velocity quaternion
	for(int i=0; i<4; i++) // loop over lines of S matrix (transpose !)
		w_q[i] = 2*s[0][i]*h2o->pq[0] + 2*s[1][i]*h2o->pq[1] + 2*s[2][i]*h2o->pq[2] + 2*s[3][i]*h2o->pq[3]; // transpose of S
	w_q[1] /= h2o->I[0][0]; w_q[2] /= h2o->I[1][1]; w_q[3] /= h2o->I[2][2]; 
	h2o->w[0] = w_q[1]; h2o->w[1] = w_q[2]; h2o->w[2] = w_q[3]; 
	/*------------------------------*/
}

double **BuildSMatrix(H2O *h2o)
{
	double **s = alloc_2darray(4,4);
	s[0][0]=h2o->q->w;  s[0][1]=-h2o->q->vx; s[0][2]=-h2o->q->vy; s[0][3]=-h2o->q->vz; 
	s[1][0]=h2o->q->vx; s[1][1]=h2o->q->w;   s[1][2]=-h2o->q->vz; s[1][3]=h2o->q->vy;
	s[2][0]=h2o->q->vy; s[2][1]=h2o->q->vz;  s[2][2]=h2o->q->w;   s[2][3]=-h2o->q->vx;
	s[3][0]=h2o->q->vz; s[3][1]=-h2o->q->vy; s[3][2]=h2o->q->vx;  s[3][3]=h2o->q->w;
	return s;
}


void NoSquish_free_rotor(System *sys, H2O *h2o)
{
	int mRot = 10;
	double small_dt = sys->dt/mRot;
	for(int i=0; i<mRot; i++)
	{
		NoSquish_rotate(sys, h2o, 3, small_dt/2);
		NoSquish_rotate(sys, h2o, 2, small_dt/2);
		NoSquish_rotate(sys, h2o, 1, small_dt);
		NoSquish_rotate(sys, h2o, 2, small_dt/2);
		NoSquish_rotate(sys, h2o, 3, small_dt/2);
	}
}

void NoSquish_rotate(System *sys, H2O *h2o, int k, double small_dt)
{
	double q[4] = {h2o->q->w, h2o->q->vx, h2o->q->vy, h2o->q->vz};  // Orientation quaternion
	double q_up[4] = {0,0,0,0};  // Updated quaternion
	double p[4] = {h2o->pq[0], h2o->pq[1], h2o->pq[2], h2o->pq[3]}; // Momenta quaternion
	double p_up[4] = {0,0,0,0};  // Updated momenta quaternion
	double zeta_dt = 0;

	/* Permutation operators*/
	double *perm_q = NoSquish_permute(k, q); // permuted quaternion
	double *perm_p = NoSquish_permute(k, p); // permuted momenta quaternion
	/* zeta*dt */
	zeta_dt = small_dt*((p[0]*perm_q[0] + p[1]*perm_q[1] + p[2]*perm_q[2] + p[3]*perm_q[3]) / (4*h2o->I[k-1][k-1]));
	for(int i=0; i<4; i++) 	// loop over quaternion component
	{
		q_up[i] = cos(zeta_dt)*q[i] + sin(zeta_dt)*perm_q[i];
		p_up[i] = cos(zeta_dt)*p[i] + sin(zeta_dt)*perm_p[i];
		h2o->pq[i] = p_up[i];
	}

	h2o->q->w = q_up[0]; h2o->q->vx = q_up[1];  h2o->q->vy = q_up[2];  h2o->q->vz = q_up[3]; 

	free(perm_q);
	free(perm_p);
}

double *NoSquish_permute(int order, double *q)
/* Permutation operator acting on quaternion */
{
	double *p = malloc(sizeof(double)*4); // permuted quaternion

	if(order == 0)
		{p[0]=q[0]; p[1]=q[1]; p[2]=q[2]; p[3]=q[3];}
	else if(order == 1)
	    {p[0]=-q[1]; p[1]=q[0]; p[2]=q[3]; p[3]=-q[2];}
	else if(order == 2)
		{p[0]=-q[2]; p[1]=-q[3]; p[2]=q[0]; p[3]=q[1];}
	else if(order == 3)
		{p[0]=-q[3]; p[1]=q[2]; p[2]=-q[1]; p[3]=q[0];}

	return p;
}

void VelocityVerlet_dt2_dt(System *sys, H2O *h2o)
{
	for(int i=0; i< sys->natoms; i++) // loop over atoms
	{
		if(sys->type[i] == 1)
		{
			for(int j=0; j<3; j++) // loop over x, y, z
				sys->v[i][j] += 0.5*sys->dt*sys->f[i][j]/sys->m_gr; 	// Half-kick
		}
	}
}

void VelocityVerlet_0_dt2(System *sys, H2O *h2o)
{
	for(int i=0; i< sys->natoms; i++) // loop over atoms
	{
		if(sys->type[i] == 1)
		{
			for(int j=0; j<3; j++) // loop over x, y, z
			{
				// Halk-kick
				sys->v[i][j] += 0.5*sys->dt*sys->f[i][j]/sys->m_gr;
				// Drift
				sys->r[i][j] += sys->v[i][j]*sys->dt;
				// PBC
				if(sys->r[i][j] > sys->sizeBox)
					sys->r[i][j] -= sys->sizeBox;
				if(sys->r[i][j] < 0)
					sys->r[i][j] += sys->sizeBox;
			}
		}
	}
}
