void InitQuat(Quaternion *q, double w, double *v)
/* Initialise un quaternion d'orientarion.
 * Seuls les quaternions représentant des orientations 
 * doivent être normalisés
 */
{
	double norm_v = sqrt(pow(v[0],2) + pow(v[1],2) + pow(v[2],2));
	if(norm_v == 0)
	{
		q->vx = 0;
		q->vy = 0;
		q->vz = 0;
		q->w = 1;
	}
	else
	{
		q->w = cos(w/2);
		q->vx = sin(w/2)*v[0];// /norm_v;
		q->vy = sin(w/2)*v[1];// /norm_v;
		q->vz = sin(w/2)*v[2];// /norm_v;
	}
	NormalizeQuat(q);
}

void InitQuatAngVel(Quaternion *q,  double *v, int dt)
{
	q->w = 0;
	q->vx = v[0]*dt;
	q->vy = v[1]*dt;
	q->vz = v[2]*dt;
}

double *CrossProd(double *v1, double *v2)
/* Retourne le produit vectoriel entr les vecteurs v1 et v2
 */
{
	double v1x = *v1;
	double v1y = *(v1+1);
	double v1z = *(v1+2);
	double v2x = *v2;
	double v2y = *(v2+1);
	double v2z = *(v2+2);

	double *crossProd = malloc(sizeof(double)*3);
//	static double crossProd[3] = {0,0,0};

	 *crossProd     = v1y*v2z - v1z*v2y; 
	 *(crossProd+1) = v1z*v2x - v1x*v2z;
	 *(crossProd+2) = v1x*v2y - v1y*v2x;

	 return crossProd;
}

double DotProd(double *v1, double *v2)
{
	double dotProd = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
	return dotProd;
}

void RotateVecQuat(double *r, Quaternion *q)
/* Rotate le vecteur r selon le quaternion q.
    Paramètres : r : vecteur a rotater
                 q : quaternion
*/
{
	double v[3] = {q->vx, q->vy, q->vz}; // partie vectorielle du quaternion
//	printf("r ini : %.15f %.15f %.15f\n", r[0], r[1], r[2]);
//	printf("norme : %.15f\n", sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]));
//	printf("v: %f %f %f\n\n", v[0], v[1], v[2]);
	double r_ini[3] = {r[0], r[1], r[2]};

	// Cross Products
	double *cp1 = CrossProd(v,r_ini);
	double *cp2 = CrossProd(v,cp1);

	for(int i=0; i<3; i++)
		r[i] += 2*q->w*cp1[i] + 2*cp2[i];
//	printf("r fini : %.15f %.15f %.15f\n", r[0], r[1], r[2]);
//	printf("norme : %.15f\n", sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]));

	free(cp1);
	free(cp2);
}

void MultQuat(Quaternion *q3, Quaternion *p, Quaternion *q)
/*    Multiplie deux quaternions.
    La multiplication de 2 quaternions n'est pas commutative : p*q != q*p
    Paramètres : q3 : résultat
				 p : premier quaternion
                 q : deuxième quaternion
    Retourne : le produit des deux quaternions q3 = p*q != q*p
	*/
{
	q3->w  = p->w*q->w - p->vx*q->vx - p->vy*q->vy - p->vz*q->vz; 
	q3->vx = p->w*q->vx + p->vx*q->w + p->vy*q->vz - p->vz*q->vy;
	q3->vy = p->w*q->vy - p->vx*q->vz +p->vy*q->w + p->vz*q->vx;
	q3->vz = p->w*q->vz + p->vx*q->vy - p->vy*q->vx + p->vz*q->w;
}

double **InvMat(double **m)
/* Inversion de la matrice m
 */
{
	double m_tmp[3][3];
	for(int i=0; i<3; i++)
	{
		for(int j =0; j<3; j++)
			m_tmp[i][j] = m[i][j];
	}

	double **m3 = alloc_2darray(3,3);

	double t1 = m_tmp[0][0]*m_tmp[1][1];
	double t2 = m_tmp[0][0]*m_tmp[1][2];
	double t3 = m_tmp[0][1]*m_tmp[1][0];
	double t4 = m_tmp[0][2]*m_tmp[1][0];
	double t5 = m_tmp[0][1]*m_tmp[2][1];
	double t6 = m_tmp[0][2]*m_tmp[2][1];

	double invDet = 1/(t1*m_tmp[2][2] - t2*m_tmp[2][1] - t3*m_tmp[2][2] + t4*m_tmp[2][1] + t5*m_tmp[1][2] - t6*m_tmp[1][1]);

	m3[0][0] = (m_tmp[1][1]*m_tmp[2][2] - m_tmp[1][2]*m_tmp[2][1]) * invDet;
	m3[0][1] = -(m_tmp[0][1]*m_tmp[2][2] - m_tmp[0][2]*m_tmp[2][1])*invDet;
	m3[0][2] = (m_tmp[0][1]*m_tmp[1][2] - m_tmp[0][2]*m_tmp[1][1])*invDet;
	m3[1][0] = -(m_tmp[1][0]*m_tmp[2][2] - m_tmp[1][2]*m_tmp[2][0])*invDet;
	m3[1][1] = (m_tmp[0][0]*m_tmp[2][2] - t6) * invDet;
	m3[1][2] = -(t2-t4)*invDet;
	m3[2][0] = (m_tmp[1][0]*m_tmp[2][1] - m_tmp[1][1]*m_tmp[2][0]) * invDet;
	m3[2][1] = -(m_tmp[0][0]*m_tmp[2][1] - t5) * invDet;
	m3[2][2] = (t1-t3) * invDet;

	return m3;
}

void NormalizeQuat(Quaternion *q)
{
	double norm = sqrt(pow(q->w,2) + pow(q->vx,2) + pow(q->vy,2) + pow(q->vz,2));
	q->w /= norm;
	q->vx /= norm;
	q->vy /= norm;
	q->vz /= norm;
}

void Quat2RotMat(double **m, Quaternion *q)
/* Forme une matrice de rotation à partir des composantes
 * du quaternion q
 */
{
	double q0 = q->w;
	double q1 = q->vx;
	double q2 = q->vy;
	double q3 = q->vz;

	if(fabs(q0) < 1e-6) q0=0;
	if(fabs(q1) < 1e-6) q1=0;
	if(fabs(q2) < 1e-6) q2=0;
	if(fabs(q3) < 1e-6) q3=0;

//	printf("q (q2mat) : %.10f %.10f %.10f %.10f\n", q0, q1, q2, q3);

	m[0][0] = q0*q0 +q1*q1 -q2*q2 -q3*q3;
	m[0][1] = 2*(q1*q2 + q0*q3);
	m[0][2] = 2*(q1*q3 - q0*q2);
	m[1][0] = 2*(q1*q2-q0*q3);
	m[1][1] = q0*q0 -q1*q1 +q2*q2 -q3*q3;
	m[1][2] = 2*(q2*q3 + q0*q1);
	m[2][0] = 2*(q1*q3 + q0*q2);
	m[2][1] = 2*(q2*q3 - q0*q1);
	m[2][2] = q0*q0 -q1*q1 -q2*q2 +q3*q3;
}

double **MultMat(double **m1, double **m2)
/* Multiplie les matrice m1 et m2
 */
{
	double **m3 = alloc_2darray(3,3);
	m3[0][0] = m1[0][0]*m2[0][0] + m1[0][1]*m2[1][0] + m1[0][2]*m2[2][0];
	m3[0][1] = m1[0][0]*m2[0][1] + m1[0][1]*m2[1][1] + m1[0][2]*m2[2][1];
	m3[0][2] = m1[0][0]*m2[0][2] + m1[0][1]*m2[1][2] + m1[0][2]*m2[2][2];
	m3[1][0] = m1[1][0]*m2[0][0] + m1[1][1]*m2[1][0] + m1[1][2]*m2[2][0];
	m3[1][1] = m1[1][0]*m2[0][1] + m1[1][1]*m2[1][1] + m1[1][2]*m2[2][1];
	m3[1][2] = m1[1][0]*m2[0][2] + m1[1][1]*m2[1][2] + m1[1][2]*m2[2][2];
	m3[2][0] = m1[2][0]*m2[0][0] + m1[2][1]*m2[1][0] + m1[2][2]*m2[2][0];
	m3[2][1] = m1[2][0]*m2[0][1] + m1[2][1]*m2[1][1] + m1[2][2]*m2[2][1];
	m3[2][2] = m1[2][0]*m2[0][2] + m1[2][1]*m2[1][2] + m1[2][2]*m2[2][2];

	return m3;
}

void MultMatVec(double *v2, double *v, double **m)
/* Multiplie le vecteur v par la matrice m.
 * Stocke le résultats dans v2.
 */
{
	v2[0] = m[0][0]*v[0] + m[0][1]*v[1] + m[0][2]*v[2];
	v2[1] = m[1][0]*v[0] + m[1][1]*v[1] + m[1][2]*v[2];
	v2[2] = m[2][0]*v[0] + m[2][1]*v[1] + m[2][2]*v[2];
}

void SFF2MFF(Quaternion *q, double *v)
/* Transforme the v vector from SFF to MFF
 * using the quaternion q. q is transformed 
 * to the rotation matrix first
 */
{
	double **m = alloc_2darray(3,3);
	zero_2darray(m, 3,3);
	Quat2RotMat(m, q);
//	printf("m %.10f %.10f %.10f\n", m[0][0], m[0][1], m[0][2]);
//	printf("m %.10f %.10f %.10f\n", m[1][0], m[1][1], m[1][2]);
//	printf("m %.10f %.10f %.10f\n", m[2][0], m[2][1], m[2][2]);
	double v2[3] = {0,0,0};
	MultMatVec(v2, v, m); 
	v[0] = v2[0]; v[1]=v2[1]; v[2]=v2[2];
	free2dArray(m, 3);
}

void MFF2SFF(Quaternion *q, double *v)
/* Transforme the v vector from SFF to MFF
 * using the quaternion q. q is transformed 
 * to the rotation matrix first
 */
{
	double **m  = alloc_2darray(3,3);
	double **m2 = alloc_2darray(3,3); // transpose of m
	Quat2RotMat(m, q);
	Transpose3x3Matrix(m2, m);
	double v2[3] = {0,0,0};
	MultMatVec(v2, v, m2); 
	v[0] = v2[0]; v[1]=v2[1]; v[2]=v2[2];
	free2dArray(m, 3);
	free2dArray(m2, 3);
}

void Transpose3x3Matrix(double **m2, double **m1)
/* Compute the transpose of matrix m1
 * and store the result in matrix m2.
 * Work only for 3x3 matrices.
 */
{
	m2[0][0] = m1[0][0]; m2[0][1] = m1[1][0]; m2[0][2] = m1[2][0]; 
	m2[1][0] = m1[0][1]; m2[1][1] = m1[1][1]; m2[1][2] = m1[2][1]; 
	m2[2][0] = m1[0][2]; m2[2][1] = m1[1][2]; m2[2][2] = m1[2][2]; 
}

void Transpose4x4Matrix(double **m2, double **m1)
/* Compute the transpose of matrix m1
 * and store the result in matrix m2.
 * Work only for 4x4 matrices.
 */
{
	m2[0][0] = m1[0][0]; m2[0][1] = m1[1][0]; m2[0][2] = m1[2][0]; m2[0][3] = m1[3][0];
	m2[1][0] = m1[0][1]; m2[1][1] = m1[1][1]; m2[1][2] = m1[2][1]; m2[1][3] = m1[3][1];
	m2[2][0] = m1[0][2]; m2[2][1] = m1[1][2]; m2[2][2] = m1[2][2]; m2[2][3] = m1[3][2];  
	m2[3][0] = m1[0][3]; m2[3][1] = m1[1][3]; m2[3][2] = m1[2][3]; m2[3][3] = m1[3][3];  
}

