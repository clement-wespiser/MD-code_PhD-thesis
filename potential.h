double Potential(double r, double theta, double phi);
double *SpheCoord(H2O *h2o, double *pos);
double **Univec(double *spheCoord);
double **GradSphe(H2O *h2o, double *pos);

double Potential(double r, double theta, double phi)
{
	double x = r*sin(theta)*cos(phi);
	double y = r*sin(theta)*sin(phi);
	double z = r*cos(theta);

	double rh = 0.9576;
	double ah = 0.912022;
	double xh = rh*cos(ah);
	double yh = rh*sin(ah);

	double **coords = alloc_2darray(3,3);
	zero_2darray(coords, 3, 3);
	coords[1][0] = xh;
	coords[1][1] = yh;
	coords[2][0] = xh;
	coords[2][1] = -yh;

//   --------------Oxygen Parameters       
    double re	=	3.8673;
    double ae	=	2.8135;
    double C0	=	-155.0498;
    double V0	=	570;
    double B3	=	0.6278;
    double B4	=	1.1452;
    double B5	=	0.5173;
    double B6	=	0.1105;
    double B7	=	-0.005;
    double B8	=	0;
    double B12	=	-1.0868;
	double B112	=	0.5207;
	double B122	=	-3.0264;
	double B1112	=	6.8625;
	double B1122	=	-0.9652;
	double B1222	=	-1.9038;
	double B11112	=	-12.928;
	double B11122	=	12.014;
	double B11222	=	-3.5689;
	double B12222	=	-0.5123;
	double B111112	=	10.045;
	double B111122	=	-15.281;
	double B111222	=	10.296;
	double B112222	=	-2.6911;
//   -------------Hydrogen Parameters   
    double B122222	=	0;
    double reH	=	2.965;
    double aeH	=	1.5495;
    double V0H	=	370;
    double C3	=	-0.2625;
    double C4	=	-4.2494;
    double C5	=	5.853;
    double C6	=	-3.4309;
    double C7	=	0;
    double C12	=	4.3595;
    double C112	=	-10.227;
    double C1122	=	15.846;
    double C1112	=	7.2077;
    double C11112	=	2.0236;
    double C11122	=	-8.4306;
    double C111112	=	0;
    double C111122	=	0;
    double C111222	=	0;
    double D123	=	16.793;
    double D1233	=	12.55;
    double D12333	=	3.4706;
    double D1123	=	-20.41;
    double D11233	=	-7.341;
    double D11223	=	22.07;
    double D112233	=	0;
    double D6	=	0.759;
    double De	=	142.8;
    double re0	=	3.66;

//======== Hydrogens-Oxygen========
   
    double u12=0;
    double u112=0;
    double u122=0;
    double u1112=0;
    double u1122=0;
    double u1222=0;
    double u11112=0;
    double u11122=0;
    double u11222=0;
    double u12222=0;
    double u111112=0;
    double u111122=0;
    double u111222=0;
    double u112222=0;
    double u122222=0;
//========== Hydrogens========== 
    double vh11=0;       
    double v12=0 ;
    double v112=0 ; 
    double v1122=0 ;
    double v1112=0;
    double v11112=0; 
    double v11122=0;
    double v111112=0;   
    double v111122=0;
    double v111222=0;
    double t123=0;
    double t1233=0;
    double t12333=0;
    double t1123=0;
    double t11233=0	;                                     
    double t11223=0;	
    double t112233=0;

//===Ar-Oxygen=======================
    
    double x1=x-coords[0][0];
    double y1=y-coords[0][1];
    double z1=z-coords[0][2];
    double xx=pow(x1,2) + pow(y1,2) + pow(z1,2);
    
    double r1 = sqrt(xx);
    double r11 = r1-re;
    
    double V1 = 1-exp(-r11*ae/re);
    double PO = pow(V1,2)*(1+V1*(B3+V1*(B4+V1*(B5+V1*(B6+V1*(B7+V1*B8))))));


//===Ar-Hydrogens====================
    for(int iii=0; iii<2; iii++)
	{
        double xk = x-coords[iii+1][0];
        double yk = y-coords[iii+1][1];
        double zk = z-coords[iii+1][2];
        double xx = pow(xk,2) + pow(yk,2) + pow(zk,2);
        double rk = sqrt(xx);
        double rkk = rk-reH;
        
        double Vk = 1-exp(-rkk*aeH/reH);
        
        double Wk = pow(Vk,2)*(1+Vk*(C3+Vk*(C4+Vk*(C5+Vk*(C6+Vk*C7)))));

        vh11=vh11+Wk;

        //Ar-Oxygen-Hydrogens
        double Vkk=Vk*Vk;
        double V11=V1*V1;
        double Vk1=Vk*V1;
        
        u12=u12+Vk1;
        u112=u112+Vk*Vk1;
        u122=u122+Vk1*V1;
        u1112=u1112+Vkk*Vk1;
        u1122=u1122+Vk1*Vk1;
        u1222=u1222+Vk1*V11;
        
        u11112=u11112+Vk*Vkk*Vk1;
        u11122=u11122+Vk*Vk1*Vk1;
        u11222=u11222+Vk1*Vk1*V1;
        u12222=u12222+Vk1*V11*V1;
        
        u111112=u111112+Vkk*Vkk*Vk1;
        u111122=u111122+Vkk*Vk1*Vk1;
        u111222=u111222+Vk1*Vk1*Vk1;
        u112222=u112222+Vk1*Vk1*V11;
        u122222=u122222+Vk1*V11*V11;


        //Ar-Hydrogen-Hydrogen
        if(iii == 0)
		{
            double xl = x-coords[2][0];
            double yl = y-coords[2][1];
            double zl = z-coords[2][2];
            double rl = sqrt(pow(xl,2) + pow(yl,2) + pow(zl,2));
            
            double rll=rl-reH;
            double Vl=1-exp(-rll*aeH/reH);
            
            double Vll=Vl*Vl;
            double Vkl=Vk*Vl;

            v12=v12+Vkl;
            v112=v112+(Vk+Vl)*Vkl;
            v1122=v1122+Vkl*Vkl;
            v1112=v1112+(Vkk+Vll)*Vkl;
            v11112=v11112+(Vkk*Vk+Vll*Vl)*Vkl;
            v11122=v11122+(Vk+Vl)*Vkl*Vkl;
            v111112=v111112+(Vkk*Vkk+Vll*Vll)*Vkl;
            v111122=v111122+(Vkk+Vll)*Vkl*Vkl;
            v111222=v111222+Vkl*Vkl*Vkl;

            //Ar-Oxygen-Hydrogen-Hydrogen
            t123=  t123+ Vkl*V1;
            t1233= t1233+Vkl*V11;
            t12333=t1233+Vkl*V11*V1 ;

            t1123=  t1123+ (Vk+Vl)*Vkl*V1;
            t11233= t11233+(Vk+Vl)*Vkl*V11;
            t11223= t11223+ Vkl*Vkl*V1;
            t112233=t112233+Vkl*Vkl*V11;
		}
	}
    
	double PHO=PO+u12*B12+u112*B112+u122*B122+u1122*B1122+u1112*B1112;
    PHO=PHO+u1222*B1222;
    PHO=PHO+u11112*B11112+u11122*B11122+u11222*B11222+u12222*B12222;
    PHO=PHO+u111112*B111112+u111122*B111122+u111222*B111222+u112222*B112222+u122222*B122222;
    double PH=vh11+v12*C12+v112*C112+v1122*C1122+v1112*C1112+v11112*C11112;
    PH=PH+v11122*C11122+v111112*C111112+v111122*C111122+v111222*C111222;
    PH=PH+t123*D123+t1233*D1233+t12333*D12333;
    PH=PH+t1123*D1123+t11233*D11233+t11223*D11223+t112233*D112233;

//-----Asymptotics    

    int gama=4;
    double h_short=1/(1+exp( gama*(r1-re0-re0)));
    double h_long= 1/(1+exp(-gama*(r1-re0-re0)));
    double DE6=D6*De;
    double V=(V0*PHO+V0H*PH+C0)*h_short-h_long*De*D6*pow(re0/r1,6);

	V *= 1.196266e-6; // conversion cm^-1 -> uma.A^2.fs^-2

	free2dArray(coords, 3);
	return V;
}

double *SpheCoord(H2O *h2o, double *pos)
/* pos : position de l'atome pour lequel calculer r, theta et phi
 */
{
	// === MFF ===
	double b[3] = {1,0,0};
	double a[3] = {0,1,0};
	double c[3] = {0,0,1};
	double dcm = ((h2o->m_h2 + h2o->m_h2) * h2o->d_oh *cos(h2o->a_oh/2))/h2o->m_h2o;
    double o[3]  = {-dcm, 0., 0.};  // centre de masse a l'origine du MFF
	// === SFF ===
	RotateVecQuat(b, h2o->q); // b
	RotateVecQuat(a, h2o->q); // a
	RotateVecQuat(c, h2o->q); // c
	RotateVecQuat(o, h2o->q); // position de l'atome d'O dans le SFF

	o[0] += h2o->r[0][0]; // translation
	o[1] += h2o->r[0][1];
	o[2] += h2o->r[0][2];

	double coordRel[3] = {0,0,0};
	// Position de l'atome par rapport a l'oxygène
	for(int i=0; i<3; i++)
		coordRel[i] = pos[i] - o[i];

	// coordonnées cartésienne du vecteur pos dans le MFF
	double y = DotProd(coordRel, a);
	double x = DotProd(coordRel, b);
	double z = DotProd(coordRel, c);
//	printf("x, y, z : %f %f %f\n", x, y, z);

	double r = sqrt(pow(y,2) + pow(x,2) + pow(z,2));
	double theta = acos(z/r);
	double phi = atan2(y,x);

	double *spheCoord = malloc(sizeof(double)*3);
	spheCoord[0] = r; spheCoord[1] = theta; spheCoord[2] = phi;

	return spheCoord;
}

double **Univec(double *spheCoord)
/* Retourne un array 2d avec chaque ligne : vecteurs unitaire r, theta et phi
 */
{
	double **uniVec = alloc_2darray(3,3);
	double r = spheCoord[0];
	double theta = spheCoord[1];
	double phi = spheCoord[2];

	uniVec[0][0] = cos(phi)*sin(theta);
	uniVec[0][1] = sin(phi)*sin(theta);
	uniVec[0][2] = cos(theta);

	uniVec[1][0] = cos(phi)*cos(theta);
	uniVec[1][1] = sin(phi)*cos(theta);
	uniVec[1][2] = -sin(theta);

	uniVec[2][0] = -sin(phi);
	uniVec[2][1] = cos(phi);
	uniVec[2][2] = 0;

	return uniVec;
}

double **GradSphe(H2O *h2o, double *pos)
{
	double dr = 0.00001; double dtheta = 0.00001; double dphi = 0.00001;
	double **gradSphe = alloc_2darray(4,3);
	double *spheCoord = SpheCoord(h2o, pos);
	double r = *spheCoord;	double theta = *(spheCoord+1);	double phi = *(spheCoord+2);

	double **univecs = Univec(spheCoord);
	double *r_uni = univecs[0];
	double *theta_uni = univecs[1];
	double *phi_uni = univecs[2];

	double norm_r = (Potential(r+dr, theta, phi) - Potential(r-dr, theta, phi))/(2*dr);
	double norm_theta = (Potential(r, theta+dtheta, phi) - Potential(r, theta-dtheta, phi))/(2*dtheta);
	double norm_phi = 0;
	if(theta == 0)
		norm_phi = 0;
	else
		norm_phi =(Potential(r, theta, phi+dphi) - Potential(r, theta, phi-dphi))/(2*dphi);

	gradSphe[0][0] = norm_r;
	gradSphe[0][1] = norm_theta;
	gradSphe[0][2] = norm_phi;
	gradSphe[1][0] = r_uni[0]; 
	gradSphe[1][1] = r_uni[1];
	gradSphe[1][2] = r_uni[2];
	gradSphe[2][0] = theta_uni[0]; 
	gradSphe[2][1] = theta_uni[1]; 
	gradSphe[2][2] = theta_uni[2]; 
	gradSphe[3][0] = phi_uni[0]; 
	gradSphe[3][1] = phi_uni[1]; 
	gradSphe[3][2] = phi_uni[2]; 

	free2dArray(univecs,3);
	free(spheCoord);

	return gradSphe;
}

void HarmonicPotential(H2O *h2o, System *sys)
{
	double k_x     = 0.05;  // constante de force pour translation (uma.A^-1.fs^-2)
	double k_theta = 0.01;  // constante de force pour rotation (uma.A^-2.fs^-2)
	double x_0     = 4;    // position d'équilibre (A)
	double theta_0 = 0;    // angle d'équilibre (rad)

	double x = h2o->r[0][0]; // position de l'orgine du MFF
	double x_mff[3] = {1,0,0};
	double x_sff[3] = {1,0,0};
	RotateVecQuat(x_mff, h2o->q);
	double theta = acos(DotProd(x_mff, x_sff));

	double v = 0.5*k_x*pow(x-x_0,2) + 0.5*k_theta*pow(theta-theta_0,2);
	sys->epot = v;
}

double *ForceHarmonicPotential(H2O *h2o)
{
	double k_x     = 0.05;  // constante de force pour translation (uma.A^-1.fs^-2)
	double k_theta = 0.01;  // constante de force pour rotation (uma.A^-2.fs^-2)
	double x_0     = 4;    // position d'équilibre (A)
	double theta_0 = 0;    // angle d'équilibre (rad)

	double x = h2o->r[0][0]; // position de l'orgine du MFF
	double x_mff[3] = {1,0,0};
	double x_sff[3] = {1,0,0};
	RotateVecQuat(x_mff, h2o->q);
	double theta = acos(DotProd(x_mff, x_sff));
	if(x_mff[1] < 0)
		theta *= -1;

	double fx     = -k_x*(x-x_0);
	double ftheta = -k_theta*(theta-theta_0); 

	double *f = malloc(sizeof(double)*2);
	f[0] = fx; f[1] = ftheta;

	return f;
}







