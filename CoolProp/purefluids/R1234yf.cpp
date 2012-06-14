/*
Properties for R1234yf
by Ian Bell

Thermo props from
M. Richter and M.O. McLinden and E.W. Lemmon, "Thermodynamic Properties of 2,3,3,3-Tetrafluoroprop-1-ene
(R1234yf): Vapor Pressure and p-rho-T Measurements and an Equation of State"
by Reiner Tillner-Roth and Hans Dieter Baehr, J. Chem. Eng. Data, v. 56, 2011, pp 3254-3264
*/

#if defined(_MSC_VER)
#define _CRTDBG_MAP_ALLOC
#define _CRT_SECURE_NO_WARNINGS
#include <stdlib.h>
#include <crtdbg.h>
#else
#include <stdlib.h>
#endif

#include <math.h>
#include "string.h"
#include "stdio.h"
#include "CoolProp.h"
#include "R1234yf.h"
#include "PropMacros.h"

static const double a[]={
	 0, //[0]
	 0.04592563, //[1]
	 1.546958, //[2]
	-2.355237, //[3]
	-0.4827835, //[4]
	 0.1758022, //[5]
	-1.210006, //[6]
	-0.6177084, //[7]
	 0.6805262, //[8]
	-0.6968555, //[9]
	-0.02695779, //[10]
	 1.389966, //[11]
	-0.4777136, //[12]
	-0.1975184, //[13]
	-1.147646, //[14]
	 0.0003428541 //[15]
};

static const double d[]={
	0, //[0]
	4, //[1]
	1, //[2]
	1, //[3]
	2, //[4]
	3, //[5]
	1, //[6]
	3, //[7]
	2, //[8]
	2, //[9]
	7, //[10]
	1, //[11]
	1, //[12]
	3, //[13]
	3, //[14]
	2, //[15]
};

static const double t[]={
	0, //[0]
	1, //[1]
	0.32, //[2]
	0.929, //[3]
	0.94, //[4]
	0.38, //[5]
	2.28, //[6]
	1.76, //[7]
	0.97, //[8]
	2.44, //[9]
	1.05, //[10]
	1.4, //[11]
	3, //[12]
	3.5, //[13]
	1, //[14]
	3.5 //[15]
};

static const double c[]={
0, //[0]
0, //[1]
0, //[2]
0, //[3]
0, //[4]
0, //[5]
2, //[6]
2, //[7]
1, //[8]
2, //[9]
1, //[10]
};

static const double alpha[]={
0, //[0]
0, //[1]
0, //[2]
0, //[3]
0, //[4]
0, //[5]
0, //[6]
0, //[7]
0, //[8]
0, //[9]
0, //[10]
1.02, //[11]
1.336, //[12]
1.055, //[13]
5.84, //[14]
16.2, //[15]
};

static const double beta[]={
0, //[0]
0, //[1]
0, //[2]
0, //[3]
0, //[4]
0, //[5]
0, //[6]
0, //[7]
0, //[8]
0, //[9]
0, //[10]
1.42, //[11]
2.31, //[12]
0.89, //[13]
80, //[14]
108, //[15]
};

static const double gamm[]={
0, //[0]
0, //[1]
0, //[2]
0, //[3]
0, //[4]
0, //[5]
0, //[6]
0, //[7]
0, //[8]
0, //[9]
0, //[10]
1.13, //[11]
0.67, //[12]
0.46, //[13]
1.28, //[14]
1.2 //[15]
};

static const double epsilon[]={
0, //[0]
0, //[1]
0, //[2]
0, //[3]
0, //[4]
0, //[5]
0, //[6]
0, //[7]
0, //[8]
0, //[9]
0, //[10]
0.712, //[11]
0.91, //[12]
0.677, //[13]
0.718, //[14]
1.64 //[15]
};

static const double v0[]={
	0.0,	//[0]
	7.549,	//[1]
	1.537,	//[2]
	2.030,	//[3]
	7.455,	//[4]
};
static const double u0[]={
	// each of the ui terms are divided by the critical temp {fie on you Mr. EOS-maker}
	0.0,		//[0]
	718/367.85,	//[1]
	877/367.85,	//[2]
	4465/367.85,//[3]
	1755/367.85,//[4]
};

R1234yfClass::R1234yfClass()
{
	std::vector<double> n_v(a,a+sizeof(a)/sizeof(double));
	std::vector<double> d_v(d,d+sizeof(d)/sizeof(double));
	std::vector<double> t_v(t,t+sizeof(t)/sizeof(double));
	std::vector<double> l_v(c,c+sizeof(c)/sizeof(double));
    std::vector<double> alpha_v(alpha,alpha+sizeof(alpha)/sizeof(double));
	std::vector<double> beta_v(beta,beta+sizeof(beta)/sizeof(double));
	std::vector<double> gamm_v(gamm,gamm+sizeof(gamm)/sizeof(double));
	std::vector<double> epsilon_v(epsilon,epsilon+sizeof(epsilon)/sizeof(double));
	std::vector<double> u0_v(u0,u0+sizeof(u0)/sizeof(double));
	std::vector<double> v0_v(v0,v0+sizeof(v0)/sizeof(double));

	phi_BC * phir_ = new phir_power(n_v,d_v,t_v,l_v,1,10);
    phi_BC * phirg_ = new phir_gaussian(n_v,d_v,t_v,alpha_v,epsilon_v,beta_v,gamm_v,11,15);
	
	phirlist.push_back(phir_);
    phirlist.push_back(phirg_);

	phi_BC * phi0_lead_ = new phi0_lead(-12.837928,8.042605);
	phi_BC * phi0_logtau_ = new phi0_logtau(4.944);
	phi_BC * phi0_Planck_Einstein_ = new phi0_Planck_Einstein(v0_v,u0_v,1,4);

	phi0list.push_back(phi0_lead_);
	phi0list.push_back(phi0_logtau_);
	phi0list.push_back(phi0_Planck_Einstein_);

	// Critical parameters
	crit.rho = 475.553441976;
	crit.p = 3382.2;
	crit.T = 367.85;
	crit.v = 1.0/crit.rho;

	// Other fluid parameters
	params.molemass = 114.0415928;
	params.Ttriple = 220;
	params.accentricfactor = 0.276;
	params.R_u = 8.314472;

	// Limits of EOS
	limits.Tmin = 220;
	limits.Tmax = 410.0;
	limits.pmax = 30000.0;
	limits.rhomax = 11.64*params.molemass;
	
	EOSReference.assign("Richter, M. and M.O. McLinden and E.W. Lemmon, \"Thermodynamic Properties of 2,3,3,3-Tetrafluoroprop-1-ene"
						"(R1234yf): Vapor Pressure and p-rho-T Measurements and an Equation of State\""
						", J. Chem. Eng. Data, v. 56, 2011, pp 3254-3264");
	TransportReference.assign("Surface Tension: Katsuyuki Tanaka, Yukihiro Higashi, \"Thermodynamic properties of HFO-1234yf (2,3,3,3-tetrafluoropropene)\", International Journal of Refrigeration 33 (2010) 474-479");

	name.assign("R1234yf");
}
double R1234yfClass::conductivity_Trho(double T, double rho)
{
	return _HUGE;
}
double R1234yfClass::viscosity_Trho(double T, double rho)
{
	return _HUGE;
}
double R1234yfClass::psat(double T)
{
	double THETA = 1-T/crit.T;
	double RHS = 0.0015481-3.29049*pow(THETA,3.68086)-6.51048*pow(THETA,0.972422);
	return exp(crit.T/T*RHS)*crit.p;
}
double R1234yfClass::rhosatL(double T)
{
	double THETA = 1-T/crit.T;
	return 5.03314154e+02+1.15871385e+03*pow(THETA,4.13963252e-01)+1.68771255e+02*pow(THETA,1.93705217e+00);
}
double R1234yfClass::rhosatV(double T)
{
	double THETA = 1-T/crit.T;
	double RHS = -0.036061-4.16439*pow(THETA,0.499635)-9.33292*pow(THETA,1.64976)-35.4903*pow(THETA,4.36857);
	return exp(RHS)*crit.rho;
}
double R1234yfClass::surface_tension_T(double T)
{
	return 0.05983*pow(1-T/reduce.T,1.367);
}