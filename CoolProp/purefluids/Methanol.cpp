#include "CoolProp.h"
#include <vector>
#include "CPExceptions.h"
#include "FluidClass.h"
#include "Methanol.h"

MethanolClass::MethanolClass()
{
	double n[] = {0, 0.96352729792779e-1, -0.10848826325874e1, 0.29919647090261e-1, -0.17963419593895e-2, 0.47354317752015e-4, 0.10013578850486e+1, -0.12555691488591e+1, 0.85469725717500e0, -0.58295570793694e-1, 0.26935675584229e-1, 0.11504892676606e0, -0.51081766133636e-2, 0.19167368789348e-2, -0.28618221186953e0, 0.48168213019845e0, -0.33081091251828e0, 0.92842083313630e-1, -0.35936470747247e-1};
	double d[] = {0, 1, 1, 4, 5, 7, 1, 1, 3, 4, 5, 1, 7, 9, 2, 3, 4, 6, 7};
	double t[] = {0, -0.125, 1.500, 0.000, -0.875, 1.250, 0.250, 2.000, 1.750, 2.500, 2.375, 6.875, 5.875, 5.000, 18.500, 19.000, 17.500, 14.000, 12.000};
	double l[] = {0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 3};

    // Critical parameters
    crit.rho = 273;
    crit.p = PressureUnit(8103.5, UNIT_KPA);
    crit.T = 512.50;
    crit.v = 1.0/crit.rho;

    // Other fluid parameters
    params.molemass = 32.04216;
    params.Ttriple = 175.61;
	params.ptriple = 0.00018629;
    params.accentricfactor = 0.5625;
    params.R_u = 8.314472;

    // Limits of EOS
	limits.Tmin = params.Ttriple;
    limits.Tmax = 500.0;
    limits.pmax = 100000.0;
    limits.rhomax = 1000000.0*params.molemass;

	// Residual part
    phirlist.push_back(new phir_power(n,d,t,l,1,18,19));
	double m = 0.977118832, epsilonbar = 5.46341463, vbarn = 0.204481952, kappabar = 0.148852832e-2;
	phirlist.push_back(new phir_SAFT_associating_2B(m, epsilonbar, vbarn, kappabar));

	// Ideal-gas part
	phi0list.push_back(new phi0_lead(13.9864114647, 3.2006369296e3));
	phi0list.push_back(new phi0_logtau(3.1950423807804));
	phi0list.push_back(new phi0_power(-1.14289818828912e-3*crit.T, -1));
	phi0list.push_back(new phi0_power(-2.62687155181005e-7*pow(crit.T,2), -2));
	phi0list.push_back(new phi0_power(6.42610441977784e-11*pow(crit.T,3), -3));
	phi0list.push_back(new phi0_Planck_Einstein( 4.70118076896145, 3.7664265756));

    name.assign("Methanol");
    aliases.push_back("methanol");
    aliases.push_back(std::string("METHANOL"));
    REFPROPname.assign("METHANOL");

	BibTeXKeys.EOS = "Piazza-FPE-2013";
	BibTeXKeys.SURFACE_TENSION = "Mulero-JPCRD-2012";
}

double MethanolClass::psat(double T)
{
	// Max error is  0.0991069278587 % between 175.61 and 512.5999 K
    const double t[] = {0, 0.077, 0.3505, 0.3515, 0.3555, 0.3645, 0.371, 0.383, 4.833333333333333, 5.0, 10.166666666666666};
    const double N[] = {0, -9.7405903858600631, 3260402084.5156689, -4777161923.1118717, 2016877532.0827882, -811371244.62298751, 338648069.09319597, -27394516.566426165, -46.888287588943143, 46.324153446741605, -3.5670972200405648};
    double summer = 0, theta;
    theta = 1 - T/crit.T;
    for (int i=1; i<=10; i++)
    {
        summer += N[i]*pow(theta, t[i]);
    }
    return crit.p.Pa*exp(crit.T/T*summer);
}

double MethanolClass::rhosatL(double T)
{
	// Max error is  0.211291380987 % between 175.61 and 512.5999 K
    const double t[] = {0, 0.058, 0.08600000000000001, 0.101, 0.11800000000000001, 0.12000000000000001, 0.376, 0.39049999999999996, 3.1666666666666665, 19.5, 19.666666666666668};
    const double N[] = {0, -4643.5288082809275, 64964.325423474293, -173441.65066216467, 695926.85475472361, -582903.36698389531, 642.65365377040234, -542.98076795251507, 0.82561848492061929, 3861.8624982134825, -4141.7016509247669};
    double summer=0,theta;
    theta = 1 - T/crit.T;
	for (int i=1; i<=10; i++)
	{
		summer += N[i]*pow(theta,t[i]);
	}
	return crit.rho*(summer+1);
}

double MethanolClass::rhosatV(double T)
{
    // Maximum absolute error is 0.370753 % between 175.610001 K and 513.379990 K
    const double t[] = {0, 0.16666666666666666, 0.3333333333333333, 0.5, 0.6666666666666666, 0.8333333333333334, 1.0, 1.1666666666666667, 1.3333333333333333, 1.5, 2.0, 2.3333333333333335, 2.8333333333333335, 3.3333333333333335, 4.166666666666667, 4.833333333333333, 5.833333333333333, 6.333333333333333, 7.333333333333333, 8.5, 9.833333333333334, 11.333333333333334, 12.833333333333334, 15.666666666666666, 18.666666666666668, 21.666666666666668};
    const double N[] = {0, 7597.1931157028366, -302294.66965590609, 4388396.8619273659, -33603221.401948757, 154437053.95209551, -445246674.15078759, 795873412.58213127, -810870583.94155431, 357750028.59211421, 87423503.387399897, -330273149.09998244, 627733578.63101625, -833167128.85553205, 1457124162.2768552, -2420510495.7178898, 6120215016.4571819, -7389022540.083499, 4611943493.9085159, -3397415663.398644, 2439723597.0041828, -1604131016.1799223, 755743289.06624866, -215370091.3594358, 95044773.686583579, -29267867.091564715};
    double summer=0,theta;
    theta=1-T/reduce.T;	
	for (int i=1; i<=25; i++)
	{
		summer += N[i]*pow(theta,t[i]);
	}
	return reduce.rho*exp(reduce.T/T*summer);
}
double MethanolClass::surface_tension_T(double T)
{
	// Mulero, JPCRD 2012
	return 0.22421*pow(1-T/reduce.T,1.3355) - 0.21408*pow(1-T/reduce.T,1.677) + 0.083233*pow(1-T/reduce.T,4.4402);
}
