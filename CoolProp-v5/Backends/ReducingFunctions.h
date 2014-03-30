#ifndef DEPARTURE_FUNCTIONS_H
#define DEPARTURE_FUNCTIONS_H

#include <vector>
#include "../Fluids/CoolPropFluid.h"

namespace CoolProp{

typedef std::vector<std::vector<double> > STLMatrix;

/// A container for the mixing parameters for CoolProp mixtures
/**

*/
class MixingParameterLibrary
{
public:
    /// Map from sorted pair of CAS numbers to reducing parameter map.  The reducing parameter map is a map from key (string) to value (double)
    std::map<std::vector<std::string>, std::vector<Dictionary> > reducing_map;
    MixingParameterLibrary();

    /// Parse a term from GERG 2008
    void parse_Kunz_JCED_2012(Dictionary &d, rapidjson::Value &val, bool swapped)
    {
        d.add_number("gammaV", cpjson::get_double(val, "gammaV"));
        d.add_number("gammaT", cpjson::get_double(val, "gammaT"));
    
        double betaV = cpjson::get_double(val, "betaV");
        double betaT = cpjson::get_double(val, "betaT");
        if (swapped){
            d.add_number("betaV", 1/betaV);
            d.add_number("betaT", 1/betaT);
        }
        else{
            d.add_number("betaV", betaV);
            d.add_number("betaT", betaT);
        }
    };

    /// Parse a term from HFC mixtures
    void parse_Lemmon_JPCRD_2004(Dictionary &d, rapidjson::Value &val, bool swapped)
    {
        d.add_number("xi", cpjson::get_double(val, "xi"));
        d.add_number("zeta", cpjson::get_double(val, "zeta"));
    };

    /// Parse a term from Air
    void parse_Lemmon_JPCRD_2000(Dictionary &d, rapidjson::Value &val, bool swapped)
    {
        d.add_number("xi", cpjson::get_double(val, "xi"));
        d.add_number("zeta", cpjson::get_double(val, "zeta"));
    };
};

/*! 
An abstract base class for the reducing function to allow for
Lemmon-Jacobsen, GERG, or other reducing function to yield the
reducing parameters \f$ \bar\rho_r \f$ and \f$ T_r \f$
*/
class ReducingFunction
{
protected:
	unsigned int N;
public:
	ReducingFunction(){};
	virtual ~ReducingFunction(){};

    /// A factory function to generate the required reducing function
    static ReducingFunction *factory(const std::vector<CoolPropFluid*> &components);

	/// The reduced temperature
	virtual double Tr(const std::vector<double> &x) = 0;
	/// The derivative of reduced temperature with respect to component i mole fraction
	virtual double dTrdxi__constxj(const std::vector<double> &x, int i) = 0;
	/// The molar reducing density
	virtual double rhormolar(const std::vector<double> &x) = 0;
	///Derivative of the molar reducing density with respect to component i mole fraction
	virtual double drhormolardxi__constxj(const std::vector<double> &x, int i) = 0;

	virtual double d2rhormolardxi2__constxj(const std::vector<double> &x, int i) = 0;
	virtual double d2rhormolardxidxj(const std::vector<double> &x, int i, int j) = 0;
	virtual double d2Trdxi2__constxj(const std::vector<double> &x, int i) = 0;
	virtual double d2Trdxidxj(const std::vector<double> &x, int i, int j) = 0;

	/*! GERG 2004 Monograph equation 7.56:
	\f[
	\left(\frac{\partial}{\partial x_j}\left(n\left(\frac{\partial T_r}{\partial n_i} \right)_{n_j}\right)\right)_{x_i} = \left(\frac{\partial^2T_r}{\partial x_j \partial x_i}\right)-\left(\frac{\partial T_r}{\partial x_j}\right)_{x_i}-\sum_{k=1}^Nx_k\left(\frac{\partial^2T_r}{\partial x_j \partial x_k}\right)
	\f]
	*/
	double d_ndTrdni_dxj__constxi(const std::vector<double> &x, int i, int j);
	/*! GERG 2004 Monograph equation 7.55:
	\f[
	\left(\frac{\partial}{\partial x_j}\left(n\left(\frac{\partial \rho_r}{\partial n_i} \right)_{n_j}\right)\right)_{x_i} = \left(\frac{\partial^2\rho_r}{\partial x_j \partial x_i}\right)-\left(\frac{\partial \rho_r}{\partial x_j}\right)_{x_i}-\sum_{k=1}^Nx_k\left(\frac{\partial^2\rho_r}{\partial x_j \partial x_k}\right)
	\f]
	*/
	double d_ndrhorbardni_dxj__constxi(const std::vector<double> &x, int i, int j);

	double ndrhorbardni__constnj(const std::vector<double> &x, int i);
	double ndTrdni__constnj(const std::vector<double> &x, int i);
};

/*! 
The Reducing parameter model used by the GERG-2008 formulation to yield the
reducing parameters \f$ \rho_r \f$ and \f$ T_r \f$ and derivatives thereof
*/
class GERG2008ReducingFunction : public ReducingFunction
{
protected:
	STLMatrix v_c; ///< \f$ v_{c,ij} = \frac{1}{8}\left(v_{c,i}^{1/3}+v_{c,j}^{1/3}\right)^{3}\f$ from GERG-2008
	STLMatrix T_c; ///< \f$ T_{c,ij} = \sqrt{T_{c,i}T_{c,j}} \f$ from GERG=2008
	STLMatrix beta_v; ///< \f$ \beta_{v,ij} \f$ from GERG-2008
	STLMatrix gamma_v; ///< \f$ \gamma_{v,ij} \f$ from GERG-2008
	STLMatrix beta_T; ///< \f$ \beta_{T,ij} \f$ from GERG-2008
	STLMatrix gamma_T; ///< \f$ \gamma_{T,ij} \f$ from GERG-2008
	std::vector<CoolPropFluid *> pFluids; ///< List of pointers to fluids
public:
	GERG2008ReducingFunction(std::vector<CoolPropFluid *> pFluids, STLMatrix beta_v, STLMatrix gamma_v, STLMatrix beta_T, STLMatrix gamma_T)
	{
		this->pFluids = pFluids;
		this->beta_v = beta_v;
		this->gamma_v = gamma_v;
		this->beta_T = beta_T;
		this->gamma_T = gamma_T;
		this->N = pFluids.size();
		T_c.resize(N,std::vector<double>(N,0));
		v_c.resize(N,std::vector<double>(N,0));
		for (unsigned int i = 0; i < N; ++i)
		{
			for (unsigned int j = 0; j < N; j++)
			{
				T_c[i][j] = sqrt(pFluids[i]->pEOS->reduce.T*pFluids[j]->pEOS->reduce.T);
				v_c[i][j] = 1.0/8.0*pow(pow(pFluids[i]->pEOS->reduce.rhomolar, -1.0/3.0)+pow(pFluids[j]->pEOS->reduce.rhomolar, -1.0/3.0),(int)3);
			}
		}
	};

	/// Default destructor
	~GERG2008ReducingFunction(){};
	/// The reduced temperature
	double Tr(const std::vector<double> &x);
	/// The derivative of reduced temperature with respect to component i mole fraction
	double dTrdxi__constxj(const std::vector<double> &x, int i);
	/// The molar reducing density
	double rhormolar(const std::vector<double> &x);
	///Derivative of the molar reducing density with respect to component i mole fraction
	double drhormolardxi__constxj(const std::vector<double> &x, int i);
	double dvrmolardxi__constxj(const std::vector<double> &x, int i);

	double d2vrmolardxi2__constxj(const std::vector<double> &x, int i);
	double d2rhormolardxi2__constxj(const std::vector<double> &x, int i);
	double d2vrmolardxidxj(const std::vector<double> &x, int i, int j);
	double d2rhormolardxidxj(const std::vector<double> &x, int i, int j);
	double d2Trdxi2__constxj(const std::vector<double> &x, int i);
	double d2Trdxidxj(const std::vector<double> &x, int i, int j);

	double c_Y_ij(int i, int j, std::vector< std::vector< double> > &beta, std::vector< std::vector< double> > &gamma, std::vector< std::vector< double> > &Y_c);
	double c_Y_ji(int j, int i, std::vector< std::vector< double> > &beta, std::vector< std::vector< double> > &gamma, std::vector< std::vector< double> > &Y_c);
	double f_Y_ij(const std::vector<double> &x, int i, int j, std::vector< std::vector< double> > &beta);

	double dfYkidxi__constxk(const std::vector<double> &x, int k, int i,std::vector< std::vector< double> > &beta);
	double dfYikdxi__constxk(const std::vector<double> &x, int i, int k, std::vector< std::vector< double> > &beta);
	double d2fYkidxi2__constxk(const std::vector<double> &x, int k, int i, std::vector< std::vector< double> > &beta);
	double d2fYikdxi2__constxk(const std::vector<double> &x, int i, int k, std::vector< std::vector< double> > &beta);
	double d2fYijdxidxj(const std::vector<double> &x, int i, int k, std::vector< std::vector< double> > &beta);
};

/*! From Lemmon, JPCRD, 2000 for the properties of Dry Air, and also from Lemmon, JPCRD, 2004 for the properties of R404A, R410A, etc.	
\f[
\rho_r(\bar x) = \left[ \sum_{i=1}^m\frac{x_i}{\rho_{c_i}}+\sum_{i=1}^{m-1}\sum_{j=i+1}^{m}x_ix_j\zeta_{ij}\right]^{-1}
\f]
\f[
T_r(\bar x) = \sum_{i=1}^mx_iT_{c_i}+\sum_{i=1}^{m-1}\sum_{j=i+1}^mx_ix_j\xi_{ij}
\f]

These can be converted to the form of GERG by the following equations:
\f[
\beta_T = 1\ \ \ \ \beta_v = 1 
\f]
and
\f[
    \boxed{\gamma_T = \dfrac{T_{c0}+T_{c1}+\xi_{01}}{2\sqrt{T_{c0}T_{c1}}}}
\f]
and
\f[
    \boxed{\gamma_v = \dfrac{v_{c0}+v_{c1}+\zeta_{01}}{\frac{1}{4}\left(\frac{1}{\rho_{c,i}^{1/3}}+\frac{1}{\rho_{c,j}^{1/3}}\right)^{3}}}
\f]
*/
class LemmonAirHFCReducingFunction
{
public:
	/// Set the coefficients based on reducing parameters loaded from JSON
	static void convert_to_GERG(const std::vector<CoolPropFluid*> &pFluids, int i, int j, Dictionary d, double &beta_T, double &beta_v, double &gamma_T, double &gamma_v)
    {
        double xi_ij = d.get_number("xi");
        double zeta_ij = d.get_number("zeta");
	    beta_T = 1;
	    beta_v = 1;
	    gamma_T = (pFluids[i]->pEOS->reduce.T + pFluids[j]->pEOS->reduce.T + xi_ij)/(2*sqrt(pFluids[i]->pEOS->reduce.T*pFluids[j]->pEOS->reduce.T));
	    double v_i = 1/pFluids[i]->pEOS->reduce.rhomolar;
	    double v_j = 1/pFluids[j]->pEOS->reduce.rhomolar;
	    double one_third = 1.0/3.0;
	    gamma_v = (v_i + v_j + zeta_ij)/(0.25*pow(pow(v_i, one_third)+pow(v_j, one_third),(int)3));
    };
};


} /* namespace CoolProp */
#endif