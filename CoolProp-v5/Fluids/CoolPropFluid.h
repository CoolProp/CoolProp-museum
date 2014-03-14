/*
 * CoolPropFluid.h
 *
 *  Created on: 20 Dec 2013
 *      Author: jowr
 */

#ifndef COOLPROPFLUID_H_
#define COOLPROPFLUID_H_


#include <string>
#include <vector>
#include <map>
#include <assert.h>
#include <iterator>
#include "../DataStructures.h"

namespace CoolProp {

struct BibTeXKeysStruct
{
	std::string EOS,
	            CP0,
	            VISCOSITY,
	            CONDUCTIVITY,
	            ECS_LENNARD_JONES,
	            ECS_FITS,
	            SURFACE_TENSION;
};

struct EnvironmentalFactorsStruct
{
	double GWP20, GWP100, GWP500, ODP, HH, PH, FH;
	std::string ASHRAE34;
};

/// The base class class for the Helmholtz energy terms
/**
The copying is a bit complicated, we use a pure-virtual clone() function that must be 
implemented by all subclasses so that the vector that contains the HE terms can be 
copied properly

\sa http://www.parashift.com/c++-faq-lite/copy-of-abc-via-clone.html
\sa http://www.parashift.com/c++-faq-lite/virtual-ctors.html
\sa http://www.learncpp.com/cpp-tutorial/912-shallow-vs-deep-copying/

Terms explicit in Helmholtz energy:

Term                               | Helmholtz Energy Contribution
----------                         | ------------------------------
alphar_power                       | \f$ \alpha_r=\left\lbrace\begin{array}{cc}\displaystyle\sum_i n_i \delta^{d_i} \tau^{t_i} & l_i=0\\ \displaystyle\sum_i n_i \delta^{d_i} \tau^{t_i} \exp(-\delta^{l_i}) & l_i\neq 0\end{array}\right.\f$
alphar_gaussian                    | \f$ \alpha_r=\displaystyle\sum_i n_i \delta^{d_i} \tau^{t_i} \exp(-\eta_i(\delta-\epsilon_i)^2-\beta_i(\tau-\gamma_i)^2)\f$
alphar_exponential                 | \f$ \alpha_r=\displaystyle\sum_i n_i \delta^{d_i} \tau^{t_i} \exp(-\gamma_i\delta^{l_i}) \f$
alphar_GERG2008_gaussian           | \f$ \alpha_r=\displaystyle\sum_i n_i \delta^{d_i} \tau^{t_i} \exp(-\eta_i(\delta-\epsilon_i)^2-\beta_i(\delta-\gamma_i))\f$
alphar_Lemmon2005                  | \f$ \alpha_r=\displaystyle\sum_i n_i \delta^{d_i} \tau^{t_i} \exp(-\delta^{l_i}) \exp(-\tau^{m_i})\f$
*/
class BaseHelmholtzTerm{
public:
	virtual BaseHelmholtzTerm* clone() const = 0;   // The Virtual (Copy) Constructor
};


class alphar_power : public BaseHelmholtzTerm
{
public:
	virtual alphar_power* clone() const { return new alphar_power(*this); };
};

/// A set of limits for the eos parameters
struct EOSLimits
{
	double Tmin, Tmax, rhomax, pmax;
};

class ViscosityCorrelation
{
public:
	double dilute(double T, double rhomolar);
	double residual(double T, double rhomolar);
	double critical(double T, double rhomolar);
};
class ThermalConductivityCorrelation
{
public:
	double dilute(double T, double rhomolar);
	double residual(double T, double rhomolar);
	double critical(double T, double rhomolar);
};
class SurfaceTensionCorrelation
{

};

/// A wrapper around the HE vector in order to be able to cleanly handle the 
/// copying and destruction without requiring a full deep copy of EquationOfState
struct HelmholtzVector{
	std::vector<BaseHelmholtzTerm*> v;
	/// Default constructor
	HelmholtzVector(){};
	/// Copy constructor (see http://stackoverflow.com/questions/6775394/copying-a-vector-of-pointers)
	HelmholtzVector(const HelmholtzVector &other)
	{
		std::transform(other.v.begin(), other.v.end(), std::back_inserter(v), std::mem_fun(&BaseHelmholtzTerm::clone));
	}

	void push_back(BaseHelmholtzTerm* alpha){ v.push_back(alpha); };
	unsigned int size(){ return v.size(); };
};

/// The core class for an equation of state
/** 
 This class holds the absolute minimum information to evaluate the equation 
 of state.  This includes the reducing state, limits on the equation of state,
 the coefficients for the Helmholtz derivative terms.

 It does NOT include derived parameters like specific heat, enthalpy, etc.
*/
class EquationOfState{
public:
	EquationOfState(){};
	SimpleState reduce; ///< Reducing state used for the EOS (usually, but not always, the critical point)
	EOSLimits limits; ///< Limits on the EOS
	double R_u; ///< The universal gas constant used for this EOS (usually, but not always, 8.314472 J/mol/K)
	HelmholtzVector alphar; ///< The residual Helmholtz energy
	HelmholtzVector alpha0; ///< The ideal-gas Helmholtz energy

	/// Validate the EOS that was just constructed
	void validate()
	{
		assert(R_u < 9 && R_u > 8);
		assert(alphar.size() > 0);
		assert(alpha0.size() > 0);
	};
};

/// A thermophysical property provider for critical and reducing values as well as derivatives of Helmholtz energy
/**
This fluid instance is populated using an entry from a JSON file
*/
class CoolPropFluid {
	protected:
		// Transport property data
		std::string ECSReferenceFluid; ///< A string that gives the name of the fluids that should be used for the ECS method for transport properties
		double ECS_qd; ///< The critical qd parameter for the Olchowy-Sengers cross-over term
    public:
		CoolPropFluid(){};
		virtual ~CoolPropFluid(){};

		std::string name; ///< The name of the fluid
		std::string REFPROPname; ///< The REFPROP-compliant name if REFPROP-"name" is not a compatible fluid name.  If not included, "name" is assumed to be a valid name for REFPROP
		std::string CAS; ///< The CAS number of the fluid
		std::vector <std::string> aliases; ///< A vector of aliases of names for the fluid

		std::vector<EquationOfState> EOSVector; ///< The equations of state that could be used for this fluid
		std::vector<ViscosityCorrelation*> ViscosityVector; ///< The viscosity correlations that could be used for this fluid
		std::vector<ThermalConductivityCorrelation*> ThermalConductivityVector; ///< The thermal conductivity correlations that could be used for this fluid
		std::vector<SurfaceTensionCorrelation*> SurfaceTensionVector; ///< The surface tension correlations that could be used for this fluid

		BibTeXKeysStruct BibTeXKeys;
		EnvironmentalFactorsStruct environment;
	
};

//#include "../Backends/HelmholtzEOSBackend.h"
//#include "../Backends/REFPROPBackend.h"
//
///// The base class for interpolators
//class AbstractInterpolator{};
//
///// BicubicInterpolator mixin
//class BicubicInterpolator{};
//
///// TTSE Interpolator mixin
//class TTSEInterpolator{};
////
//class HelmholtzTTSEInterpolator : public HelmholtzEOSBackend, TTSEInterpolator
//{
//};
//
//class HelmholtzBicubicInterpolator : public HelmholtzEOSBackend, BicubicInterpolator
//{
//};
//
//class REFPROPTTSEInterpolator : public REFPROPBackend, TTSEInterpolator
//{
//};
//
//class REFPROPBicubicInterpolator : public REFPROPBackend, BicubicInterpolator
//{
//};


} /* namespace CoolProp */
#endif /* COOLPROPFLUID_H_ */
