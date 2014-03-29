/*
 * CoolPropFluid.h
 *
 *  Created on: 20 Dec 2013
 *      Author: jowr
 */

#ifndef COOLPROPFLUID_H_
#define COOLPROPFLUID_H_

#include <numeric>
#include <string>
#include <vector>
#include <map>
#include <assert.h>
#include <iterator>

#include "../DataStructures.h"
#include "../Helmholtz.h"

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
/**
*/
class AncillaryFunction
{
private:
    std::vector<double> n, t, s;
    bool using_tau_r;
    double Tmax, Tmin, reducing_value, T_r;
    int type;
    enum ancillaryfunctiontypes{TYPE_NOT_EXPONENTIAL = 0, TYPE_EXPONENTIAL = 1};
    std::size_t N;
public:
    
    AncillaryFunction(){};
    AncillaryFunction(rapidjson::Value &json_code)
    {
        n = cpjson::get_double_array(json_code["n"]);
        t = cpjson::get_double_array(json_code["t"]);
        Tmin = cpjson::get_double(json_code,"Tmin");
        Tmax = cpjson::get_double(json_code,"Tmax");
        reducing_value = cpjson::get_double(json_code,"reducing_value");
        using_tau_r = cpjson::get_bool(json_code,"using_tau_r");
        T_r = cpjson::get_double(json_code,"T_r");
        std::string type = cpjson::get_string(json_code,"type");
        if (!type.compare("rhoLnoexp"))
            this->type = TYPE_NOT_EXPONENTIAL;
        else
            this->type = TYPE_EXPONENTIAL;
        this->N = n.size();
        s = n;
    };
    double evaluate(double T)
    {
        double THETA = 1-T/T_r;

        for (std::size_t i = 0; i < N; ++i)
        {
            s[i] = n[i]*pow(THETA, t[i]);
        }
        double summer = std::accumulate(s.begin(), s.end(), 0.0);
        
        if (type == TYPE_NOT_EXPONENTIAL)
        {    
            return reducing_value*(1+summer);
        }
        else
        {
            double tau_r_value;
            if (using_tau_r)
                tau_r_value = T_r/T;
            else
                tau_r_value = 1.0;
            return reducing_value*exp(tau_r_value*summer);
        }
    }
};

struct Ancillaries
{
    AncillaryFunction p,rhoL,rhoV;
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
    ~EquationOfState()
    {
	    while (!alpha0_vector.empty()){ 
            delete alpha0_vector.back();  alpha0_vector.pop_back(); 
        }
    };
    SimpleState reduce; ///< Reducing state used for the EOS (usually, but not always, the critical point)
    EOSLimits limits; ///< Limits on the EOS
    double R_u; ///< The universal gas constant used for this EOS (usually, but not always, 8.314472 J/mol/K)
    double molar_mass;
    double accentric; ///< The accentric factor \f$ \omega = -log_{10}\left(\frac{p_s(T/T_c=0.7)}{p_c}\right)-1\f$
    ResidualHelmholtzContainer alphar; ///< The residual Helmholtz energy
    std::vector<BaseHelmholtzTerm*> alpha0_vector; ///< The ideal-gas Helmholtz energy

    /// Validate the EOS that was just constructed
    void validate()
    {
        assert(R_u < 9 && R_u > 8);
        assert(molar_mass > 0.001 && molar_mass < 1);
    };
    long double baser(const double tau, const double delta) throw()
    {
        return alphar.base(tau, delta);
    };
    // First partials
    long double dalphar_dDelta(const double tau, const double delta) throw()
    {
        return alphar.dDelta(tau, delta);
    };
    long double dalphar_dTau(const double tau, const double delta) throw()
    {
        return alphar.dTau(tau, delta);
    };
    // Second partials
    long double d2alphar_dDelta2(const double tau, const double delta) throw()
    {
        return alphar.dDelta2(tau, delta);
    };
    long double d2alphar_dDelta_dTau(const double tau, const double delta) throw()
    {
        return alphar.dDelta_dTau(tau, delta);
    };
    long double d2alphar_dTau2(const double tau, const double delta) throw()
    {
        return alphar.dTau2(tau, delta);
    };
    // Third partials
    long double d3alphar_dDelta3(const double tau, const double delta) throw()
    {
        return alphar.dDelta3(tau, delta);
    };
    long double d3alphar_dDelta2_dTau(const double tau, const double delta) throw()
    {
        return alphar.dDelta2_dTau(tau, delta);
    };
    long double d3alphar_dDelta_dTau2(const double tau, const double delta) throw()
    {
        return alphar.dDelta_dTau2(tau, delta);
    };
    long double d3alphar_dTau3(const double tau, const double delta) throw()
    {
        return alphar.dTau3(tau, delta);
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
        EquationOfState *pEOS; ///< A pointer to the currently used EOS
        std::vector<EquationOfState> EOSVector; ///< The equations of state that could be used for this fluid

        std::string name; ///< The name of the fluid
        std::string REFPROPname; ///< The REFPROP-compliant name if REFPROP-"name" is not a compatible fluid name.  If not included, "name" is assumed to be a valid name for REFPROP
        std::string CAS; ///< The CAS number of the fluid
        std::vector <std::string> aliases; ///< A vector of aliases of names for the fluid

        std::vector<ViscosityCorrelation*> viscosity_vector; ///< The viscosity correlations that could be used for this fluid
        std::vector<ThermalConductivityCorrelation*> thermal_conductivity_vector; ///< The thermal conductivity correlations that could be used for this fluid
        std::vector<SurfaceTensionCorrelation*> surface_tension_vector; ///< The surface tension correlations that could be used for this fluid

        BibTeXKeysStruct BibTeXKeys;
        EnvironmentalFactorsStruct environment;
        Ancillaries ancillaries;

        double gas_constant(){ return pEOS->R_u; };
        double molar_mass(){ return pEOS->molar_mass; };
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
