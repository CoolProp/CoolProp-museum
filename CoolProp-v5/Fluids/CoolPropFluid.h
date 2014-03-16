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

/// A wrapper around the HE vector in order to be able to cleanly handle the 
/// copying and destruction without requiring a full deep copy of EquationOfState
struct HelmholtzVector{
    std::vector<BaseHelmholtzTerm*> v;
    /// Default constructor
    HelmholtzVector(){};

    void push_back(BaseHelmholtzTerm* alpha){ v.push_back(alpha); };
    unsigned int size(){ return v.size(); };
    std::vector<BaseHelmholtzTerm*>::iterator begin(){ return v.begin();};
    std::vector<BaseHelmholtzTerm*>::iterator end(){ return v.end();};
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
    double molar_mass;
    HelmholtzVector alphar_vector, ///< The residual Helmholtz energy
                    alpha0_vector; ///< The ideal-gas Helmholtz energy

    /// Validate the EOS that was just constructed
    void validate()
    {
        assert(R_u < 9 && R_u > 8);
        assert(molar_mass > 0.001 && molar_mass < 1);
        assert(alphar_vector.size() > 0);
        assert(alpha0_vector.size() > 0);
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

        double gas_constant(){ return EOSVector[0].R_u; };
        double molar_mass(){ return EOSVector[0].molar_mass; };

        /// The derivative of the residual Helmholtz energy
        double alphar(double tau, double delta)
        {
            double summer = 0;
            for (std::vector<BaseHelmholtzTerm*>::iterator it = EOSVector[0].alphar_vector.begin(); it != EOSVector[0].alphar_vector.end(); ++it)
                summer += (*it)->base(tau,delta);
            return summer;
        }
        // First derivative
        double dalphar_dDelta(double tau, double delta)
        {
            double summer = 0;
            for (std::vector<BaseHelmholtzTerm*>::iterator it = EOSVector[0].alphar_vector.begin(); it != EOSVector[0].alphar_vector.end(); ++it)
                summer += (*it)->dDelta(tau,delta);
            return summer;
        };
        double dalphar_dTau(double tau, double delta);
        // Second derivative
        double dalphar_dDelta2(double tau, double delta);
        double d2alphar_dDelta_dTau(double tau, double delta);
        double d2alphar_dTau2(double tau, double delta);
        // Third derivative
        double d3alphar_dDelta3(double tau, double delta);
        double d3alphar_dDelta2_dTau(double tau, double delta);
        double d3alphar_dDelta_dTau2(double tau, double delta);
        double d3alphar_dTau3(double tau, double delta);
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
