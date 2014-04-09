/*
 * AbstractState.h
 *
 *  Created on: 21 Dec 2013
 *      Author: jowr
 */

#ifndef ABSTRACTSTATE_H_
#define ABSTRACTSTATE_H_

#include "Cache/CachedElement.h"
#include "Exceptions.h"
#include "DataStructures.h"
#include "l10n/english.h"

#include "Fluids/CoolPropFluid.h"

namespace CoolProp {

//! The mother of all state classes
/*!
This class provides the basic properties based on interrelations of the
properties, their derivatives and the Helmholtz energy terms. It does not
provide the mechanism to update the values. This has to be implemented in
a subclass. Most functions are defined as virtual functions allowing us
redefine them later, for example to implement the TTSE technique. The
functions defined here are always used as a fall-back.

This base class does not perform any checks on the two-phase conditions and
alike. Most of the functions defined here only apply to compressible single
state substances. Make sure you are aware of all the assumptions we made
when using this class.

Add build table function to Abstract State
Interpolator inherit AS implemented by TTSE BICUBIC

*/
class AbstractState {
protected:

    /// Some administrative variables
    long _fluid_type;
    long _phase;
    bool _forceSinglePhase, _forceTwoPhase;

    bool isCompressibleFluid(void){
        return !(_fluid_type == FLUID_TYPE_INCOMPRESSIBLE_LIQUID
              || _fluid_type == FLUID_TYPE_INCOMPRESSIBLE_SOLUTION);
    }

    bool checkCompressible(void){
        if (!this->isCompressibleFluid()){throw ValueError(ERR_NOT_COMPRESSIBLE);}
        return true;
    }

    bool isHomogeneousPhase(void){
        return (this->_phase==iphase_liquid || this->_phase==iphase_gas || this->_phase == iphase_supercritical);
    }

    bool isTwoPhase(void){
        return (this->_phase==iphase_twophase);
    }

    bool checkTwoPhase(void){
        if (!this->isCompressibleFluid()){throw ValueError(ERR_NOT_A_TWO_PHASE_FLUID);}
        if (!this->isTwoPhase()&&!_forceTwoPhase){throw ValueError(ERR_NOT_A_TWO_PHASE_STATE);}
        return true;
    }

    bool checkSinglePhase(void){
        if (!this->isHomogeneousPhase()||!_forceSinglePhase){throw ValueError(ERR_NOT_A_TWO_PHASE_FUNCTION);}
        return true;
    }

    /// Two important points
    SimpleState _critical, _reducing;

    /// Molar mass [mol/kg]
    CachedElement _molar_mass;
    
    /// Universal gas constant [J/mol/K]
    CachedElement _gas_constant;

    /// Bulk values
    double _rhomolar, _T, _p, _Q, _R;

    /// Transport properties
    CachedElement _viscosity, _conductivity, _surface_tension;

    CachedElement _hmolar, _smolar, _logp, _logrhomolar, _cpmolar, _cvmolar, _speed_sound;

    CachedElement _fugacity_coefficient;

    /// Smoothing values
    double _rhospline, _dsplinedp, _dsplinedh;

    /// Cached low-level elements for in-place calculation of other properties
    CachedElement _alpha0, _dalpha0_dTau, _dalpha0_dDelta, _d2alpha0_dTau2, _d2alpha0_dDelta_dTau,
            _d2alpha0_dDelta2, _d3alpha0_dTau3, _d3alpha0_dDelta_dTau2, _d3alpha0_dDelta2_dTau,
            _d3alpha0_dDelta3, _alphar, _dalphar_dTau, _dalphar_dDelta, _d2alphar_dTau2, _d2alphar_dDelta_dTau,
            _d2alphar_dDelta2, _d3alphar_dTau3, _d3alphar_dDelta_dTau2, _d3alphar_dDelta2_dTau,
            _d3alphar_dDelta3;

    CachedElement _dalphar_dDelta_lim, _d2alphar_dDelta2_lim,
            _d2alphar_dDelta_dTau_lim, _d3alphar_dDelta2_dTau_lim;

    /// Two-Phase variables
    CachedElement _rhoLmolar, _rhoVmolar;

    // ----------------------------------------
    // Property accessors to be optionally implemented by the backend
    // for properties that are not always calculated
    // ----------------------------------------
    /// Using this backend, calculate the molar enthalpy in J/mol
    virtual long double calc_hmolar(void){throw NotImplementedError("calc_hmolar is not implemented for this backend");};
    /// Using this backend, calculate the molar entropy in J/mol/K
    virtual long double calc_smolar(void){throw NotImplementedError("calc_smolar is not implemented for this backend");};
    /// Using this backend, calculate the molar constant-pressure specific heat in J/mol/K
    virtual long double calc_cpmolar(void){throw NotImplementedError("calc_cpmolar is not implemented for this backend");};
    /// Using this backend, calculate the molar constant-volume specific heat in J/mol/K
    virtual long double calc_cvmolar(void){throw NotImplementedError("calc_cvmolar is not implemented for this backend");};
    /// Using this backend, calculate the speed of sound in m/s
    virtual long double calc_speed_sound(void){throw NotImplementedError("calc_speed_sound is not implemented for this backend");};
    /// Using this backend, calculate the isothermal compressibility \f$ \kappa = -\frac{1}{v}\left.\frac{\partial v}{\partial p}\right|_T=\frac{1}{\rho}\left.\frac{\partial \rho}{\partial p}\right|_T\f$  in 1/Pa
    virtual long double calc_isothermal_compressibility(void){throw NotImplementedError("calc_isothermal_compressibility is not implemented for this backend");};
    /// Using this backend, calculate the isobaric expansion coefficient \f$ \beta = \frac{1}{v}\left.\frac{\partial v}{\partial T}\right|_p = -\frac{1}{\rho}\left.\frac{\partial \rho}{\partial T}\right|_p\f$  in 1/K
    virtual long double calc_isobaric_expansion_coefficient(void){throw NotImplementedError("calc_isobaric_expansion_coefficient is not implemented for this backend");};
    /// Using this backend, calculate the viscosity in Pa-s
    virtual long double calc_viscosity(void){throw NotImplementedError("calc_viscosity is not implemented for this backend");};
    /// Using this backend, calculate the thermal conductivity in W/m/K
    virtual long double calc_conductivity(void){throw NotImplementedError("calc_conductivity is not implemented for this backend");};
    /// Using this backend, calculate the surface tension in N/m
    virtual long double calc_surface_tension(void){throw NotImplementedError("calc_surface_tension is not implemented for this backend");};
    /// Using this backend, calculate the molar mass in kg/mol
    virtual long double calc_molar_mass(void){throw NotImplementedError("calc_molar_mass is not implemented for this backend");};
    /// Using this backend, calculate the pressure in Pa
    virtual long double calc_pressure(void){throw NotImplementedError("calc_pressure is not implemented for this backend");};
    /// Using this backend, calculate the universal gas constant \f$R_u\f$ in J/mol/K
    virtual long double calc_gas_constant(void){throw NotImplementedError("calc_gas_constant is not implemented for this backend");};
    /// Using this backend, calculate the fugacity coefficient (dimensionless)
    virtual long double calc_fugacity_coefficient(int i){throw NotImplementedError("calc_fugacity_coefficient is not implemented for this backend");};


    // Derivatives of residual helmholtz energy
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r\f$ (dimensionless)
    virtual long double calc_alphar(void){throw NotImplementedError("calc_alphar is not implemented for this backend");};
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r_{\delta}\f$ (dimensionless)
    virtual long double calc_dalphar_dDelta(void){throw NotImplementedError("calc_dalphar_dDelta is not implemented for this backend");};
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r_{\tau}\f$ (dimensionless)
    virtual long double calc_dalphar_dTau(void){throw NotImplementedError("calc_dalphar_dTau is not implemented for this backend");};
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r_{\delta\delta}\f$ (dimensionless)
    virtual long double calc_d2alphar_dDelta2(void){throw NotImplementedError("calc_d2alphar_dDelta2 is not implemented for this backend");};
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r_{\delta\tau}\f$ (dimensionless)
    virtual long double calc_d2alphar_dDelta_dTau(void){throw NotImplementedError("calc_d2alphar_dDelta_dTau is not implemented for this backend");};
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r_{\tau\tau}\f$ (dimensionless)
    virtual long double calc_d2alphar_dTau2(void){throw NotImplementedError("calc_d2alphar_dTau2 is not implemented for this backend");};
    
    // Derivatives of ideal-gas helmholtz energy
    /// Using this backend, calculate the ideal-gas Helmholtz energy term \f$\alpha^0\f$ (dimensionless)
    virtual long double calc_alpha0(void){throw NotImplementedError("calc_alpha0 is not implemented for this backend");};
    /// Using this backend, calculate the ideal-gas Helmholtz energy term \f$\alpha^0_{\delta}\f$ (dimensionless)
    virtual long double calc_dalpha0_dDelta(void){throw NotImplementedError("calc_dalpha0_dDelta is not implemented for this backend");};
    /// Using this backend, calculate the ideal-gas Helmholtz energy term \f$\alpha^0_{\tau}\f$ (dimensionless)
    virtual long double calc_dalpha0_dTau(void){throw NotImplementedError("calc_dalpha0_dTau is not implemented for this backend");};
    /// Using this backend, calculate the ideal-gas Helmholtz energy term \f$\alpha^0_{\delta\delta}\f$ (dimensionless)
    virtual long double calc_d2alpha0_dDelta_dTau(void){throw NotImplementedError("calc_d2alpha0_dDelta_dTau is not implemented for this backend");};
    /// Using this backend, calculate the ideal-gas Helmholtz energy term \f$\alpha^0_{\delta\tau}\f$ (dimensionless)
    virtual long double calc_d2alpha0_dDelta2(void){throw NotImplementedError("calc_d2alpha0_dDelta2 is not implemented for this backend");};
    /// Using this backend, calculate the ideal-gas Helmholtz energy term \f$\alpha^0_{\tau\tau}\f$ (dimensionless)
    virtual long double calc_d2alpha0_dTau2(void){throw NotImplementedError("calc_d2alpha0_dTau2 is not implemented for this backend");};

    virtual void calc_reducing_state(void){throw NotImplementedError("calc_reducing_state is not implemented for this backend");};
    
public:
    AbstractState();
    virtual ~AbstractState();

    double _tau, _delta;
    
    /// A factory function to return a pointer to a new-allocated instance of one of the backends.
    /**
    Very Important!! : You must ensure to delete the backend instance that is created, otherwise there will be a memory leak

    The backend that is selected is based on the string passed in:
    
    1. If it starts with "REFPROP-", the REFPROP backend will be used.  The remaining part of the string should then 
       either be
       1. A pure or pseudo-pure fluid name (eg. "PROPANE" or "R410A"), yielding a REFPROPBackend instance.
       2. A string that encodes the components of the mixture with a vertical bar between them (e.g. "R32|R125"), yielding a REFPROPMixtureBackend instance.
    2. If it starts with "TTSE", the TTSE backend will be used, yielding a TTSEBackend instance
    3. If it starts with "BICUBIC", the BICUBIC backend will be used, yielding a BICUBICBackend instance

    */
    static AbstractState * factory(std::string fluid_string);
    
    bool clear();
    virtual void update(long input_pair, double Value1, double Value2) = 0;
    virtual void set_mole_fractions(const std::vector<double> &mole_fractions) = 0;
    virtual void set_mass_fractions(const std::vector<double> &mass_fractions) = 0;

    const SimpleState & get_reducing(){return _reducing;};

    // ----------------------------------------
    // Bulk properties - temperature and density are directly calculated every time
    // All other parameters are calculated on an as-needed basis
    // ----------------------------------------
    double T(void)  {return _T;};
    double rhomolar(void){return _rhomolar;};
    double p(void)  {return _p;};
    double Q(void)  {return _Q;};

    double molar_mass(void);
    double gas_constant(void);

    double hmolar(void);
    double smolar(void);
    double cpmolar(void);
    double cvmolar(void);
    double speed_sound(void);
    double isothermal_compressibility(void);
    double isobaric_expansion_coefficient(void);
    double fugacity_coefficient(int i);

    // ----------------------------------------
    // Transport properties
    // ----------------------------------------
    double viscosity(void);
    double conductivity(void);
    double surface_tension(void);

    //// ----------------------------------------
    //// Derivatives of properties
    //// ----------------------------------------
    //virtual double dvdp_constT(void);
    //virtual double dvdT_constp(void);

    //// Density
    //virtual double drhodh_constp(void);
    //virtual double drhodp_consth(void);
    //virtual double drhodp_constT(void);
    //virtual double drhodT_constp(void);
    //virtual double d2rhodh2_constp(void);
    //virtual double d2rhodhdp(void);
    //virtual double d2rhodhdQ(void);
    //virtual double d2rhodp2_constT(void);
    //virtual double d2rhodpdQ(void);
    //virtual double d2rhodT2_constp(void);
    //virtual double d2rhodTdp(void);

    //// Pressure
    //virtual double dpdrho_consth(void);
    //virtual double dpdrho_constT(void);
    //virtual double dpdT_consth(void);
    //virtual double dpdT_constrho(void);
    //virtual double d2pdrho2_constT(void);
    //virtual double d2pdrhodT(void);
    //virtual double d2pdT2_constrho(void);

    //// Enthalpy
    //virtual double dhdp_constrho(void);
    //virtual double dhdp_constT(void);
    //virtual double dhdrho_constp(void);
    //virtual double dhdrho_constT(void);
    //virtual double dhdT_constp(void);
    //virtual double dhdT_constrho(void);
    //virtual double d2hdp2_constT(void);
    //virtual double d2hdrho2_constT(void);
    //virtual double d2hdrhodT(void);
    //virtual double d2hdT2_constp(void);
    //virtual double d2hdT2_constrho(void);
    //virtual double d2hdTdp(void);

    //// Entropy
    //virtual double dsdp_constT(void);
    //virtual double dsdrho_constp(void);
    //virtual double dsdrho_constT(void);
    //virtual double dsdT_constp(void);
    //virtual double dsdT_constrho(void);
    //virtual double d2sdp2_constT(void);
    //virtual double d2sdrho2_constT(void);
    //virtual double d2sdrhodT(void);
    //virtual double d2sdT2_constp(void);
    //virtual double d2sdT2_constrho(void);
    //virtual double d2sdTdp(void);

    //// Fundamental derivative of gas dynamics
    //virtual double fundamental_derivative_of_gas_dynamics(void);
    //virtual double d2pdv2_consts(void);

    //// Other functions and derivatives
    //virtual double A(void);
    //virtual double B(void);
    //virtual double C(void);
    //virtual double Z(void);

    //virtual double dAdT_constrho(void);
    //virtual double dAdrho_constT(void);
    //// TODO: Add constXX qualifier
    //virtual double dBdT(void);
    //virtual double dCdT(void);
    //virtual double dZdDelta(void);
    //virtual double dZdTau(void);

    // ----------------------------------------
    // Helmholtz energy and derivatives
    // ----------------------------------------
    /// Return the derivative \f$ \alpha^0 \f$
    long double alpha0(void){
        if (!_alpha0) _alpha0 = calc_alpha0();
        return _alpha0;
    };
    long double dalpha0_dDelta(void){
        if (!_dalpha0_dDelta) _dalpha0_dDelta = calc_dalpha0_dDelta();
        return _dalpha0_dDelta;
    };
    long double dalpha0_dTau(void){
        if (!_dalpha0_dTau) _dalpha0_dTau = calc_dalpha0_dTau();
        return _dalpha0_dTau;
    };
    long double d2alpha0_dDelta2(void){
        if (!_d2alpha0_dDelta2) _d2alpha0_dDelta2 = calc_d2alpha0_dDelta2();
        return _d2alpha0_dDelta2;
    };
    long double d2alpha0_dDelta_dTau(void){
        if (!_d2alpha0_dDelta_dTau) _d2alpha0_dDelta_dTau = calc_d2alpha0_dDelta_dTau();
        return _d2alpha0_dDelta_dTau;
    };
    long double d2alpha0_dTau2(void){
        if (!_d2alpha0_dTau2) _d2alpha0_dTau2 = calc_d2alpha0_dTau2();
        return _d2alpha0_dTau2;
    };
    /*virtual double d3alpha0_dDelta3(void) = 0;
    virtual double d3alpha0_dDelta2_dTau(void) = 0;
    virtual double d3alpha0_dDelta_dTau2(void) = 0;
    virtual double d3alpha0_dTau3(void) = 0;
*/

    long double alphar(void){
        if (!_alphar) _alphar = calc_alphar();
        return _alphar;
    };
    long double dalphar_dDelta(void){
        if (!_dalphar_dDelta) _dalphar_dDelta = calc_dalphar_dDelta();
        return _dalphar_dDelta;
    };
    long double dalphar_dTau(void){
        if (!_dalphar_dTau) _dalphar_dTau = calc_dalphar_dTau();
        return _dalphar_dTau;
    };
    long double d2alphar_dDelta2(void){
        if (!_d2alphar_dDelta2) _d2alphar_dDelta2 = calc_d2alphar_dDelta2();
        return _d2alphar_dDelta2;
    };
    long double d2alphar_dDelta_dTau(void){
        if (!_d2alphar_dDelta_dTau) _d2alphar_dDelta_dTau = calc_d2alphar_dDelta_dTau();
        return _d2alphar_dDelta_dTau;
    };
    long double d2alphar_dTau2(void){
        if (!_d2alphar_dTau2) _d2alphar_dTau2 = calc_d2alphar_dTau2();
        return _d2alphar_dTau2;
    };
    /*
    virtual double d3alphar_dDelta3(void) = 0;
    virtual double d3alphar_dDelta2_dTau(void) = 0;
    virtual double d3alphar_dDelta_dTau2(void) = 0;
    virtual double d3alphar_dTau3(void) = 0;

    virtual double dalphar_dDelta_lim(void) = 0;
    virtual double d2alphar_dDelta2_lim(void) = 0;
    virtual double d2alphar_dDelta_dTau_lim(void) = 0;
    virtual double d3alphar_dDelta2_dTau_lim(void) = 0;
    */
};

} /* namespace CoolProp */
#endif /* ABSTRACTSTATE_H_ */
