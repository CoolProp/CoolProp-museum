/*
 * AbstractBackend.h
 *
 *  Created on: 20 Dec 2013
 *      Author: jowr
 */

#ifndef HELMHOLTZEOSMIXTUREBACKEND_H_
#define HELMHOLTZEOSMIXTUREBACKEND_H_

#include "../AbstractState.h"
#include "../Fluids/CoolPropFluid.h"
#include <vector>
#include "ReducingFunctions.h"
#include "ExcessHEFunction.h"

namespace CoolProp {

class SaturationSolvers;

class HelmholtzEOSMixtureBackend : public AbstractState  {
protected:
    std::vector<CoolPropFluid*> components; ///< The components that are in use
    
    bool is_pure_or_pseudopure; ///< A flag for whether the substance is a pure or pseudo-pure fluid (true) or a mixture (false)
    std::vector<double> mole_fractions; ///< The mole fractions of the components
    std::vector<double> mole_fractions_liq, ///< The mole fractions of the saturated liquid 
                        mole_fractions_vap; ///< The mole fractions of the saturated vapor
    SimpleState _crit;
    int imposed_phase_index;
public:
    HelmholtzEOSMixtureBackend(){SatL = NULL; SatV = NULL; imposed_phase_index = -1;};
    HelmholtzEOSMixtureBackend(std::vector<CoolPropFluid*> components, bool generate_SatL_and_SatV = true);
    HelmholtzEOSMixtureBackend(std::vector<std::string> component_names, bool generate_SatL_and_SatV = true);
    virtual ~HelmholtzEOSMixtureBackend(){};
    ReducingFunctionContainer Reducing;
    ExcessTerm Excess;

    const std::vector<CoolPropFluid*> &get_components(){return components;};

    HelmholtzEOSMixtureBackend *SatL, *SatV; ///< 

    void update(long input_pair, double value1, double value2);

    /// Set the components of the mixture
    /**
    @param components The components that are to be used in this mixture
    */
    void set_components(std::vector<CoolPropFluid*> components, bool generate_SatL_and_SatV = true);

    /**
    \brief Specify the phase - this phase will always be used in calculations
    @param phase_index The index from CoolProp::phases
    */
    void specify_phase(int phase_index){imposed_phase_index = phase_index;};

    void set_reducing_function();
    void set_excess_term();

    /// Set the mole fractions
    /** 
    @param mole_fractions The vector of mole fractions of the components
    */
    void set_mole_fractions(const std::vector<double> &mole_fractions);

    const std::vector<double> &get_mole_fractions(){return mole_fractions;};
    
    /// Set the mass fractions
    /** 
    @param mass_fractions The vector of mass fractions of the components
    */
    void set_mass_fractions(const std::vector<double> &mass_fractions){throw std::exception();};
    
    long double calc_molar_mass(void);
    long double calc_gas_constant(void);

    long double calc_alphar(void);
    long double calc_dalphar_dDelta(void);
    long double calc_dalphar_dTau(void);
    long double calc_d2alphar_dDelta2(void);
    long double calc_d2alphar_dDelta_dTau(void);
    long double calc_d2alphar_dTau2(void);

    long double calc_alpha0(void);
    long double calc_dalpha0_dDelta(void);
    long double calc_dalpha0_dTau(void);
    long double calc_d2alpha0_dDelta2(void);
    long double calc_d2alpha0_dDelta_dTau(void);
    long double calc_d2alpha0_dTau2(void);

    long double calc_alphar_deriv_nocache(const int nTau, const int nDelta, const std::vector<double> & mole_fractions, const long double &tau, const long double &delta);
    
    /**
    \brief Take derivatives of the ideal-gas part of the Helmholtz energy, don't use any cached values, or store any cached values

    @param nTau How many derivatives with respect to \f$\tau\f$ to take
    @param nDelta How many derivatives with respect to \f$\delta\f$ to take
    @param mole_fractions Mole fractions
    @param tau Reciprocal reduced temperature where \f$\tau=T_r / T\f$
    @param delta Reduced density where \f$\delta = \rho / \rho_r \f$
    @param Tr Reducing temperature of the mixture [K]
    @param rhor Reducing molar density of the mixture [mol/m^3]

    \f[
    \alpha^0 = \displaystyle\sum_{i=1}^{N}x_i[\alpha^0_{oi}(\rho,T) + \ln x_i]
    \f]
    where in this case, we use the \f$\alpha^0\f$ for the given fluid, which uses the inputs \f$\tau_i\f$ and \f$\delta_i\f$, so we do the conversion between mixture and component reduced states with
    \f[
    \tau_i = \frac{T_{c,i}}{T} = \frac{\tau T_{c,i}}{T_r}
    \f]
    \f[
    \delta_i = \frac{\rho}{\rho_{c,i}} = \frac{\delta\rho_r}{\rho_{c,i}}
    \f]

    \sa Table B5, GERG 2008 from Kunz Wagner, JCED, 2012
    */
    
    long double calc_alpha0_deriv_nocache(const int nTau, const int nDelta, const std::vector<double> & mole_fractions, const long double &tau, const long double &delta, const long double &Tr, const long double &rhor);
    
    void calc_reducing_state(void);
    void calc_reducing_state_nocache(const std::vector<double> & mole_fractions);

    long double calc_pressure(void);
    long double calc_cvmolar(void);
    long double calc_cpmolar(void);
    long double calc_hmolar(void);
    long double calc_smolar(void);
    long double calc_speed_sound(void);

    double p_rhoT(long double rhomolar, long double T);

    // ***************************************************************
    // *************  PHASE DETERMINATION ROUTINES  ******************
    // ***************************************************************
    void DmolarT_phase_determination_pure_or_pseudopure();
    void DmolarP_phase_determination();
    void PT_phase_determination();

    // ***************************************************************
    // *******************  FLASH ROUTINES  **************************
    // ***************************************************************
    void DmolarT_flash();
    void DmolarP_flash();
    void PT_flash();
    void QT_flash();

    // ***************************************************************
    // *******************  SOLVER ROUTINES  *************************
    // ***************************************************************
        
    void solver_rho_Tp();
    long double solver_rho_Tp_SRK();
};


struct saturation_T_pure_Akasaka_options{
    bool use_guesses;
    long double omega, rhoL, rhoV, pL, pV;
};
class SaturationSolvers
{
public:
    static void saturation_T_pure_Akasaka(HelmholtzEOSMixtureBackend *HEOS, long double T, saturation_T_pure_Akasaka_options &options);
};


} /* namespace CoolProp */
#endif /* HELMHOLTZEOSMIXTUREBACKEND_H_ */
