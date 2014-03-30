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

class HelmholtzEOSMixtureBackend : public AbstractState  {
protected:
    std::vector<CoolPropFluid*> components; ///< The components that are in use
    
    bool is_pure_or_pseudopure; ///< A flag for whether the substance is a pure or pseudo-pure fluid (true) or a mixture (false)
    std::vector<double> mole_fractions; ///< The mole fractions of the components
    std::vector<double> mole_fractions_liq, ///< The mole fractions of the saturated liquid 
                        mole_fractions_vap; ///< The mole fractions of the saturated vapor
    SimpleState _crit;
public:
    HelmholtzEOSMixtureBackend(){pReducing = NULL;};
    HelmholtzEOSMixtureBackend(std::vector<CoolPropFluid*> components);
    HelmholtzEOSMixtureBackend(std::vector<std::string> components);
    virtual ~HelmholtzEOSMixtureBackend(){delete(pReducing);};
    ReducingFunction *pReducing;
    ExcessTerm Excess;

    void update(long input_pair, double value1, double value2);

    /// Set the components of the mixture
    /**
    @param components The components that are to be used in this mixture
    */
    void set_components(std::vector<CoolPropFluid*> components);

    void set_reducing_function();
    void set_excess_term();

    /// Set the mole fractions
    /** 
    @param mole_fractions The vector of mole fractions of the components
    */
    void set_mole_fractions(const std::vector<double> &mole_fractions){this->mole_fractions = mole_fractions;};
    
    /// Set the mass fractions
    /** 
    @param mass_fractions The vector of mass fractions of the components
    */
    void set_mass_fractions(const std::vector<double> &mass_fractions){throw std::exception();};
    
    double calc_molar_mass(void);
    double calc_gas_constant(void);

    double calc_dalphar_dDelta(void);

    double calc_alphar_deriv_nocache(const int nTau, const int nDelta, const std::vector<double> & mole_fractions, double tau, double delta);
    
    void calc_reducing_state(void);
    void calc_reducing_state_nocache(const std::vector<double> & mole_fractions);

    void calc_pressure(void);

    double p_rhoT(long double rhomolar, long double T);

    // ***************************************************************
    // *************  PHASE DETERMINATION ROUTINES  ******************
    // ***************************************************************
    void DmolarT_phase_determination();
    void DmolarP_phase_determination();
    void PT_phase_determination();

    // ***************************************************************
    // *******************  FLASH ROUTINES  **************************
    // ***************************************************************
    void DmolarT_flash();
    void DmolarP_flash();
    void PT_flash();

    // ***************************************************************
    // *******************  SOLVER ROUTINES  *************************
    // ***************************************************************
    void solver_rho_Tp();
    long double solver_rho_Tp_SRK();
};

} /* namespace CoolProp */
#endif /* HELMHOLTZEOSMIXTUREBACKEND_H_ */
