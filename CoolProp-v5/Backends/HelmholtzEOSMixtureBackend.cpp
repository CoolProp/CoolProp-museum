/*
 * AbstractBackend.cpp
 *
 *  Created on: 20 Dec 2013
 *      Author: jowr
 */

#if defined(_MSC_VER)
#define _CRTDBG_MAP_ALLOC
#define _CRT_SECURE_NO_WARNINGS
#include <crtdbg.h>
#include <sys/stat.h>
#else
#include <sys/stat.h>
#endif

#include <string>
//#include "CoolProp.h"

#include "HelmholtzEOSMixtureBackend.h"

namespace CoolProp {

HelmholtzEOSMixtureBackend::HelmholtzEOSMixtureBackend(std::vector<CoolPropFluid*> components) {
    
    /// Set the components and associated flags
    set_components(components);
}

void HelmholtzEOSMixtureBackend::set_components(std::vector<CoolPropFluid*> components) {

    // Copy the components
    this->components = components;

    if (components.size() == 1){
        is_pure_or_pseudopure = true;
        mole_fractions = std::vector<double>(1,1);
        c = components[0];
    }
    else{
        is_pure_or_pseudopure = false;
    }
    
}
double HelmholtzEOSMixtureBackend::calc_gas_constant(void)
{
    double summer = 0;
    for (unsigned int i = 0; i < components.size(); ++i)
    {
        summer += mole_fractions[i]*components[i]->gas_constant();
    }
    return summer;
}

double HelmholtzEOSMixtureBackend::calc_molar_mass(void)
{
    double summer = 0;
    for (unsigned int i = 0; i < components.size(); ++i)
    {
        summer += mole_fractions[i]*components[i]->molar_mass();
    }
    return summer;
}

void HelmholtzEOSMixtureBackend::update(long input_pair, double value1, double value2 )
{
    clear();

    // Set the mole-fraction weighted gas constant for the mixture 
    // (or the pure/pseudo-pure fluid) if it hasn't been set yet
    gas_constant();

    // Set the cache value for the molar mass if it hasn't been set yet
    molar_mass();

    // Reducing state
    calc_reducing_state();
    
    switch(input_pair)
    {
        case DmolarT_INPUTS:
        {
            _rhomolar = value1; _T = value2;
            DmolarT_flash();
            break;
        }
        case DmassT_INPUTS:
        {
            // Call again, but this time with molar units [kg/m^3] * [mol/kg] -> [mol/m^3]
            update(DmolarT_INPUTS, value1 / (double)_molar_mass, value2);
            return;
        }
    }
    // Check the values that must always be set
    if (!ValidNumber(_p)){ throw ValueError("p is not a valid number");}
    if (!ValidNumber(_T)){ throw ValueError("T is not a valid number");}
    if (!ValidNumber(_rhomolar)){ throw ValueError("rhomolar is not a valid number");}
    if (!ValidNumber(_Q)){ throw ValueError("Q is not a valid number");}
}
void HelmholtzEOSMixtureBackend::DmolarT_phase_determination()
{
    if (_T < _crit.T)
	{
		// Start to think about the saturation stuff
		// First try to use the ancillary equations if you are far enough away
		// Ancillary equations are good to within 1% in pressure in general
		// Some industrial fluids might not be within 3%
        if (_rhomolar < 0.95*c->ancillaries.rhoV.evaluate(_T)){
            this->_phase = iphase_gas; return;
		}
        else if (_rhomolar > 1.05*c->ancillaries.rhoL.evaluate(_T)){
			this->_phase = iphase_liquid; return;
		}
		else{
			//// Actually have to use saturation information sadly
			//// For the given temperature, find the saturation state
			//// Run the saturation routines to determine the saturation densities and pressures
			//// Use the passed in variables to save calls to the saturation routine so the values can be re-used again
			//saturation_T(T, enabled_TTSE_LUT, pL, pV, rhoL, rhoV);
			//double Q = (1/rho-1/rhoL)/(1/rhoV-1/rhoL);
			//if (Q < -100*DBL_EPSILON){
			//	return iLiquid;
			//}
			//else if (Q > 1+100*DBL_EPSILON){
			//	return iGas;
			//}
			//else{
			//	return iTwoPhase;
			//}
		}
	}
	// Now check the states above the critical temperature.

    // Calculate the pressure if it is not already cached.
	calc_pressure();

    if (_T > _crit.T && _p > _crit.p){
        this->_phase = iphase_supercritical; return;
	}
	else if (_T > _crit.T && _p < _crit.p){
		this->_phase = iphase_gas; return;
	}
	else if (_T < _crit.T && _p > _crit.p){
		this->_phase = iphase_liquid; return;
	}
	/*else if (p < params.ptriple){
		return iphase_gas;
	}*/
	else{
		throw ValueError(format("phase cannot be determined"));
	}
}
void HelmholtzEOSMixtureBackend::DmolarT_flash()
{
    // Find the phase, while updating all internal variables possible
    DmolarT_phase_determination();

    if (!isHomogeneousPhase())
    {
        throw ValueError("twophase not implemented yet");
    }
    else
    {
        if (!ValidNumber(_p))
        {
            calc_pressure();
            _Q = -1;
        }
    }
}
void HelmholtzEOSMixtureBackend::calc_pressure(void)
{
    _delta = _rhomolar/_reducing.rhomolar;
    _tau = _reducing.T/_T;
    
    // Calculate derivative if needed
    dalphar_dDelta();

    // Get pressure
    _p = _rhomolar*(double)_gas_constant*_T*(1+_delta*(double)_dalphar_dDelta);
}
void HelmholtzEOSMixtureBackend::calc_reducing_state(void)
{
    if (is_pure_or_pseudopure){
        _reducing = components[0]->EOSVector[0].reduce;
    }
    else{
        throw ValueError();
    }
    _crit = _reducing;
}

double HelmholtzEOSMixtureBackend::calc_dalphar_dDelta(void)
{
    return c->EOSVector[0].dalphar_dDelta(_tau, _delta);
}

} /* namespace CoolProp */
