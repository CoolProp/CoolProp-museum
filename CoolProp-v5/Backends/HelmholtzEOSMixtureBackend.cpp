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
    }
    else{
        is_pure_or_pseudopure = false;
    }

    // Set the mole-fraction weighted gas constant for the mixture 
    // (or the pure/pseudo-pure fluid) if it hasn't been set yet
    gas_constant();

    // Set the cache value for the molar mass if it hasn't been set yet
    molar_mass();

    // Reducing state
    calc_reducing_state();
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
    switch(input_pair)
    {
        case DmolarT_INPUTS:
        {
            _rhomolar = value1; _T = value2;

            _delta = _rhomolar/_reducing.rhomolar;
            _tau = _reducing.T/_T;

            // Calculate derivative if needed
            dalphar_dDelta();

            // Get pressure
            _p = _rhomolar*(double)_gas_constant*_T*(1+_delta*(double)_dalphar_dDelta);

            break;
        }
        case DmassT_INPUTS:
        {
            // Call again, but this time with molar units [kg/m^3] * [mol/kg] -> [mol/m^3]
            update(DmolarT_INPUTS, value1 / (double)_molar_mass, value2);
            return;
        }
    }
    
    throw std::exception();
}

void HelmholtzEOSMixtureBackend::calc_reducing_state(void)
{
    if (is_pure_or_pseudopure){
        _reducing = components[0]->EOSVector[0].reduce;
    }
    else{
        throw ValueError();
    }
}

double HelmholtzEOSMixtureBackend::calc_dalphar_dDelta(void)
{
    if (is_pure_or_pseudopure){
        return components[0]->dalphar_dDelta(_tau, _delta);
    }
    else{
        throw ValueError();
    }
}

} /* namespace CoolProp */
