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
            // Use molar units [kg/m^3] / [kg/mol] -> [mol/m^3]
            _rhomolar = value1 / (double)_molar_mass; _T = value2;
            DmolarT_flash();
            break;
        }
        case PT_INPUTS:
        {
            _p = value1; _T = value2;
            PT_flash();
            break;
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
            throw NotImplementedError("potentially two phase inputs not possible yet");
			//// Actually have to use saturation information sadly
			//// For the given temperature, find the saturation state
			//// Run the saturation routines to determine the saturation densities and pressures
			//// Use the passed in variables to save calls to the saturation routine so the values can be re-used again
			//saturation_T(T, enabled_TTSE_LUT, pL, pV, rhoL, rhoV);
			//_Q = (1/rho-1/rhoL)/(1/rhoV-1/rhoL);
			//if (Q < -100*DBL_EPSILON){
			//	this->_phase = iphase_liquid; return;
			//}
			//else if (Q > 1+100*DBL_EPSILON){
			//	this->_phase = iphase_gas; return;
			//}
			//else{
			//	this->_phase = iphase_twophase; return;
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

void HelmholtzEOSMixtureBackend::PT_phase_determination()
{
    if (_T < _crit.T)
	{
		// Start to think about the saturation stuff
		// First try to use the ancillary equations if you are far enough away
		// Ancillary equations are good to within 1% in pressure in general
		// Some industrial fluids might not be within 3%
        if (_p > 1.05*c->ancillaries.p.evaluate(_T)){
            this->_phase = iphase_liquid; return;
		}
        else if (_p < 0.95*c->ancillaries.p.evaluate(_T)){
			this->_phase = iphase_gas; return;
		}
		else{
            throw NotImplementedError("potentially two phase inputs not possible yet");
			//// Actually have to use saturation information sadly
			//// For the given temperature, find the saturation state
			//// Run the saturation routines to determine the saturation densities and pressures
			//// Use the passed in variables to save calls to the saturation routine so the values can be re-used again
			//saturation_T(T, enabled_TTSE_LUT, pL, pV, rhoL, rhoV);
			//double Q = (1/rho-1/rhoL)/(1/rhoV-1/rhoL);
			//if (Q < -100*DBL_EPSILON){
			//	this->_phase = iphase_liquid; return;
			//}
			//else if (Q > 1+100*DBL_EPSILON){
			//	this->_phase = iphase_gas; return;
			//}
			//else{
			//	this->_phase = iphase_twophase; return;
			//}
		}
	}
	// Now check the states above the critical temperature.
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
double HelmholtzEOSMixtureBackend::p_rhoT(long double rhomolar, long double T)
{
    _delta = rhomolar/_reducing.rhomolar;
    _tau = _reducing.T/T;
    
    // Calculate derivative if needed
    long double dalphar_dDelta = c->EOSVector[0].dalphar_dDelta(_tau, _delta);

    // Get pressure
    return rhomolar*(double)_gas_constant*T*(1+_delta*dalphar_dDelta);
}
void HelmholtzEOSMixtureBackend::solver_rho_Tp()
{
    long double molar_mass = (double)_molar_mass;
    long double rhomolar = solver_rho_Tp_SRK();
    double res = p_rhoT(rhomolar, _T);
    double rhomass = rhomolar*molar_mass;
}
long double HelmholtzEOSMixtureBackend::solver_rho_Tp_SRK()
{
    long double rhomolar, R = (double)_gas_constant, T = _T;
    long double accentric = c->EOSVector[0].accentric;

    // Use SRK to get preliminary guess for the density
    long double m = 0.480+1.574*accentric-0.176*pow(accentric,2);
    long double b = 0.08664*R*_reducing.T/_reducing.p;
    long double a = 0.42747*pow(R*_reducing.T,2)/_reducing.p*pow(1+m*(1-sqrt(_T/_reducing.T)),2);
    long double A = a*_p/pow(R*T,2);
    long double B = b*_p/(R*T);

    //Solve the cubic for solutions for Z = p/(rho*R*T)
    double Z0, Z1, Z2; int N;
    solve_cubic(1, -1, A-B-B*B, -A*B, N, Z0, Z1, Z2);

    // Determine the guess value
    if (N == 1){
        rhomolar = _p/(Z0*R*T);
    }
    else{
        long double rhomolar0 = _p/(Z0*R*T);
        long double rhomolar1 = _p/(Z1*R*T);
        long double rhomolar2 = _p/(Z2*R*T);
        switch(_phase)
        {
        case iphase_liquid:
            rhomolar = max3(rhomolar0, rhomolar1, rhomolar2); break;
        case iphase_gas:
            rhomolar = min3(rhomolar0, rhomolar1, rhomolar2); break;
        default:
            throw ValueError();
        };
    }
    return rhomolar;
}
void HelmholtzEOSMixtureBackend::PT_flash()
{
    // Find the phase, while updating all internal variables possible
    PT_phase_determination();

    if (!isHomogeneousPhase())
    {
        throw ValueError("twophase not implemented yet");
    }
    else
    {
        // Find density
        solver_rho_Tp();
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
