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
#include "../Fluids/FluidLibrary.h"

namespace CoolProp {

HelmholtzEOSMixtureBackend::HelmholtzEOSMixtureBackend(std::vector<std::string> component_names) {
    std::vector<CoolPropFluid*> components;
    components.resize(component_names.size());

    for (unsigned int i = 0; i < components.size(); ++i)
    {
        components[i] = &(get_library().get(component_names[i]));
    }

    pReducing = NULL;

    /// Set the components and associated flags
    set_components(components);
}
HelmholtzEOSMixtureBackend::HelmholtzEOSMixtureBackend(std::vector<CoolPropFluid*> components) {
    pReducing = NULL;

    /// Set the components and associated flags
    set_components(components);
}
void HelmholtzEOSMixtureBackend::set_components(std::vector<CoolPropFluid*> components) {

    // Copy the components
    this->components = components;

    if (components.size() == 1){
        is_pure_or_pseudopure = true;
        mole_fractions = std::vector<double>(1, 1);
    }
    else{
        is_pure_or_pseudopure = false;
    }

    // Set the reducing model
    set_reducing_function();

    // Set the excess Helmholtz energy if a mixture
    if (!is_pure_or_pseudopure)
    {
        set_excess_term();
    }
}
void HelmholtzEOSMixtureBackend::set_mole_fractions(const std::vector<double> &mole_fractions)
{
    if (mole_fractions.size() != components.size())
    {
        throw ValueError(format("size of mole fraction vector [%d] does not equal that of component vector [%d]",mole_fractions.size(), components.size()));
    }
    this->mole_fractions = mole_fractions;
};
void HelmholtzEOSMixtureBackend::set_reducing_function()
{
    pReducing = ReducingFunction::factory(components);
}
void HelmholtzEOSMixtureBackend::set_excess_term()
{
    Excess.construct(components);
}
long double HelmholtzEOSMixtureBackend::calc_gas_constant(void)
{
    double summer = 0;
    for (unsigned int i = 0; i < components.size(); ++i)
    {
        summer += mole_fractions[i]*components[i]->gas_constant();
    }
    return summer;
}
long double HelmholtzEOSMixtureBackend::calc_molar_mass(void)
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

    if (is_pure_or_pseudopure == false && mole_fractions.size() == 0) { 
        throw ValueError("Mole fractions must be set"); 
    }

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
        case QT_INPUTS:
        {
            _Q = value1; _T = value2;
            QT_flash();
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
        if (_rhomolar < 0.95*components[0]->ancillaries.rhoV.evaluate(_T)){
            this->_phase = iphase_gas; return;
		}
        else if (_rhomolar > 1.05*components[0]->ancillaries.rhoL.evaluate(_T)){
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
        if (_p > 1.05*components[0]->ancillaries.p.evaluate(_T)){
            this->_phase = iphase_liquid; return;
		}
        else if (_p < 0.95*components[0]->ancillaries.p.evaluate(_T)){
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
void HelmholtzEOSMixtureBackend::QT_flash()
{
    if (is_pure_or_pseudopure)
    {

    }
    else
    {

    }
}
void HelmholtzEOSMixtureBackend::DmolarT_flash()
{
    if (is_pure_or_pseudopure)
    {
        // Find the phase, while updating all internal variables possible
        DmolarT_phase_determination();
    }
    else
    {
        // TODO: hard coded phase
        _phase  = iphase_gas;
    }

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
// TODO: no cache
double HelmholtzEOSMixtureBackend::p_rhoT(long double rhomolar, long double T)
{
    _delta = rhomolar/_reducing.rhomolar;
    _tau = _reducing.T/T;
    
    // Calculate derivative if needed
    long double dalphar_dDelta = components[0]->pEOS->dalphar_dDelta(_tau, _delta);

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
    long double accentric = components[0]->pEOS->accentric;

    // Use SRK to get preliminary guess for the density
    long double m = 0.480+1.574*accentric-0.176*pow(accentric,2);
    long double b = 0.08664*R*_reducing.T/_reducing.p;
    long double a = 0.42747*pow(R*_reducing.T,2)/_reducing.p*pow(1+m*(1-sqrt(_T/_reducing.T)),2);
    long double A = a*_p/pow(R*T,2);
    long double B = b*_p/(R*T);

    //Solve the cubic for solutions for Z = p/(rho*R*T)
    double Z0, Z1, Z2; int Nsolns;
    solve_cubic(1, -1, A-B-B*B, -A*B, Nsolns, Z0, Z1, Z2);

    // Determine the guess value
    if (Nsolns == 1){
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
long double HelmholtzEOSMixtureBackend::calc_pressure(void)
{    
    // Calculate the reducing parameters
    _delta = _rhomolar/_reducing.rhomolar;
    _tau = _reducing.T/_T;

    // Calculate derivative if needed
    long double dar_dDelta = dalphar_dDelta();
    long double R_u = static_cast<double>(_gas_constant);

    // Get pressure
    _p = _rhomolar*R_u*_T*(1+_delta*dar_dDelta);

    return static_cast<double>(_p);
}
long double HelmholtzEOSMixtureBackend::calc_cvmolar(void)
{
    // Calculate the reducing parameters
    _delta = _rhomolar/_reducing.rhomolar;
    _tau = _reducing.T/_T;

    // Calculate derivatives if needed, or just use cached values
    long double d2ar_dTau2 = d2alphar_dTau2();
    long double d2a0_dTau2 = d2alpha0_dTau2();
    long double R_u = static_cast<double>(_gas_constant);

    // Get cv
    _cvmolar = -R_u*pow(_tau,2)*(d2ar_dTau2 + d2a0_dTau2);

    return static_cast<double>(_cvmolar);
}
long double HelmholtzEOSMixtureBackend::calc_cpmolar(void)
{
    // Calculate the reducing parameters
    _delta = _rhomolar/_reducing.rhomolar;
    _tau = _reducing.T/_T;

    // Calculate derivatives if needed, or just use cached values
    long double d2a0_dTau2 = d2alpha0_dTau2();
    long double dar_dDelta = dalphar_dDelta();
    long double d2ar_dDelta2 = d2alphar_dDelta2();
    long double d2ar_dDelta_dTau = d2alphar_dDelta_dTau();
    long double d2ar_dTau2 = d2alphar_dTau2();
    long double R_u = static_cast<double>(_gas_constant);

    // Get cp
    _cpmolar = R_u*(-pow(_tau,2)*(d2ar_dTau2 + d2a0_dTau2)+pow(1+_delta*dar_dDelta-_delta*_tau*d2ar_dDelta_dTau,2)/(1+2*_delta*dar_dDelta+pow(_delta,2)*d2ar_dDelta2));

    return static_cast<double>(_cpmolar);
}

void HelmholtzEOSMixtureBackend::calc_reducing_state_nocache(const std::vector<double> & mole_fractions)
{
    if (is_pure_or_pseudopure){
        _reducing = components[0]->pEOS->reduce;
        _crit = _reducing;
    }
    else{
        _reducing.T = pReducing->Tr(mole_fractions);
        _reducing.rhomolar = pReducing->rhormolar(mole_fractions);
    }
}
void HelmholtzEOSMixtureBackend::calc_reducing_state(void)
{
    calc_reducing_state_nocache(mole_fractions);
}
long double HelmholtzEOSMixtureBackend::calc_alphar_deriv_nocache(const int nTau, const int nDelta, const std::vector<double> &mole_fractions, const long double &tau, const long double &delta)
{
    if (is_pure_or_pseudopure)
    {
        if (nTau == 0 && nDelta == 0){
            return components[0]->pEOS->baser(tau, delta);
        }
        else if (nTau == 0 && nDelta == 1){
            return components[0]->pEOS->dalphar_dDelta(tau, delta);
        }
        else if (nTau == 1 && nDelta == 0){
            return components[0]->pEOS->dalphar_dTau(tau, delta);
        }
        else if (nTau == 0 && nDelta == 2){
            return components[0]->pEOS->d2alphar_dDelta2(tau, delta);
        }
        else if (nTau == 1 && nDelta == 1){
            return components[0]->pEOS->d2alphar_dDelta_dTau(tau, delta);
        }
        else if (nTau == 2 && nDelta == 0){
            return components[0]->pEOS->d2alphar_dTau2(tau, delta);
        }
        else if (nTau == 0 && nDelta == 3){
            return components[0]->pEOS->d3alphar_dDelta3(tau, delta);
        }
        else if (nTau == 1 && nDelta == 2){
            return components[0]->pEOS->d3alphar_dDelta2_dTau(tau, delta);
        }
        else if (nTau == 2 && nDelta == 1){
            return components[0]->pEOS->d3alphar_dDelta_dTau2(tau, delta);
        }
        else if (nTau == 3 && nDelta == 0){
            return components[0]->pEOS->d3alphar_dTau3(tau, delta);
        }
        else 
        {
            throw ValueError();
        }
    }
    else{
        unsigned int N = mole_fractions.size();
        long double summer = 0;
        if (nTau == 0 && nDelta == 0){
            for (unsigned int i = 0; i < N; ++i){ summer += mole_fractions[i]*components[i]->pEOS->baser(tau, delta); }
            return summer + Excess.alphar(tau, delta, mole_fractions);
        }
        else if (nTau == 0 && nDelta == 1){
            for (unsigned int i = 0; i < N; ++i){ summer += mole_fractions[i]*components[i]->pEOS->dalphar_dDelta(tau, delta); }
            return summer + Excess.dalphar_dDelta(tau, delta, mole_fractions);
        }
        else if (nTau == 1 && nDelta == 0){
            for (unsigned int i = 0; i < N; ++i){ summer += mole_fractions[i]*components[i]->pEOS->dalphar_dTau(tau, delta); }
            return summer + Excess.dalphar_dTau(tau, delta, mole_fractions);
        }
        else if (nTau == 0 && nDelta == 2){
            for (unsigned int i = 0; i < N; ++i){ summer += mole_fractions[i]*components[i]->pEOS->d2alphar_dDelta2(tau, delta); }
            return summer + Excess.d2alphar_dDelta2(tau, delta, mole_fractions);
        }
        else if (nTau == 1 && nDelta == 1){
            for (unsigned int i = 0; i < N; ++i){ summer += mole_fractions[i]*components[i]->pEOS->d2alphar_dDelta_dTau(tau, delta); }
            return summer + Excess.d2alphar_dDelta_dTau(tau, delta, mole_fractions);
        }
        else if (nTau == 2 && nDelta == 0){
            for (unsigned int i = 0; i < N; ++i){ summer += mole_fractions[i]*components[i]->pEOS->d2alphar_dTau2(tau, delta); }
            return summer + Excess.d2alphar_dTau2(tau, delta, mole_fractions);
        }
        /*else if (nTau == 0 && nDelta == 3){
            for (unsigned int i = 0; i < N; ++i){ summer += mole_fractions[i]*components[i]->pEOS->d3alphar_dDelta3(tau, delta); }
            return summer + pExcess.d3alphar_dDelta3(tau, delta);
        }
        else if (nTau == 1 && nDelta == 2){
            for (unsigned int i = 0; i < N; ++i){ summer += mole_fractions[i]*components[i]->pEOS->d3alphar_dDelta2_dTau(tau, delta); }
            return summer + pExcess.d3alphar_dDelta2_dTau(tau, delta);
        }
        else if (nTau == 2 && nDelta == 1){
            for (unsigned int i = 0; i < N; ++i){ summer += mole_fractions[i]*components[i]->pEOS->d3alphar_dDelta_dTau2(tau, delta); }
            return summer + pExcess.d3alphar_dDelta_dTau2(tau, delta);
        }
        else if (nTau == 3 && nDelta == 0){
            for (unsigned int i = 0; i < N; ++i){ summer += mole_fractions[i]*components[i]->pEOS->d3alphar_dTau3(tau, delta); }
            return summer + pExcess.d3alphar_dTau3(tau, delta);
        }*/
        else 
        {
            throw ValueError();
        }
    }
}
long double HelmholtzEOSMixtureBackend::calc_alpha0_deriv_nocache(const int nTau, const int nDelta, const std::vector<double> &mole_fractions, 
                                                                  const long double &tau, const long double &delta, const long double &Tr, const long double &rhor)
{
    if (is_pure_or_pseudopure)
    {
        if (nTau == 0 && nDelta == 0){
            return components[0]->pEOS->base0(tau, delta);
        }
        else if (nTau == 0 && nDelta == 1){
            return components[0]->pEOS->dalpha0_dDelta(tau, delta);
        }
        else if (nTau == 1 && nDelta == 0){
            return components[0]->pEOS->dalpha0_dTau(tau, delta);
        }
        else if (nTau == 0 && nDelta == 2){
            return components[0]->pEOS->d2alpha0_dDelta2(tau, delta);
        }
        else if (nTau == 1 && nDelta == 1){
            return components[0]->pEOS->d2alpha0_dDelta_dTau(tau, delta);
        }
        else if (nTau == 2 && nDelta == 0){
            return components[0]->pEOS->d2alpha0_dTau2(tau, delta);
        }
        else if (nTau == 0 && nDelta == 3){
            return components[0]->pEOS->d3alpha0_dDelta3(tau, delta);
        }
        else if (nTau == 1 && nDelta == 2){
            return components[0]->pEOS->d3alpha0_dDelta2_dTau(tau, delta);
        }
        else if (nTau == 2 && nDelta == 1){
            return components[0]->pEOS->d3alpha0_dDelta_dTau2(tau, delta);
        }
        else if (nTau == 3 && nDelta == 0){
            return components[0]->pEOS->d3alpha0_dTau3(tau, delta);
        }
        else 
        {
            throw ValueError();
        }
    }
    else{
        // See Table B5, GERG 2008 from Kunz Wagner, JCED, 2012
        unsigned int N = mole_fractions.size();
        long double summer = 0;
        long double tau_i, delta_i, rho_ci, T_ci;
        for (unsigned int i = 0; i < N; ++i){ 
            rho_ci = components[i]->pEOS->reduce.rhomolar; 
            T_ci = components[i]->pEOS->reduce.T;
            tau_i = T_ci*tau/Tr;
            delta_i = delta*rhor/rho_ci;

            if (nTau == 0 && nDelta == 0){    
                summer += mole_fractions[i]*(components[i]->pEOS->base0(tau_i, delta_i)+log(mole_fractions[i])); 
            }
            else if (nTau == 0 && nDelta == 1){
                summer += mole_fractions[i]*rhor/rho_ci*components[i]->pEOS->dalpha0_dDelta(tau_i, delta_i); 
            }
            else if (nTau == 1 && nDelta == 0){
                summer += mole_fractions[i]*T_ci/Tr*components[i]->pEOS->dalpha0_dTau(tau_i, delta_i); 
            }
            else if (nTau == 0 && nDelta == 2){
                summer += mole_fractions[i]*pow(rhor/rho_ci,2)*components[i]->pEOS->d2alpha0_dDelta2(tau_i, delta_i); 
            }
            else if (nTau == 1 && nDelta == 1){
                summer += mole_fractions[i]*rhor/rho_ci*T_ci/Tr*components[i]->pEOS->d2alpha0_dDelta_dTau(tau_i, delta_i); 
            }
            else if (nTau == 2 && nDelta == 0){
                summer += mole_fractions[i]*pow(T_ci/Tr,2)*components[i]->pEOS->d2alpha0_dTau2(tau_i, delta_i); 
            }
            else 
            {
                throw ValueError();
            }
        }
        return summer;
    }
}
long double HelmholtzEOSMixtureBackend::calc_alphar(void)
{
    const int nTau = 0, nDelta = 0;
    return calc_alphar_deriv_nocache(nTau, nDelta, mole_fractions, _tau, _delta);
}
long double HelmholtzEOSMixtureBackend::calc_dalphar_dDelta(void)
{
    const int nTau = 0, nDelta = 1;
    return calc_alphar_deriv_nocache(nTau, nDelta, mole_fractions, _tau, _delta);
}
long double HelmholtzEOSMixtureBackend::calc_dalphar_dTau(void)
{
    const int nTau = 1, nDelta = 0;
    return calc_alphar_deriv_nocache(nTau, nDelta, mole_fractions, _tau, _delta);
}
long double HelmholtzEOSMixtureBackend::calc_d2alphar_dTau2(void)
{
    const int nTau = 2, nDelta = 0;
    return calc_alphar_deriv_nocache(nTau, nDelta, mole_fractions, _tau, _delta);
}
long double HelmholtzEOSMixtureBackend::calc_d2alphar_dDelta_dTau(void)
{
    const int nTau = 1, nDelta = 1;
    return calc_alphar_deriv_nocache(nTau, nDelta, mole_fractions, _tau, _delta);
}
long double HelmholtzEOSMixtureBackend::calc_d2alphar_dDelta2(void)
{
    const int nTau = 0, nDelta = 2;
    return calc_alphar_deriv_nocache(nTau, nDelta, mole_fractions, _tau, _delta);
}

long double HelmholtzEOSMixtureBackend::calc_alpha0(void)
{
    const int nTau = 0, nDelta = 0;
    return calc_alpha0_deriv_nocache(nTau, nDelta, mole_fractions, _tau, _delta, _reducing.T, _reducing.rhomolar);
}
long double HelmholtzEOSMixtureBackend::calc_dalpha0_dDelta(void)
{
    const int nTau = 0, nDelta = 1;
    return calc_alpha0_deriv_nocache(nTau, nDelta, mole_fractions, _tau, _delta, _reducing.T, _reducing.rhomolar);
}
long double HelmholtzEOSMixtureBackend::calc_dalpha0_dTau(void)
{
    const int nTau = 1, nDelta = 0;
    return calc_alpha0_deriv_nocache(nTau, nDelta, mole_fractions, _tau, _delta, _reducing.T, _reducing.rhomolar);
}
long double HelmholtzEOSMixtureBackend::calc_d2alpha0_dDelta2(void)
{
    const int nTau = 0, nDelta = 2;
    return calc_alpha0_deriv_nocache(nTau, nDelta, mole_fractions, _tau, _delta, _reducing.T, _reducing.rhomolar);
}
long double HelmholtzEOSMixtureBackend::calc_d2alpha0_dDelta_dTau(void)
{
    const int nTau = 1, nDelta = 1;
    return calc_alpha0_deriv_nocache(nTau, nDelta, mole_fractions, _tau, _delta, _reducing.T, _reducing.rhomolar);
}
long double HelmholtzEOSMixtureBackend::calc_d2alpha0_dTau2(void)
{
    const int nTau = 2, nDelta = 0;
    return calc_alpha0_deriv_nocache(nTau, nDelta, mole_fractions, _tau, _delta, _reducing.T, _reducing.rhomolar);
}

} /* namespace CoolProp */
