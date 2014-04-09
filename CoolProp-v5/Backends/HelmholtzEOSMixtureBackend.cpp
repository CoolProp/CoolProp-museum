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
#include "../Solvers.h"

namespace CoolProp {

HelmholtzEOSMixtureBackend::HelmholtzEOSMixtureBackend(std::vector<std::string> component_names, bool generate_SatL_and_SatV) {
    std::vector<CoolPropFluid*> components;
    components.resize(component_names.size());

    for (unsigned int i = 0; i < components.size(); ++i)
    {
        components[i] = &(get_library().get(component_names[i]));
    }

    /// Set the components and associated flags
    set_components(components, generate_SatL_and_SatV);
}
HelmholtzEOSMixtureBackend::HelmholtzEOSMixtureBackend(std::vector<CoolPropFluid*> components, bool generate_SatL_and_SatV) {

    /// Set the components and associated flags
    set_components(components, generate_SatL_and_SatV);
}
void HelmholtzEOSMixtureBackend::set_components(std::vector<CoolPropFluid*> components, bool generate_SatL_and_SatV) {

    // Copy the components
    this->components = components;

    if (components.size() == 1){
        is_pure_or_pseudopure = true;
        mole_fractions = std::vector<double>(1, 1);
    }
    else{
        is_pure_or_pseudopure = false;
    }

    // Set the excess Helmholtz energy if a mixture
    if (!is_pure_or_pseudopure)
    {
        // Set the reducing model
        set_reducing_function();
        set_excess_term();
    }

    imposed_phase_index = -1;

    // Top-level class can hold copies of the base saturation classes, 
    // saturation classes cannot hold copies of the saturation classes
    if (generate_SatL_and_SatV)
    {
        SatL = new HelmholtzEOSMixtureBackend(components, false);
        SatL->specify_phase(iphase_liquid);
        SatV = new HelmholtzEOSMixtureBackend(components, false);
        SatV->specify_phase(iphase_gas);
    }
    else
    {
        SatL = NULL; SatV = NULL;
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
    Reducing.set(ReducingFunction::factory(components));
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
        default:
            throw ValueError(format("This input pair index [%d] is invalid", input_pair));
    }
    // Check the values that must always be set
    if (!ValidNumber(_p)){ throw ValueError("p is not a valid number");}
    if (!ValidNumber(_T)){ throw ValueError("T is not a valid number");}
    if (!ValidNumber(_rhomolar)){ throw ValueError("rhomolar is not a valid number");}
    if (!ValidNumber(_Q)){ throw ValueError("Q is not a valid number");}
}
void HelmholtzEOSMixtureBackend::DmolarT_phase_determination_pure_or_pseudopure()
{
    if (_T < _crit.T)
	{
		// Start to think about the saturation stuff
		// First try to use the ancillary equations if you are far enough away
		// You know how accurate the ancillary equations are thanks to using CoolProp code to refit them
        if (_rhomolar < 0.95*components[0]->ancillaries.rhoV.evaluate(_T)){
            this->_phase = iphase_gas; return;
		}
        else if (_rhomolar > 1.05*components[0]->ancillaries.rhoL.evaluate(_T)){
			this->_phase = iphase_liquid; return;
		}
		else{
			// Actually have to use saturation information sadly
			// For the given temperature, find the saturation state
			// Run the saturation routines to determine the saturation densities and pressures
            HelmholtzEOSMixtureBackend HEOS(components);
            SaturationSolvers::saturation_T_pure_options options;
            SaturationSolvers::saturation_T_pure(&HEOS, _T, options);

            long double Q = (1/_rhomolar-1/HEOS.SatL->rhomolar())/(1/HEOS.SatV->rhomolar()-1/HEOS.SatL->rhomolar());
			if (Q < -100*DBL_EPSILON){
				this->_phase = iphase_liquid;
			}
			else if (Q > 1+100*DBL_EPSILON){
				this->_phase = iphase_gas; 
			}
			else{
                this->_phase = iphase_twophase; 
            }
            _Q = Q;
             // Load the outputs
            _p = _Q*HEOS.SatV->p() + (1-_Q)*HEOS.SatL->p();
            _rhomolar = 1/(_Q/HEOS.SatV->rhomolar() + (1-_Q)/HEOS.SatL->rhomolar());
            return;
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
        if (!(components[0]->pEOS->pseudo_pure))
        {
            // Set some imput options
            SaturationSolvers::saturation_T_pure_options options;
            options.omega = 1.0;
            options.use_guesses = false;
            // Actually call the solver
            SaturationSolvers::saturation_T_pure(this, _T, options);
            // Load the outputs
            _p = _Q*SatV->p() + (1-_Q)*SatL->p();
            _rhomolar = 1/(_Q/SatV->rhomolar() + (1-_Q)/SatL->rhomolar());
        }
        else{
            //saturation_T_pseudopure();
        }
    }
    else
    {

    }
}
//void HelmholtzEOSMixtureBackend::PQ_flash()
//{
//    if (is_pure_or_pseudopure)
//    {
//        if (!(components[0]->pEOS->pseudo_pure))
//        {
//            // Set some imput options
//            SaturationSolvers::saturation_p_pure_Akasaka_options options;
//            options.omega = 1.0;
//            options.use_guesses = false;
//            // Actually call the solver
//            SaturationSolvers::saturation_p_pure_Akasaka(this, _T, options);
//            // Load the outputs
//            _p = _Q*SatV->p() + (1-_Q)*SatL->p();
//            _rhomolar = 1/(_Q/SatV->rhomolar() + (1-_Q)/SatL->rhomolar());
//        }
//        else{
//            //saturation_p_pseudopure();
//        }
//    }
//    else
//    {
//
//    }
//}
void HelmholtzEOSMixtureBackend::DmolarT_flash()
{
    if (is_pure_or_pseudopure)
    {
        if (imposed_phase_index > -1) 
        {
            // Use the phase defined by the imposed phase
            _phase = imposed_phase_index;
        }
        else
        {
            // Find the phase, while updating all internal variables possible
            DmolarT_phase_determination_pure_or_pseudopure();
        }
    }
    else
    {
        // TODO: hard coded phase to assume homogeneous phase
        _phase  = iphase_gas;
    }

    if (isHomogeneousPhase() && !ValidNumber(_p))
    {
        calc_pressure();
        _Q = -1;
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

    return static_cast<long double>(_p);
}
long double HelmholtzEOSMixtureBackend::calc_hmolar(void)
{
    // Calculate the reducing parameters
    _delta = _rhomolar/_reducing.rhomolar;
    _tau = _reducing.T/_T;

    // Calculate derivatives if needed, or just use cached values
    long double da0_dTau = dalpha0_dTau();
    long double dar_dTau = dalphar_dTau();
    long double dar_dDelta = dalphar_dDelta();
    long double R_u = static_cast<long double>(_gas_constant);

    // Get molar enthalpy
    _hmolar = R_u*_T*(1 + _tau*(da0_dTau+dar_dTau) + _delta*dar_dDelta);

    return static_cast<long double>(_hmolar);
}
long double HelmholtzEOSMixtureBackend::calc_smolar(void)
{
    // Calculate the reducing parameters
    _delta = _rhomolar/_reducing.rhomolar;
    _tau = _reducing.T/_T;

    // Calculate derivatives if needed, or just use cached values
    long double da0_dTau = dalpha0_dTau();
    long double ar = alphar();
    long double a0 = alpha0();
    long double dar_dTau = dalphar_dTau();
    long double dar_dDelta = dalphar_dDelta();
    long double R_u = static_cast<long double>(_gas_constant);

    // Get molar entropy
    _smolar = R_u*(_tau*(da0_dTau+dar_dTau) - a0 - ar);

    return static_cast<long double>(_smolar);
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
long double HelmholtzEOSMixtureBackend::calc_speed_sound(void)
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
    long double R_u = static_cast<long double>(_gas_constant);
    long double mm = static_cast<long double>(_molar_mass);

    // Get speed of sound
    _speed_sound = sqrt(R_u*_T/mm*(1+2*_delta*dar_dDelta+pow(_delta,2)*d2ar_dDelta2 - pow(1+_delta*dar_dDelta-_delta*_tau*d2ar_dDelta_dTau,2)/(pow(_tau,2)*(d2ar_dTau2 + d2a0_dTau2))));

    return static_cast<double>(_speed_sound);
}

long double HelmholtzEOSMixtureBackend::calc_fugacity_coefficient(int i)
{
    return exp(mixderiv_ln_fugacity_coefficient(i));
}

void HelmholtzEOSMixtureBackend::calc_reducing_state_nocache(const std::vector<double> & mole_fractions)
{
    if (is_pure_or_pseudopure){
        _reducing = components[0]->pEOS->reduce;
        _crit = _reducing;
    }
    else{
        _reducing.T = Reducing.p->Tr(mole_fractions);
        _reducing.rhomolar = Reducing.p->rhormolar(mole_fractions);
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
        std::size_t N = mole_fractions.size();
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
        std::size_t N = mole_fractions.size();
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


long double HelmholtzEOSMixtureBackend::mixderiv_dalphar_dxi(int i)
{
    return components[i]->pEOS->baser(_tau, _delta) + Excess.dalphar_dxi(_tau, _delta, mole_fractions, i);
}
long double HelmholtzEOSMixtureBackend::mixderiv_d2alphar_dxi_dTau(int i)
{
    return components[i]->pEOS->dalphar_dTau(_tau, _delta) + Excess.d2alphar_dxi_dTau(_tau, _delta, mole_fractions, i);
}
long double HelmholtzEOSMixtureBackend::mixderiv_d2alphar_dxi_dDelta(int i)
{
    return components[i]->pEOS->dalphar_dDelta(_tau, _delta) + Excess.d2alphar_dxi_dDelta(_tau, _delta, mole_fractions, i);
}
long double HelmholtzEOSMixtureBackend::mixderiv_d2alphardxidxj(int i, int j)
{
    return 0                           + Excess.d2alphardxidxj(_tau, _delta, mole_fractions, i, j);
}

long double HelmholtzEOSMixtureBackend::mixderiv_ln_fugacity_coefficient(int i)
{
	return alphar() + mixderiv_ndalphar_dni__constT_V_nj(i)-log(1+_delta*dalphar_dDelta());
}
long double HelmholtzEOSMixtureBackend::mixderiv_dln_fugacity_coefficient_dT__constrho(int i)
{
	double dtau_dT = -_tau/_T; //[1/K]
	return (dalphar_dTau() + mixderiv_d_ndalphardni_dTau(i)-1/(1+_delta*dalphar_dDelta())*(_delta*d2alphar_dDelta_dTau()))*dtau_dT;
}
long double HelmholtzEOSMixtureBackend::mixderiv_dnalphar_dni__constT_V_nj(int i)
{
	// GERG Equation 7.42
	return alphar() + mixderiv_ndalphar_dni__constT_V_nj(i);
}
long double HelmholtzEOSMixtureBackend::mixderiv_d2nalphar_dni_dT(int i)
{
	return -_tau/_T*(dalphar_dTau() + mixderiv_d_ndalphardni_dTau(i));
}
long double HelmholtzEOSMixtureBackend::mixderiv_dln_fugacity_coefficient_dT__constp_n(int i)
{
	double T = _reducing.T/_tau;
    long double R_u = static_cast<long double>(_gas_constant);
	return mixderiv_d2nalphar_dni_dT(i) + 1/T-mixderiv_partial_molar_volume(i)/(R_u*T)*mixderiv_dpdT__constV_n();
}
long double HelmholtzEOSMixtureBackend::mixderiv_partial_molar_volume(int i)
{
	return -mixderiv_ndpdni__constT_V_nj(i)/mixderiv_ndpdV__constT_n();
}

long double HelmholtzEOSMixtureBackend::mixderiv_dln_fugacity_coefficient_dp__constT_n(int i)
{
	// GERG equation 7.30
    long double R_u = static_cast<long double>(_gas_constant);
	double partial_molar_volume = mixderiv_partial_molar_volume(i); // [m^3/mol]
	double term1 = partial_molar_volume/(R_u*_T); // m^3/mol/(N*m)*mol = m^2/N = 1/Pa
	double term2 = 1.0/p();
	return term1 - term2;
}

long double HelmholtzEOSMixtureBackend::mixderiv_dln_fugacity_coefficient_dxj__constT_p_xi(int i, int j)
{
	// Gernert 3.115
    long double R_u = static_cast<long double>(_gas_constant);
	// partial molar volume is -dpdn/dpdV, so need to flip the sign here
	return mixderiv_d2nalphar_dni_dxj__constT_V(i,j) - mixderiv_partial_molar_volume(i)/(R_u*_T)*mixderiv_dpdxj__constT_V_xi(j);
}
long double HelmholtzEOSMixtureBackend::mixderiv_dpdxj__constT_V_xi(int j)
{
	// Gernert 3.130
    long double R_u = static_cast<long double>(_gas_constant);
	return _rhomolar*R_u*_T*(mixderiv_ddelta_dxj__constT_V_xi(j)*dalphar_dDelta()+_delta*mixderiv_d_dalpharddelta_dxj__constT_V_xi(j));
}

long double HelmholtzEOSMixtureBackend::mixderiv_d_dalpharddelta_dxj__constT_V_xi(int j)
{
	// Gernert Equation 3.134 (Catch test provided)
	return d2alphar_dDelta2()*mixderiv_ddelta_dxj__constT_V_xi(j)
		 + d2alphar_dDelta_dTau()*mixderiv_dtau_dxj__constT_V_xi(j)
		 + mixderiv_d2alphar_dxi_dDelta(j);
}

long double HelmholtzEOSMixtureBackend::mixderiv_dalphar_dxj__constT_V_xi(int j)
{
	//Gernert 3.119 (Catch test provided)
	return dalphar_dDelta()*mixderiv_ddelta_dxj__constT_V_xi(j)+dalphar_dTau()*mixderiv_dtau_dxj__constT_V_xi(j)+mixderiv_dalphar_dxi(j);
}
long double HelmholtzEOSMixtureBackend::mixderiv_d_ndalphardni_dxj__constT_V_xi(int i, int j)
{
	// Gernert 3.118
	return mixderiv_d_ndalphardni_dxj__constdelta_tau_xi(i,j)
		  + mixderiv_ddelta_dxj__constT_V_xi(j)*mixderiv_d_ndalphardni_dDelta(i) 
		  + mixderiv_dtau_dxj__constT_V_xi(j)*mixderiv_d_ndalphardni_dTau(i);
}
long double HelmholtzEOSMixtureBackend::mixderiv_ddelta_dxj__constT_V_xi(int j)
{
	// Gernert 3.121 (Catch test provided)
    return -_delta/_reducing.rhomolar*Reducing.p->drhormolardxi__constxj(mole_fractions,j);
}
long double HelmholtzEOSMixtureBackend::mixderiv_dtau_dxj__constT_V_xi(int j)
{
	// Gernert 3.122 (Catch test provided)
	return 1/_T*Reducing.p->dTrdxi__constxj(mole_fractions,j);
}

long double HelmholtzEOSMixtureBackend::mixderiv_dpdT__constV_n()
{
    long double R_u = static_cast<long double>(_gas_constant);
	return _rhomolar*R_u*(1+_delta*dalphar_dDelta()-_delta*_tau*d2alphar_dDelta_dTau());
}
long double HelmholtzEOSMixtureBackend::mixderiv_ndpdV__constT_n()
{
    long double R_u = static_cast<long double>(_gas_constant);
	return -pow(_rhomolar,2)*R_u*_T*(1+2*_delta*dalphar_dDelta()+pow(_delta,2)*d2alphar_dDelta2());
}
long double HelmholtzEOSMixtureBackend::mixderiv_ndpdni__constT_V_nj(int i)
{
	// Eqn 7.64 and 7.63
    long double R_u = static_cast<long double>(_gas_constant);
	double ndrhorbar_dni__constnj = Reducing.p->ndrhorbardni__constnj(mole_fractions,i);
	double ndTr_dni__constnj = Reducing.p->ndTrdni__constnj(mole_fractions,i);
	double summer = 0;
    for (unsigned int k = 0; k < mole_fractions.size(); ++k)
	{
		summer += mole_fractions[k]*mixderiv_d2alphar_dxi_dDelta(k);
	}
    double nd2alphar_dni_dDelta = _delta*d2alphar_dDelta2()*(1-1/_reducing.rhomolar*ndrhorbar_dni__constnj)+_tau*d2alphar_dDelta_dTau()/_reducing.T*ndTr_dni__constnj+mixderiv_d2alphar_dxi_dDelta(i)-summer;
    return _rhomolar*R_u*_T*(1+_delta*dalphar_dDelta()*(2-1/_reducing.rhomolar*ndrhorbar_dni__constnj)+_delta*nd2alphar_dni_dDelta);
}

long double HelmholtzEOSMixtureBackend::mixderiv_ndalphar_dni__constT_V_nj(int i)
{
    double term1 = _delta*dalphar_dDelta()*(1-1/_reducing.rhomolar*Reducing.p->ndrhorbardni__constnj(mole_fractions,i));
	double term2 = _tau*dalphar_dTau()*(1/_reducing.T)*Reducing.p->ndTrdni__constnj(mole_fractions,i);

	double s = 0;
	for (unsigned int k = 0; k < mole_fractions.size(); k++)
	{
		s += mole_fractions[k]*mixderiv_dalphar_dxi(k);
	}
	double term3 = mixderiv_dalphar_dxi(i);
	return term1 + term2 + term3 - s;
}
long double HelmholtzEOSMixtureBackend::mixderiv_ndln_fugacity_coefficient_dnj__constT_p(int i, int j)
{
    long double R_u = static_cast<long double>(_gas_constant);
	return mixderiv_nd2nalphardnidnj__constT_V(j, i) + 1 - mixderiv_partial_molar_volume(i)/(R_u*_T)*mixderiv_ndpdni__constT_V_nj(j);
}
long double HelmholtzEOSMixtureBackend::mixderiv_nddeltadni__constT_V_nj(int i)
{
    return _delta-_delta/_reducing.rhomolar*Reducing.p->ndrhorbardni__constnj(mole_fractions, i);
}
long double HelmholtzEOSMixtureBackend::mixderiv_ndtaudni__constT_V_nj(int i)
{
	return _tau/_reducing.T*Reducing.p->ndTrdni__constnj(mole_fractions, i);
}
long double HelmholtzEOSMixtureBackend::mixderiv_d_ndalphardni_dxj__constdelta_tau_xi(int i, int j)
{
    double line1 = _delta*mixderiv_d2alphar_dxi_dDelta(j)*(1-1/_reducing.rhomolar*Reducing.p->ndrhorbardni__constnj(mole_fractions, i));
    double line2 = -_delta*dalphar_dDelta()*(1/_reducing.rhomolar)*(Reducing.p->d_ndrhorbardni_dxj__constxi(mole_fractions, i, j)-1/_reducing.rhomolar*Reducing.p->drhormolardxi__constxj(mole_fractions,j)*Reducing.p->ndrhorbardni__constnj(mole_fractions,i));
	double line3 = _tau*mixderiv_d2alphar_dxi_dTau(j)*(1/_reducing.T)*Reducing.p->ndTrdni__constnj(mole_fractions, i);
	double line4 = _tau*dalphar_dTau()*(1/_reducing.T)*(Reducing.p->d_ndTrdni_dxj__constxi(mole_fractions,i,j)-1/_reducing.T*Reducing.p->dTrdxi__constxj(mole_fractions,j)*Reducing.p->ndTrdni__constnj(mole_fractions, i));
	double s = 0;
	for (unsigned int m = 0; m < mole_fractions.size(); m++)
	{
		s += mole_fractions[m]*mixderiv_d2alphardxidxj(j,m);
	}
	double line5 = mixderiv_d2alphardxidxj(i,j)-mixderiv_dalphar_dxi(j)-s;
	return line1+line2+line3+line4+line5;
}
long double HelmholtzEOSMixtureBackend::mixderiv_nd2nalphardnidnj__constT_V(int i, int j)
{	
	double line0 = mixderiv_ndalphar_dni__constT_V_nj(j); // First term from 7.46
	double line1 = mixderiv_d_ndalphardni_dDelta(i)*mixderiv_nddeltadni__constT_V_nj(j);
	double line2 = mixderiv_d_ndalphardni_dTau(i)*mixderiv_ndtaudni__constT_V_nj(j);
	double summer = 0;
    for (unsigned int k = 0; k < mole_fractions.size(); k++)
	{
		summer += mole_fractions[k]*mixderiv_d_ndalphardni_dxj__constdelta_tau_xi(i, k);
	}
	double line3 = mixderiv_d_ndalphardni_dxj__constdelta_tau_xi(i, j)-summer;
	return line0 + line1 + line2 + line3;
}
long double HelmholtzEOSMixtureBackend::mixderiv_d_ndalphardni_dDelta(int i)
{
	// The first line
    double term1 = (_delta*d2alphar_dDelta2()+dalphar_dDelta())*(1-1/_reducing.rhomolar*Reducing.p->ndrhorbardni__constnj(mole_fractions, i));

	// The second line
	double term2 = _tau*d2alphar_dDelta_dTau()*(1/_reducing.T)*Reducing.p->ndTrdni__constnj(mole_fractions, i);

	// The third line
	double term3 = mixderiv_d2alphar_dxi_dDelta(i);
	for (unsigned int k = 0; k < mole_fractions.size(); k++)
	{
		term3 -= mole_fractions[k]*mixderiv_d2alphar_dxi_dDelta(k);
	}
	return term1 + term2 + term3;
}

long double HelmholtzEOSMixtureBackend::mixderiv_d_ndalphardni_dTau(int i)
{
	// The first line
    double term1 = _delta*d2alphar_dDelta_dTau()*(1-1/_reducing.rhomolar*Reducing.p->ndrhorbardni__constnj(mole_fractions, i));

	// The second line
    double term2 = (_tau*d2alphar_dTau2()+dalphar_dTau())*(1/_reducing.T)*Reducing.p->ndTrdni__constnj(mole_fractions, i);

	// The third line
	double term3 = mixderiv_d2alphar_dxi_dTau(i);
	for (unsigned int k = 0; k < mole_fractions.size(); k++)
	{
		term3 -= mole_fractions[k]*mixderiv_d2alphar_dxi_dTau(k);
	}
	return term1 + term2 + term3;
}

void SaturationSolvers::saturation_T_pure(HelmholtzEOSMixtureBackend *HEOS, long double T, saturation_T_pure_options &options)
{
    // Set some imput options
    SaturationSolvers::saturation_T_pure_Akasaka_options _options;
    _options.omega = 1.0;
    _options.use_guesses = false;
    // Actually call the solver
    SaturationSolvers::saturation_T_pure_Akasaka(HEOS, T, _options);
}
void SaturationSolvers::saturation_T_pure_Akasaka(HelmholtzEOSMixtureBackend *HEOS, long double T, saturation_T_pure_Akasaka_options &options)
{
    // Start with the method of Akasaka

    /*
	This function implements the method of Akasaka 

	R. Akasaka,"A Reliable and Useful Method to Determine the Saturation State from 
	Helmholtz Energy Equations of State", 
	Journal of Thermal Science and Technology v3 n3,2008

	Ancillary equations are used to get a sensible starting point
	*/
    
    HEOS->calc_reducing_state();
    const SimpleState & reduce = HEOS->get_reducing();
    long double R_u = HEOS->calc_gas_constant();
    HelmholtzEOSMixtureBackend *SatL = HEOS->SatL, *SatV = HEOS->SatV;
    const std::vector<double> & mole_fractions = HEOS->get_mole_fractions();

	long double rhoL,rhoV,JL,JV,KL,KV,dJL,dJV,dKL,dKV;
	long double DELTA, deltaL=0, deltaV=0, tau=0, error, PL, PV, stepL, stepV;
	int iter=0;
	// Use the density ancillary function as the starting point for the solver
    try
	{
        if (!options.use_guesses)
		{
			// If very close to the critical temp, evaluate the ancillaries for a slightly lower temperature
            if (T > 0.99*HEOS->get_reducing().T){
				rhoL = HEOS->get_components()[0]->ancillaries.rhoL.evaluate(T-1);
				rhoV = HEOS->get_components()[0]->ancillaries.rhoV.evaluate(T-1);
			}
			else{
				rhoL = HEOS->get_components()[0]->ancillaries.rhoL.evaluate(T);
				rhoV = HEOS->get_components()[0]->ancillaries.rhoV.evaluate(T);
			}
		}
		else
		{
            rhoL = options.rhoL;
			rhoV = options.rhoV;
		}

		deltaL = rhoL/reduce.rhomolar;
		deltaV = rhoV/reduce.rhomolar;
		tau = reduce.T/T;
	}
	catch(NotImplementedError &)
	{
		/*double Tc = crit.T;
		double pc = crit.p.Pa;
		double w = 6.67228479e-09*Tc*Tc*Tc-7.20464352e-06*Tc*Tc+3.16947758e-03*Tc-2.88760012e-01;
		double q = -6.08930221451*w -5.42477887222;
		double pt = exp(q*(Tc/T-1))*pc;*/

		//double rhoL = density_Tp_Soave(T, pt, 0), rhoV = density_Tp_Soave(T, pt, 1);

		//deltaL = rhoL/reduce.rhomolar;
		//deltaV = rhoV/reduce.rhomolar;
		//tau = reduce.T/T;
	}
	//if (get_debug_level()>5){
	//		std::cout << format("%s:%d: Akasaka guess values deltaL = %g deltaV = %g tau = %g\n",__FILE__,__LINE__,deltaL, deltaV, tau).c_str();
	//	}

	do{
		/*if (get_debug_level()>8){
			std::cout << format("%s:%d: right before the derivs with deltaL = %g deltaV = %g tau = %g\n",__FILE__,__LINE__,deltaL, deltaV, tau).c_str();
		}*/

		// Calculate once to save on calls to EOS
        SatL->update(DmolarT_INPUTS, rhoL, T);
        SatV->update(DmolarT_INPUTS, rhoV, T);
        double alpharL = SatL->alphar();
		double alpharV = SatV->alphar();
        double dalphar_ddeltaL = SatL->dalphar_dDelta();
		double dalphar_ddeltaV = SatV->dalphar_dDelta();
        double d2alphar_ddelta2L = SatL->d2alphar_dDelta2();
        double d2alphar_ddelta2V = SatV->d2alphar_dDelta2();
		
		JL = deltaL * (1 + deltaL*dalphar_ddeltaL);
		JV = deltaV * (1 + deltaV*dalphar_ddeltaV);
		KL = deltaL*dalphar_ddeltaL + alpharL + log(deltaL);
		KV = deltaV*dalphar_ddeltaV + alpharV + log(deltaV);

		PL = R_u*reduce.rhomolar*T*JL;
		PV = R_u*reduce.rhomolar*T*JV;
		
		// At low pressure, the magnitude of d2alphar_ddelta2L and d2alphar_ddelta2V are enormous, truncation problems arise for all the partials
		dJL = 1 + 2*deltaL*dalphar_ddeltaL + deltaL*deltaL*d2alphar_ddelta2L;
		dJV = 1 + 2*deltaV*dalphar_ddeltaV + deltaV*deltaV*d2alphar_ddelta2V;
		dKL = 2*dalphar_ddeltaL + deltaL*d2alphar_ddelta2L + 1/deltaL;
		dKV = 2*dalphar_ddeltaV + deltaV*d2alphar_ddelta2V + 1/deltaV;
		
		DELTA = dJV*dKL-dJL*dKV;

		error = fabs(KL-KV)+fabs(JL-JV);

		//  Get the predicted step
		stepL = options.omega/DELTA*( (KV-KL)*dJV-(JV-JL)*dKV);
		stepV = options.omega/DELTA*( (KV-KL)*dJL-(JV-JL)*dKL);
		
		if (deltaL+stepL > 1 && deltaV+stepV < 1 && deltaV+stepV > 0){
            deltaL += stepL;  deltaV += stepV; 
            rhoL = deltaL*reduce.rhomolar; rhoV = deltaV*reduce.rhomolar;
		}
		else{
			throw ValueError(format("rhosatPure_Akasaka failed"));
		}
		iter++;
		if (iter > 100){
			throw SolutionError(format("Akasaka solver did not converge after 100 iterations"));
		}
	}
	while (error > 1e-10 && fabs(stepL) > 10*DBL_EPSILON*fabs(stepL) && fabs(stepV) > 10*DBL_EPSILON*fabs(stepV));
}


} /* namespace CoolProp */
