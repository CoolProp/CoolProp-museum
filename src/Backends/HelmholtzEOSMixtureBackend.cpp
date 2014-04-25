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
#include "Solvers.h"
#include "MatrixMath.h"

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
        mole_fractions = std::vector<long double>(1, 1);
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
void HelmholtzEOSMixtureBackend::set_mole_fractions(const std::vector<long double> &mole_fractions)
{
    if (mole_fractions.size() != components.size())
    {
        throw ValueError(format("size of mole fraction vector [%d] does not equal that of component vector [%d]",mole_fractions.size(), components.size()));
    }
    this->mole_fractions = mole_fractions;
    this->K.resize(mole_fractions.size());
    this->lnK.resize(mole_fractions.size());
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

void HelmholtzEOSMixtureBackend::update_TP_guessrho(long double T, long double p, long double rho_guess)
{
    double rho = solver_rho_Tp(T, p, rho_guess);
    update(DmolarT_INPUTS, rho, T);
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
        case PQ_INPUTS:
        {
            _p = value1; _Q = value2;
            PQ_flash();
            break;
        }
        case PT_INPUTS:
        {
            _p = value1; _T = value2;
            PT_flash();
            break;
        }
        default:
        {
            throw ValueError(format("This pair of inputs [%s] is not yet supported", get_input_pair_short_desc(input_pair).c_str()));
        }
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
        if (_p > 1.05*components[0]->ancillaries.pL.evaluate(_T)){
            this->_phase = iphase_liquid; return;
		}
        else if (_p < 0.95*components[0]->ancillaries.pV.evaluate(_T)){
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
            // Pseudo-pure fluid
            long double rhoLanc, rhoVanc, rhoLsat, rhoVsat;
            long double psatLanc = components[0]->ancillaries.pL.evaluate(_T); // These ancillaries are used explicitly
		    long double psatVanc = components[0]->ancillaries.pV.evaluate(_T); // These ancillaries are used explicitly
		    try{
                rhoLanc = components[0]->ancillaries.rhoL.evaluate(_T);
                rhoVanc = components[0]->ancillaries.rhoV.evaluate(_T);

			    if (!ValidNumber(rhoLanc) || !ValidNumber(rhoVanc))
			    {
				    throw ValueError("pseudo-pure failed");
			    }

                rhoLsat = solver_rho_Tp(_T, psatLanc, rhoLanc);
			    rhoVsat = solver_rho_Tp(_T, psatVanc, rhoLanc);
			    if (!ValidNumber(rhoLsat) || !ValidNumber(rhoVsat) || 
				     fabs(rhoLsat/rhoLanc-1) > 0.1 || fabs(rhoVanc/rhoVsat-1) > 0.1)
			    {
				    throw ValueError("pseudo-pure failed");
			    }
		    }
		    catch (std::exception &){
			    // Near the critical point, the behavior is not very nice, so we will just use the ancillary near the critical point
			    rhoLsat = rhoLanc;
                rhoVsat = rhoVanc;
		    }
            _p = _Q*psatVanc + (1-_Q)*psatLanc;
            _rhomolar = 1/(_Q/rhoVsat + (1-_Q)/rhoLsat);
        }
    }
    else
    {
        // Set some imput options
        SaturationSolvers::mixture_VLE_IO options;
        options.sstype = SaturationSolvers::imposed_T;
        options.Nstep_max = 5;

        // Get an extremely rough guess by interpolation of ln(p) v. T curve where the limits are mole-fraction-weighted
        long double pguess = SaturationSolvers::saturation_preconditioner(this, _T, SaturationSolvers::imposed_T, mole_fractions);

        // Use Wilson iteration to obtain updated guess for pressure
        pguess = SaturationSolvers::saturation_Wilson(this, _Q, _T, SaturationSolvers::imposed_T, mole_fractions, pguess);
        
        // Actually call the successive substitution solver
        SaturationSolvers::successive_substitution(this, _Q, _T, pguess, mole_fractions, K, options);

        
    }
}
void HelmholtzEOSMixtureBackend::PQ_flash()
{
    if (is_pure_or_pseudopure)
    {
        if (!(components[0]->pEOS->pseudo_pure))
        {
            throw NotImplementedError();
            //// Set some imput options
            //SaturationSolvers::saturation_p_pure_Akasaka_options options;
            //options.omega = 1.0;
            //options.use_guesses = false;
            //// Actually call the solver
            //SaturationSolvers::saturation_p_pure_Akasaka(this, _T, options);
            //// Load the outputs
            //_p = _Q*SatV->p() + (1-_Q)*SatL->p();
            //_rhomolar = 1/(_Q/SatV->rhomolar() + (1-_Q)/SatL->rhomolar());
        }
        else{
            throw NotImplementedError();
        }
    }
    else
    {
        // Set some imput options
        SaturationSolvers::mixture_VLE_IO io;
        io.sstype = SaturationSolvers::imposed_p;
        io.Nstep_max = 20;

        // Get an extremely rough guess by interpolation of ln(p) v. T curve where the limits are mole-fraction-weighted
        long double Tguess = SaturationSolvers::saturation_preconditioner(this, _p, SaturationSolvers::imposed_p, mole_fractions);

        // Use Wilson iteration to obtain updated guess for temperature
        Tguess = SaturationSolvers::saturation_Wilson(this, _Q, _p, SaturationSolvers::imposed_p, mole_fractions, Tguess);
        
        // Actually call the successive substitution solver
        SaturationSolvers::successive_substitution(this, _Q, Tguess, _p, mole_fractions, K, io);
        
        PhaseEnvelope::PhaseEnvelope_GV ENV_GV;
        ENV_GV.build(this,mole_fractions,K,io);
    }
}
void HelmholtzEOSMixtureBackend::DmolarT_flash()
{
    if (imposed_phase_index > -1) 
    {
        // Use the phase defined by the imposed phase
        _phase = imposed_phase_index;
    }
    else
    {
        if (is_pure_or_pseudopure)
        {
            // Find the phase, while updating all internal variables possible
            DmolarT_phase_determination_pure_or_pseudopure();
        }
        else
        {
            _phase = iphase_gas;
            // Find the phase, while updating all internal variables possible
            //DmolarT_phase_determination_pure_or_pseudopure();
        }
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
    SimpleState reducing = calc_reducing_state_nocache(mole_fractions);
    long double delta = rhomolar/reducing.rhomolar;
    long double tau = reducing.T/T;
    
    // Calculate derivative if needed
    int nTau = 0, nDelta = 1;
    long double dalphar_dDelta = calc_alphar_deriv_nocache(nTau, nDelta, mole_fractions, tau, delta);

    // Get pressure
    return rhomolar*gas_constant()*T*(1+delta*dalphar_dDelta);
}
long double HelmholtzEOSMixtureBackend::solver_rho_Tp(long double T, long double p, long double rhomolar_guess)
{
    int phase;
    if (imposed_phase_index > -1)
        phase = imposed_phase_index;
    else
        phase = _phase;
    if (rhomolar_guess < 0){
        rhomolar_guess = solver_rho_Tp_SRK(T, p, phase);
    }
    
    class solver_TP_resid : public FuncWrapper1D
    {
    public:
	    double T, p;
	    HelmholtzEOSMixtureBackend *HEOS;

	    solver_TP_resid(HelmholtzEOSMixtureBackend *HEOS, long double T, long double p){ 
            this->HEOS = HEOS; this->T = T; this->p = p;
	    };
	    double call(double rhomolar){ 
            return (HEOS->p_rhoT(rhomolar,T)-p)/p;
        };
    };
    solver_TP_resid resid(this,T,p);
    std::string errstring;
    double rhomolar = Secant(resid, rhomolar_guess, 0.0001*rhomolar_guess, 1e-12, 100, errstring);
    return rhomolar;
}
long double HelmholtzEOSMixtureBackend::solver_rho_Tp_SRK(long double T, long double p, int phase)
{
    long double rhomolar, R_u = gas_constant(), a = 0, b = 0, k_ij = 0;

    for (std::size_t i = 0; i < components.size(); ++i)
    {
        long double Tci = components[i]->pEOS->reduce.T, pci = components[i]->pEOS->reduce.p, accentric_i = components[i]->pEOS->accentric;
        long double m_i = 0.480+1.574*accentric_i-0.176*pow(accentric_i, 2);
        long double b_i = 0.08664*R_u*Tci/pci;
        b += mole_fractions[i]*b_i;

        long double a_i = 0.42747*pow(R_u*Tci,2)/pci*pow(1+m_i*(1-sqrt(T/Tci)),2);

        for (std::size_t j = 0; j < components.size(); ++j)
        {
            long double Tcj = components[j]->pEOS->reduce.T, pcj = components[j]->pEOS->reduce.p, accentric_j = components[j]->pEOS->accentric;
            long double m_j = 0.480+1.574*accentric_j-0.176*pow(accentric_j, 2);

            long double a_j = 0.42747*pow(R_u*Tcj,2)/pcj*pow(1+m_j*(1-sqrt(T/Tcj)),2);
            
            if (i == j){
                k_ij = 0;
            }
            else{
                k_ij = 0;
            }

            a += mole_fractions[i]*mole_fractions[j]*sqrt(a_i*a_j)*(1-k_ij);
        }
    }

    long double A = a*p/pow(R_u*T,2);
    long double B = b*p/(R_u*T);

    //Solve the cubic for solutions for Z = p/(rho*R*T)
    double Z0, Z1, Z2; int Nsolns;
    solve_cubic(1, -1, A-B-B*B, -A*B, Nsolns, Z0, Z1, Z2);

    // Determine the guess value
    if (Nsolns == 1){
        rhomolar = p/(Z0*R_u*T);
    }
    else{
        long double rhomolar0 = p/(Z0*R_u*T);
        long double rhomolar1 = p/(Z1*R_u*T);
        long double rhomolar2 = p/(Z2*R_u*T);
        switch(phase)
        {
        case iphase_liquid:
            rhomolar = max3(rhomolar0, rhomolar1, rhomolar2); break;
        case iphase_gas:
            rhomolar = min3(rhomolar0, rhomolar1, rhomolar2); break;
        default:
            throw ValueError("Bad phase to solver_rho_Tp_SRK");
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
        solver_rho_Tp(_T, _p);
    }
}
long double HelmholtzEOSMixtureBackend::calc_pressure(void)
{    
    // Calculate the reducing parameters
    _delta = _rhomolar/_reducing.rhomolar;
    _tau = _reducing.T/_T;

    // Calculate derivative if needed
    long double dar_dDelta = dalphar_dDelta();
    long double R_u = gas_constant();

    // Get pressure
    _p = _rhomolar*R_u*_T*(1+_delta*dar_dDelta);

    //std::cout << format("p: %13.12f %13.12f %10.9f %10.9f %10.9f %10.9f %g\n",_T,_rhomolar,_tau,_delta,mole_fractions[0],dar_dDelta,_p);
    if (_p < 0){
        throw ValueError("Pressure is less than zero");
    }

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
    long double R_u = gas_constant();

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
    long double R_u = gas_constant();

    // Get molar entropy
    _smolar = R_u*(_tau*(da0_dTau+dar_dTau) - a0 - ar);

    return static_cast<long double>(_smolar);
}
long double HelmholtzEOSMixtureBackend::calc_umolar(void)
{
    // Calculate the reducing parameters
    _delta = _rhomolar/_reducing.rhomolar;
    _tau = _reducing.T/_T;

    // Calculate derivatives if needed, or just use cached values
    long double da0_dTau = dalpha0_dTau();
    long double dar_dTau = dalphar_dTau();
    long double dar_dDelta = dalphar_dDelta();
    long double R_u = gas_constant();

    // Get molar internal energy
    _umolar = R_u*_T*_tau*(da0_dTau+dar_dTau);

    return static_cast<long double>(_umolar);
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

SimpleState HelmholtzEOSMixtureBackend::calc_reducing_state_nocache(const std::vector<long double> & mole_fractions)
{
    SimpleState reducing;
    if (is_pure_or_pseudopure){
        reducing = components[0]->pEOS->reduce;
        
    }
    else{
        reducing.T = Reducing.p->Tr(mole_fractions);
        reducing.rhomolar = Reducing.p->rhormolar(mole_fractions);
    }
    return reducing;
}
void HelmholtzEOSMixtureBackend::calc_reducing_state(void)
{
    _reducing = calc_reducing_state_nocache(mole_fractions);
    _crit = _reducing;
}
long double HelmholtzEOSMixtureBackend::calc_alphar_deriv_nocache(const int nTau, const int nDelta, const std::vector<long double> &mole_fractions, const long double &tau, const long double &delta)
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
long double HelmholtzEOSMixtureBackend::calc_alpha0_deriv_nocache(const int nTau, const int nDelta, const std::vector<long double> &mole_fractions, 
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
long double HelmholtzEOSMixtureBackend::mixderiv_dln_fugacity_coefficient_dT__constrho_n(int i)
{
	double dtau_dT = -_tau/_T; //[1/K]
	return (dalphar_dTau() + mixderiv_d_ndalphardni_dTau(i)-1/(1+_delta*dalphar_dDelta())*(_delta*d2alphar_dDelta_dTau()))*dtau_dT;
}
long double HelmholtzEOSMixtureBackend::mixderiv_dln_fugacity_coefficient_drho__constT_n(int i)
{
    double ddelta_drho = 1/_reducing.rhomolar; //[m^3/mol]
	return (dalphar_dDelta() + mixderiv_d_ndalphardni_dDelta(i)-1/(1+_delta*dalphar_dDelta())*(_delta*d2alphar_dDelta2()+dalphar_dDelta()))*ddelta_drho;
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
long double HelmholtzEOSMixtureBackend::mixderiv_dpdrho__constT_n()
{
    long double R_u = static_cast<long double>(_gas_constant);
	return R_u*_T*(1+2*_delta*dalphar_dDelta()+pow(_delta,2)*d2alphar_dDelta2());
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
	return mixderiv_nd2nalphardnidnj__constT_V(j, i) + 1 - mixderiv_partial_molar_volume(j)/(R_u*_T)*mixderiv_ndpdni__constT_V_nj(i);
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
    const std::vector<long double> & mole_fractions = HEOS->get_mole_fractions();

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

void SaturationSolvers::x_and_y_from_K(long double beta, const std::vector<long double> &K, const std::vector<long double> &z, std::vector<long double> &x, std::vector<long double> &y)
{
	for (unsigned int i=0; i < K.size(); i++)
	{
		double denominator = (1-beta+beta*K[i]); // Common denominator
		x[i] = z[i]/denominator;
		y[i] = K[i]*z[i]/denominator;
	}
	//normalize_vector(x);
	//normalize_vector(y);
}

long double SaturationSolvers::successive_substitution(HelmholtzEOSMixtureBackend *HEOS, const long double beta, long double T, long double p, const std::vector<long double> &z, 
                                                       std::vector<long double> &K, mixture_VLE_IO &options)
{
	int iter = 1;
	long double change, f, df, deriv_liq, deriv_vap;
	std::size_t N = z.size();
    std::vector<long double> x, y, ln_phi_liq, ln_phi_vap;
	ln_phi_liq.resize(N); ln_phi_vap.resize(N); x.resize(N); y.resize(N);

	x_and_y_from_K(beta, K, z, x, y);
    HelmholtzEOSMixtureBackend *SatL = new HelmholtzEOSMixtureBackend(HEOS->get_components()), 
                               *SatV = new HelmholtzEOSMixtureBackend(HEOS->get_components());
    SatL->specify_phase(iphase_liquid);
    SatV->specify_phase(iphase_gas);

    SatL->set_mole_fractions(x);
    SatV->set_mole_fractions(y);
    long double rhomolar_liq = SatL->solver_rho_Tp_SRK(T, p, iphase_liquid); // [mol/m^3]
	long double rhomolar_vap = SatV->solver_rho_Tp_SRK(T, p, iphase_gas); // [mol/m^3]
    
	do
	{
        
        SatL->update_TP_guessrho(T, p, rhomolar_liq);
        SatV->update_TP_guessrho(T, p, rhomolar_vap);

		f = 0;
		df = 0;

		for (std::size_t i = 0; i < N; ++i)
		{
			ln_phi_liq[i] = SatL->mixderiv_ln_fugacity_coefficient(i);
            ln_phi_vap[i] = SatV->mixderiv_ln_fugacity_coefficient(i);

            if (options.sstype == imposed_p){
			    deriv_liq = SatL->mixderiv_dln_fugacity_coefficient_dT__constp_n(i);
                deriv_vap = SatV->mixderiv_dln_fugacity_coefficient_dT__constp_n(i);
            }
            else if (options.sstype == imposed_T){
                deriv_liq = SatL->mixderiv_dln_fugacity_coefficient_dp__constT_n(i);
                deriv_vap = SatV->mixderiv_dln_fugacity_coefficient_dp__constT_n(i);
            }
            else {throw ValueError();}

			K[i] = exp(ln_phi_liq[i]-ln_phi_vap[i]);
			
			f += z[i]*(K[i]-1)/(1-beta+beta*K[i]);

			double dfdK = K[i]*z[i]/pow(1-beta+beta*K[i],(int)2);
			df += dfdK*(deriv_liq-deriv_vap);
		}
		
		change = -f/df;

		if (options.sstype == imposed_p){
			T += change;
        }
        else if (options.sstype == imposed_T){
            p += change;
        }

		x_and_y_from_K(beta, K, z, x, y);
        SatL->set_mole_fractions(x);
        SatV->set_mole_fractions(y);

		iter += 1;
		if (iter > 50)
		{
			return _HUGE;
			//throw ValueError(format("saturation_p was unable to reach a solution within 50 iterations"));
		}
        rhomolar_liq = SatL->rhomolar();
        rhomolar_vap = SatV->rhomolar();
	}
	while(fabs(f) > 1e-12 || iter < options.Nstep_max);

    SatL->update_TP_guessrho(T, p, rhomolar_liq);
    SatV->update_TP_guessrho(T, p, rhomolar_vap);
	
    double pL = SatL->calc_pressure();
    double pV = SatV->calc_pressure();


    options.rhomolar_liq = SatL->rhomolar();
    options.rhomolar_vap = SatV->rhomolar();
    options.p = pL;
    options.T = T;
    options.x = SatL->get_mole_fractions();
    options.y = SatV->get_mole_fractions();

    delete SatL; delete SatV;
}

void SaturationSolvers::newton_raphson_VLE_GV::resize(unsigned int N)
{
	this->N = N;
	x.resize(N); 
	y.resize(N); 
	phi_ij_liq.resize(N); 
	phi_ij_vap.resize(N);
    dlnphi_drho_liq.resize(N);
    dlnphi_drho_vap.resize(N);

	r.resize(N+2);
    negative_r.resize(N+2);
	J.resize(N+2, std::vector<long double>(N+2, 0));

	neg_dFdS.resize(N+2);
	dXdS.resize(N+2);

	// Fill the vector -dFdS with zeros (Gerg Eqn. 7.132)
	std::fill(neg_dFdS.begin(), neg_dFdS.end(), (double)0.0);
	// Last entry is 1
	neg_dFdS[N+1] = 1.0;
}
void SaturationSolvers::newton_raphson_VLE_GV::check_Jacobian(HelmholtzEOSMixtureBackend *HEOS, const std::vector<long double> &z, std::vector<long double> &K, mixture_VLE_IO &IO)
{
    // Reset all the variables and resize
	pre_call();
    std::size_t N = K.size();
	resize(N);

    SatL = new HelmholtzEOSMixtureBackend(HEOS->get_components()), 
    SatV = new HelmholtzEOSMixtureBackend(HEOS->get_components());
    SatL->specify_phase(iphase_liquid);
    SatV->specify_phase(iphase_gas);

    long double rhomolar_liq0 = IO.rhomolar_liq;
    const long double rhomolar_vap0 = IO.rhomolar_vap;
    long double T0 = IO.T;
    long double beta0 = IO.beta;

    // Build the Jacobian and residual vectors for given inputs of K_i,T,p
	build_arrays(HEOS,beta0,T0,rhomolar_liq0,rhomolar_vap0,z,K);
    
    // Make copies of the base
    std::vector<long double> r0 = r;
    STLMatrix J0 = J;
    STLMatrix Jnum = J;

    for (std::size_t i = 0; i < N+2; ++i)
    {
        for (std::size_t j = 0; j < N+2; ++j)
        {
            Jnum[i][j] = _HUGE;
        }    
    }

    for (std::size_t j = 0; j < N; ++j)
    {
        std::vector<long double> KK = K;
        KK[j] += 1e-6;
        build_arrays(HEOS,beta0,T0,rhomolar_liq0,rhomolar_vap0,z,KK);
        std::vector<long double> r1 = r;
        for (std::size_t i = 0; i < N+2; ++i)
        {
            Jnum[i][j] = (r1[i]-r0[i])/(log(KK[j])-log(K[j]));
        }
        std::cout << vec_to_string(get_col(Jnum,j),"%12.11f") << std::endl;
        std::cout << vec_to_string(get_col(J,j),"%12.11f") << std::endl;
    }
    
	build_arrays(HEOS,beta0,T0+1e-6,rhomolar_liq0,rhomolar_vap0,z,K);
    std::vector<long double> r1 = r, JN = r;
    for (std::size_t i = 0; i < r.size(); ++i)
    {
        Jnum[i][N] = (r1[i]-r0[i])/(log(T0+1e-6)-log(T0));
    }
    std::cout << vec_to_string(get_col(Jnum,N),"%12.11f") << std::endl;
    std::cout << vec_to_string(get_col(J,N),"%12.11f") << std::endl;

    // Build the Jacobian and residual vectors for given inputs of K_i,T,p
	build_arrays(HEOS,beta0,T0,rhomolar_liq0+1e-3,rhomolar_vap0,z,K);
    std::vector<long double> r2 = r, JNp1 = r;
    for (std::size_t i = 0; i < r.size(); ++i)
    {
        Jnum[i][N+1] = (r2[i]-r0[i])/(log(rhomolar_liq0+1e-3)-log(rhomolar_liq0));
    }
    std::cout << vec_to_string(get_col(Jnum, N+1),"%12.11f") << std::endl;
    std::cout << vec_to_string(get_col(J,N+1),"%12.11f") << std::endl;

    delete SatL; delete SatV;
}
void SaturationSolvers::newton_raphson_VLE_GV::call(HelmholtzEOSMixtureBackend *HEOS, const std::vector<long double> &z, std::vector<long double> &K, mixture_VLE_IO &IO)
{
	int iter = 0;

	// Reset all the variables and resize
	pre_call();
	resize(K.size());

    SatL = new HelmholtzEOSMixtureBackend(HEOS->get_components()), 
    SatV = new HelmholtzEOSMixtureBackend(HEOS->get_components());
    SatL->specify_phase(iphase_liquid); // So it will always just use single-phase solution 
    SatV->specify_phase(iphase_gas); // So it will always just use single-phase solution

	do
	{
		// Build the Jacobian and residual vectors for given inputs of K_i,T,p
		build_arrays(HEOS,IO.beta,IO.T,IO.rhomolar_liq,IO.rhomolar_vap,z,K);

		// Solve for the step; v is the step with the contents 
		// [delta(lnK0), delta(lnK1), ..., delta(lnT), delta(lnrho')]
		std::vector<long double> v = linsolve(J, negative_r);

        max_rel_change = max_abs_value(v);

		// Set the variables again, the same structure independent of the specified variable
		for (unsigned int i = 0; i < N; i++)
		{
			K[i] = exp(log(K[i]) + v[i]);
			if (!ValidNumber(K[i]))
			{
				throw ValueError(format("K[i] (%g) is invalid",K[i]).c_str());
			}
		}
		IO.T = exp(log(IO.T) + v[N]);
		IO.rhomolar_liq = exp(log(IO.rhomolar_liq) + v[N+1]);

		if (fabs(IO.T) > 1e6)
		{
			/*std::cout << "J = " << vec_to_string(J,"%16.15g");
			std::cout << "nr = " << vec_to_string(r,"%16.15g");*/
			throw ValueError("Temperature or p has bad value");
		}
		
		//std::cout << iter << " " << T << " " << p << " " << error_rms << std::endl;
		iter++;
	}
	while(this->error_rms > 1e-8 && max_rel_change > 1000*LDBL_EPSILON && iter < IO.Nstep_max);
    Nsteps = iter;
    IO.p = p;
    IO.x = x; // Mole fractions in liquid
    IO.y = y; // Mole fractions in vapor
}

void SaturationSolvers::newton_raphson_VLE_GV::build_arrays(HelmholtzEOSMixtureBackend *HEOS, long double beta, long double T, long double rhomolar_liq, const long double rhomolar_vap, const std::vector<long double> &z, std::vector<long double> &K)
{
	// Step 0:
	// --------
	// Calculate the mole fractions in liquid and vapor phases
	x_and_y_from_K(beta, K, z, x, y);

    // Set the mole fractions in the classes
    SatL->set_mole_fractions(x);
    SatV->set_mole_fractions(y);

    // Update the liquid and vapor classes
    SatL->update(DmolarT_INPUTS, rhomolar_liq, T);
    SatV->update(DmolarT_INPUTS, rhomolar_vap, T);

    // For diagnostic purposes calculate the pressures (no derivatives are evaluated)
	long double p_liq = SatL->p();
	long double p_vap = SatV->p();
    p = 0.5*(p_liq+p_vap);

	// Step 2:
	// -------
	// Build the residual vector and the Jacobian matrix

	// For the residuals F_i
	for (unsigned int i = 0; i < N; ++i)
	{
		long double ln_phi_liq = SatL->mixderiv_ln_fugacity_coefficient(i);
		long double phi_iT_liq = SatL->mixderiv_dln_fugacity_coefficient_dT__constrho_n(i);
        dlnphi_drho_liq[i] = SatL->mixderiv_dln_fugacity_coefficient_drho__constT_n(i);
		for (unsigned int j = 0; j < N; ++j)
		{
            // I think this is wrong.
            phi_ij_liq[j] = SatL->mixderiv_ndln_fugacity_coefficient_dnj__constT_p(i,j) + (SatL->mixderiv_partial_molar_volume(i)/(SatL->gas_constant()*T)-1/p)*SatL->mixderiv_ndpdni__constT_V_nj(i); // 7.126 from GERG monograph
		}

		long double ln_phi_vap = SatV->mixderiv_ln_fugacity_coefficient(i);
		long double phi_iT_vap = SatV->mixderiv_dln_fugacity_coefficient_dT__constrho_n(i);
        dlnphi_drho_vap[i] = SatV->mixderiv_dln_fugacity_coefficient_drho__constT_n(i);
		for (unsigned int j = 0; j < N; ++j)
		{
            // I think this is wrong.
			phi_ij_vap[j] = SatV->mixderiv_ndln_fugacity_coefficient_dnj__constT_p(i,j) + (SatV->mixderiv_partial_molar_volume(i)/(SatV->gas_constant()*T)-1/p)*SatV->mixderiv_ndpdni__constT_V_nj(i); ; // 7.126 from GERG monograph
		}
		
		r[i] = log(K[i]) + ln_phi_vap - ln_phi_liq;
        // dF_i/d(ln(K_j))
		for (unsigned int j = 0; j < N; ++j)
		{	
			J[i][j] = K[j]*z[j]/pow(1-beta+beta*K[j],(int)2)*((1-beta)*phi_ij_vap[j]+beta*phi_ij_liq[j])+Kronecker_delta(i,j);
		}
		// dF_{i}/d(ln(T))
		J[i][N] = T*(phi_iT_vap-phi_iT_liq);
		// dF_{i}/d(ln(rho'))
		J[i][N+1] = -rhomolar_liq*dlnphi_drho_liq[i];
	}

	double summer1 = 0;
	for (unsigned int i = 0; i < N; ++i)
	{
		// Although the definition of this term is given by 
		// y[i]-x[i], when x and y are normalized, you get 
		// the wrong values.  Why? No idea.
		summer1 += z[i]*(K[i]-1)/(1-beta+beta*K[i]); 
	}
	r[N] = summer1;

	// For the residual term F_{N}, only non-zero derivatives are with respect
	// to ln(K[i])
	for (unsigned int j = 0; j < N; ++j)
	{
		J[N][j] = K[j]*z[j]/pow(1-beta+beta*K[j],(int)2);
	}

    // For the residual term F_{N+1} = p'-p''
    r[N+1] = p_liq-p_vap;
    for (unsigned int j = 0; j < N; ++j)
	{	
        J[N+1][j] = HEOS->gas_constant()*T*K[j]*z[j]/pow(1-beta+beta*K[j],(int)2)*((1-beta)*dlnphi_drho_vap[j]+beta*dlnphi_drho_liq[j]);
	}
    // dF_{N+1}/d(ln(T))
    J[N+1][N] = T*(SatL->mixderiv_dpdT__constV_n() - SatV->mixderiv_dpdT__constV_n());
    // dF_{N+1}/d(ln(rho'))
    J[N+1][N+1] = rhomolar_liq*SatL->mixderiv_dpdrho__constT_n();

	// Flip all the signs of the entries in the residual vector since we are solving Jv = -r, not Jv=r
	// Also calculate the rms error of the residual vector at this step
	error_rms = 0;
	for (unsigned int i = 0; i < N+2; ++i)
	{
		negative_r[i] = -r[i];
		error_rms += r[i]*r[i]; // Sum the squares
	}
	error_rms = sqrt(error_rms); // Square-root (The R in RMS)
}

void PhaseEnvelope::PhaseEnvelope_GV::build(HelmholtzEOSMixtureBackend *HEOS, const std::vector<long double> &z, std::vector<long double> &K, SaturationSolvers::mixture_VLE_IO &IO)
{
    // Use the residual function based on ln(K_i), ln(T) and ln(rho') as independent variables.  rho'' is specified
    SaturationSolvers::newton_raphson_VLE_GV NRVLE;
    SaturationSolvers::mixture_VLE_IO IO_NRVLE = IO;
    bubble.resize(z.size());
    dew.resize(z.size());
        
    // HACK
    IO_NRVLE.beta = 1.0;
    IO_NRVLE.Nstep_max = 30;
    int iter = 0;
    long double factor = IO_NRVLE.rhomolar_vap*0.01;
    for (;;)
    {
        if (iter > 0){ IO_NRVLE.rhomolar_vap += factor;}
        double p = IO_NRVLE.p;
        if (iter == 2 || (factor > 2 && factor < 0.24))
        {
            long double x = log(IO_NRVLE.rhomolar_vap);
            IO_NRVLE.T = exp(LinearInterp(dew.lnrhomolar_vap,dew.lnT,iter-2,iter-1,x));
            IO_NRVLE.rhomolar_liq = exp(LinearInterp(dew.lnrhomolar_vap,dew.lnrhomolar_liq,iter-2,iter-1,x));
            for (std::size_t i = 0; i < K.size(); ++i)
            {
                K[i] = exp(LinearInterp(dew.lnrhomolar_vap,dew.lnK[i],iter-2,iter-1,x));
            }
        }
        else if (iter == 3)
        {
            long double x = log(IO_NRVLE.rhomolar_vap);
            IO_NRVLE.T = exp(QuadInterp(dew.lnrhomolar_vap,dew.lnT,iter-3,iter-2,iter-1,x));
            IO_NRVLE.rhomolar_liq = exp(QuadInterp(dew.lnrhomolar_vap,dew.lnrhomolar_liq,iter-3,iter-2,iter-1,x));
            for (std::size_t i = 0; i < K.size(); ++i)
            {
                K[i] = exp(QuadInterp(dew.lnrhomolar_vap,dew.lnK[i],iter-3,iter-2,iter-1,x));
            }
        }
        else if (iter > 3)
        {
            long double x = log(IO_NRVLE.rhomolar_vap);
            IO_NRVLE.T = exp(CubicInterp(dew.lnrhomolar_vap, dew.lnT, iter-4, iter-3, iter-2, iter-1, x));
            IO_NRVLE.rhomolar_liq = exp(CubicInterp(dew.lnrhomolar_vap, dew.lnrhomolar_liq, iter-4, iter-3, iter-2, iter-1, x));
            for (std::size_t i = 0; i < K.size(); ++i)
            {
                K[i] = exp(CubicInterp(dew.lnrhomolar_vap, dew.lnK[i], iter-4, iter-3, iter-2, iter-1, x));
            }
        }
        /*if (IO_NRVLE.T > 344)
        {
            NRVLE.check_Jacobian(HEOS,z,K,IO_NRVLE);
        }*/
        NRVLE.call(HEOS,z,K,IO_NRVLE);
        dew.store_variables(IO_NRVLE.T,IO_NRVLE.p,IO_NRVLE.rhomolar_liq,IO_NRVLE.rhomolar_vap,K,IO_NRVLE.x,IO_NRVLE.y);
        iter ++;
        std::cout << format("%g %g %g %g %g %d %g\n",IO_NRVLE.p,IO_NRVLE.rhomolar_liq,IO_NRVLE.rhomolar_vap,IO_NRVLE.T,K[0],NRVLE.Nsteps,factor);
        if (iter < 5){continue;}
        if (NRVLE.Nsteps > 10)
        {
            factor /= 5;
        }
        else if (NRVLE.Nsteps > 5)
        {
            factor /= 1.2;
        }
        else if (NRVLE.Nsteps <= 4)
        {
            factor *= 1.2;
        }
        // Min step is 0.1 mol/m^3
        factor = std::max(factor,static_cast<long double>(0.1));
    }
}

} /* namespace CoolProp */
