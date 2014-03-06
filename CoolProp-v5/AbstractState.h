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

	bool isSinglePhase(void){
		return (this->_phase==iLiquid || this->_phase==iGas);
	}

	bool isTwoPhase(void){
		return (this->_phase==iTwoPhase);
	}

	bool checkTwoPhase(void){
		if (!this->isCompressibleFluid()){throw ValueError(ERR_NOT_A_TWO_PHASE_FLUID);}
		if (!this->isTwoPhase()&&!_forceTwoPhase){throw ValueError(ERR_NOT_A_TWO_PHASE_STATE);}
		return true;
	}

	bool checkSinglePhase(void){
		if (!this->isSinglePhase()||!_forceSinglePhase){throw ValueError(ERR_NOT_A_TWO_PHASE_FUNCTION);}
		return true;
	}

	/// Two important points
	SimpleState _critical, _reducing;

	/// Bulk values
	double _rho, _T, _p, _Q, _R, _tau, _delta;

	CachedElement _h, _s, _logp, _logrho;

	/// Smoothing values
	double _rhospline, _dsplinedp, _dsplinedh;

	/// Cached low-level elements for in-place calculation of other properties
	/// These values cannot be reconstructed from the TTSE data and therefore
	/// always require a call to the EOS, hence the caching mechanism here.
	CachedElement _alpha0, _dalpha0_dTau, _dalpha0_dDelta, _d2alpha0_dTau2, _d2alpha0_dDelta_dTau,
			_d2alpha0_dDelta2, _d3alpha0_dTau3, _d3alpha0_dDelta_dTau2, _d3alpha0_dDelta2_dTau,
			_d3alpha0_dDelta3, _alphar, _dalphar_dTau, _dalphar_dDelta, _d2alphar_dTau2, _d2alphar_dDelta_dTau,
			_d2alphar_dDelta2, _d3alphar_dTau3, _d3alphar_dDelta_dTau2, _d3alphar_dDelta2_dTau,
			_d3alphar_dDelta3;

	CachedElement _dalphar_dDelta_lim, _d2alphar_dDelta2_lim,
			_d2alphar_dDelta_dTau_lim, _d3alphar_dDelta2_dTau_lim;

public:
	virtual AbstractState();
	virtual ~AbstractState();

	bool clear();
	virtual bool update(long iInput1, double Value1, long iInput2, double Value2);

	// ----------------------------------------
	// Bulk properties - temperature and density are directly calculated every time
	// All other parameters are calculated on an as-needed basis
	// ----------------------------------------
	double T(void)  {return _T;};
	double rho(void){return _rho;};
	double p(void)  {return _p;};
	double Q(void)  {return _Q;};

	double h(void);
	double s(void);
	double cp(void);
	double cv(void);
	double speed_sound(void);
	double isothermal_compressibility(void);
	double isobaric_expansion_coefficient(void);


	// ----------------------------------------
	// Transport properties
	// ----------------------------------------
	virtual double viscosity(void) = 0;
	virtual double conductivity(void) = 0;
	virtual double surface_tension(void) = 0;

	// ----------------------------------------
	// Derivatives of properties
	// ----------------------------------------
	virtual double dvdp_constT(void);
	virtual double dvdT_constp(void);

	// Density
	virtual double drhodh_constp(void);
	virtual double drhodp_consth(void);
	virtual double drhodp_constT(void);
	virtual double drhodT_constp(void);
	virtual double d2rhodh2_constp(void);
	virtual double d2rhodhdp(void);
	virtual double d2rhodhdQ(void);
	virtual double d2rhodp2_constT(void);
	virtual double d2rhodpdQ(void);
	virtual double d2rhodT2_constp(void);
	virtual double d2rhodTdp(void);

	// Pressure
	virtual double dpdrho_consth(void);
	virtual double dpdrho_constT(void);
	virtual double dpdT_consth(void);
	virtual double dpdT_constrho(void);
	virtual double d2pdrho2_constT(void);
	virtual double d2pdrhodT(void);
	virtual double d2pdT2_constrho(void);

	// Enthalpy
	virtual double dhdp_constrho(void);
	virtual double dhdp_constT(void);
	virtual double dhdrho_constp(void);
	virtual double dhdrho_constT(void);
	virtual double dhdT_constp(void);
	virtual double dhdT_constrho(void);
	virtual double d2hdp2_constT(void);
	virtual double d2hdrho2_constT(void);
	virtual double d2hdrhodT(void);
	virtual double d2hdT2_constp(void);
	virtual double d2hdT2_constrho(void);
	virtual double d2hdTdp(void);

	// Entropy
	virtual double dsdp_constT(void);
	virtual double dsdrho_constp(void);
	virtual double dsdrho_constT(void);
	virtual double dsdT_constp(void);
	virtual double dsdT_constrho(void);
	virtual double d2sdp2_constT(void);
	virtual double d2sdrho2_constT(void);
	virtual double d2sdrhodT(void);
	virtual double d2sdT2_constp(void);
	virtual double d2sdT2_constrho(void);
	virtual double d2sdTdp(void);

	// Fundamental derivative of gas dynamics
	virtual double fundamental_derivative_of_gas_dynamics(void);
	virtual double d2pdv2_consts(void);

	// Other functions and derivatives
	virtual double A(void);
	virtual double B(void);
	virtual double C(void);
	virtual double Z(void);

	virtual double dAdT_constrho(void);
	virtual double dAdrho_constT(void);
	// TODO: Add constXX qualifier
	virtual double dBdT(void);
	virtual double dCdT(void);
	virtual double dZdDelta(void);
	virtual double dZdTau(void);


	// ----------------------------------------
	// Helmholtz energy and derivatives
	// ----------------------------------------
	virtual double alpha0(void) = 0;
	virtual double dalpha0_dDelta(void) = 0;
	virtual double dalpha0_dTau(void) = 0;
	virtual double d2alpha0_dDelta2(void) = 0;
	virtual double d2alpha0_dDelta_dTau(void) = 0;
	virtual double d2alpha0_dTau2(void) = 0;
	virtual double d3alpha0_dDelta3(void) = 0;
	virtual double d3alpha0_dDelta2_dTau(void) = 0;
	virtual double d3alpha0_dDelta_dTau2(void) = 0;
	virtual double d3alpha0_dTau3(void) = 0;

	virtual double alphar(void) = 0;
	virtual double dalphar_dDelta(void) = 0;
	virtual double dalphar_dTau(void) = 0;
	virtual double d2alphar_dDelta2(void) = 0;
	virtual double d2alphar_dDelta_dTau(void) = 0;
	virtual double d2alphar_dTau2(void) = 0;
	virtual double d3alphar_dDelta3(void) = 0;
	virtual double d3alphar_dDelta2_dTau(void) = 0;
	virtual double d3alphar_dDelta_dTau2(void) = 0;
	virtual double d3alphar_dTau3(void) = 0;

	virtual double dalphar_dDelta_lim(void) = 0;
	virtual double d2alphar_dDelta2_lim(void) = 0;
	virtual double d2alphar_dDelta_dTau_lim(void) = 0;
	virtual double d3alphar_dDelta2_dTau_lim(void) = 0;
};

} /* namespace CoolProp */
#endif /* ABSTRACTSTATE_H_ */
