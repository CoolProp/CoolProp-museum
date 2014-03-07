/*
 * AbstractState.cpp
 *
 *  Created on: 21 Dec 2013
 *      Author: jowr
 */

#include "math.h"
#include "AbstractState.h"

namespace CoolProp {

AbstractState::AbstractState() {
	// TODO Auto-generated constructor stub
}

AbstractState::~AbstractState() {
	// TODO Auto-generated destructor stub
}

bool AbstractState::clear() {
	// Reset all instances of CachedElement and overwrite
	// the internal double values with -_HUGE
	this->_fluid_type = FLUID_TYPE_UNDEFINED;
	this->_phase = iUnknown;
	this->_forceSinglePhase = false;
	this->_forceTwoPhase = false;

	this->_critical.T = -_HUGE;
	this->_critical.h = -_HUGE;
	this->_critical.p = -_HUGE;
	this->_critical.rho = -_HUGE;
	this->_critical.s = -_HUGE;

	this->_reducing.T = -_HUGE;
	this->_reducing.h = -_HUGE;
	this->_reducing.p = -_HUGE;
	this->_reducing.rho = -_HUGE;
	this->_reducing.s = -_HUGE;

	/// Bulk values
	this->_rhomolar = -_HUGE;
	this->_T = -_HUGE;
	this->_p = -_HUGE;
	this->_Q = -_HUGE;
	this->_tau = -_HUGE;
	this->_delta = -_HUGE;
	this->_hmolar.clear();
	this->_smolar.clear();
	this->_logp.clear();
	this->_logrhomolar.clear();

	///// Smoothing values
	//this->rhospline = -_HUGE;
	//this->dsplinedp = -_HUGE;
	//this->dsplinedh = -_HUGE;

	/// Cached low-level elements for in-place calculation of other properties
	this->_alpha0.clear();
	this->_dalpha0_dTau.clear();
	this->_dalpha0_dDelta.clear();
	this->_d2alpha0_dTau2.clear();
	this->_d2alpha0_dDelta_dTau.clear();
	this->_d2alpha0_dDelta2.clear();
	this->_d3alpha0_dTau3.clear();
	this->_d3alpha0_dDelta_dTau2.clear();
	this->_d3alpha0_dDelta2_dTau.clear();
	this->_d3alpha0_dDelta3.clear();
	this->_alphar.clear();
	this->_dalphar_dTau.clear();
	this->_dalphar_dDelta.clear();
	this->_d2alphar_dTau2.clear();
	this->_d2alphar_dDelta_dTau.clear();
	this->_d2alphar_dDelta2.clear();
	this->_d3alphar_dTau3.clear();
	this->_d3alphar_dDelta_dTau2.clear();
	this->_d3alphar_dDelta2_dTau.clear();
	this->_d3alphar_dDelta3.clear();

	this->_dalphar_dDelta_lim.clear();
	this->_d2alphar_dDelta2_lim.clear();
	this->_d2alphar_dDelta_dTau_lim.clear();
	this->_d3alphar_dDelta2_dTau_lim.clear();

	return true;
}

double AbstractState::hmolar(void){
	if (!_hmolar) _hmolar = calc_hmolar();
	return _hmolar;
}
double AbstractState::smolar(void){
	if (!_smolar) _smolar = calc_smolar();
	return _smolar;
}
double AbstractState::cpmolar(void){
	if (!_cpmolar) _cpmolar = calc_cpmolar();
	return _cpmolar;
}
double AbstractState::cvmolar(void){
	if (!_cvmolar) _cvmolar = calc_cvmolar();
	return _cvmolar;
}
double AbstractState::speed_sound(void){
	if (!_speed_sound) _speed_sound = calc_speed_sound();
	return _speed_sound;
}
double AbstractState::viscosity(void){
	if (!_viscosity) _viscosity = calc_viscosity();
	return _viscosity;
}
double AbstractState::conductivity(void){
	if (!_conductivity) _conductivity = calc_conductivity();
	return _conductivity;
}
double AbstractState::surface_tension(void){
	if (!_surface_tension) _surface_tension = calc_surface_tension();
	return _surface_tension;
}

//virtual double AbstractState::isothermal_compressibility(void){
//	return 1.0/(_rho*dpdrho_constT());
//}
//virtual double AbstractState::isobaric_expansion_coefficient(void){
//	return -1.0/(_rho*_rho)*drhodT_constp();
//}

//	// ----------------------------------------
//	// Smoothing functions for density
//	// ----------------------------------------
//	/// A smoothed version of the derivative using a spline curve in the region of x=0 to x=xend
//	virtual double AbstractState::drhodh_constp_smoothed(double xend);
//	/// A smoothed version of the derivative using a spline curve in the region of x=0 to x=xend
//	virtual double AbstractState::drhodp_consth_smoothed(double xend);
//	/// Density corresponding to the smoothed derivatives in the region of x=0 to x=xend
//	virtual void AbstractState::rho_smoothed(double xend, double *rho_spline, double *dsplinedh, double *dsplinedp);
//
//
//	// ----------------------------------------
//	// Transport properties // TODO: Fix it!
//	// ----------------------------------------
//
//	virtual double AbstractState::surface_tension(void);
//
//
	// ----------------------------------------
	// Derivatives of properties
	// ----------------------------------------
//	virtual double AbstractState::dvdp_constT(void){
//		return -1/(_rho*_rho)/dpdrho_constT();
//	}
//	virtual double AbstractState::dvdT_constp(void){
//		return -1/(_rho*_rho)*drhodT_constp();
//	}
//
//	// Density
//	virtual double AbstractState::drhodh_constp(void){
//		return 1.0/(dhdrho_constT()-dhdT_constrho()*dpdrho_constT()/dpdT_constrho());
//	}
//	virtual double AbstractState::drhodp_consth(void){
//		return 1.0/(dpdrho_constT()-dpdT_constrho()*dhdrho_constT()/dhdT_constrho());
//	}
//	virtual double AbstractState::drhodp_constT(void){
//		return -dpdT_constrho()/dpdrho_constT();
//	}
//	virtual double AbstractState::drhodT_constp(void){
//		return -dpdT_constrho()/dpdrho_constT();
//	}
//	virtual double AbstractState::d2rhodh2_constp(void){
//		//double A = dpdT_constrho()*dhdrho_constT()-dpdrho_constT()*dhdT_constrho();
//		//double dAdT_constrho = d2pdT2_constrho()*dhdrho_constT()+dpdT_constrho()*d2hdrhodT()-d2pdrhodT()*dhdT_constrho()-dpdrho_constT()*d2hdT2_constrho();
//		//double dAdrho_constT = d2pdrhodT()*dhdrho_constT()+dpdT_constrho()*d2hdrho2_constT()-d2pdrho2_constT()*dhdT_constrho()-dpdrho_constT()*d2hdrhodT();
//		double ddT_drhodh_p_constrho = 1.0/A()*d2pdT2_constrho()-1.0/(A()*A())*dAdT_constrho()*dpdT_constrho();
//		double ddrho_drhodh_p_constT = 1.0/A()*d2pdrhodT()      -1.0/(A()*A())*dAdrho_constT()*dpdT_constrho();
//		return ddT_drhodh_p_constrho/dhdT_constp()+ddrho_drhodh_p_constT/dhdrho_constp();
//	}
//	virtual double AbstractState::d2rhodhdp(void){
//		//double A = dpdT_constrho()*dhdrho_constT()-dpdrho_constT()*dhdT_constrho();
//		//double dAdT_constrho = d2pdT2_constrho()*dhdrho_constT()+dpdT_constrho()*d2hdrhodT()-d2pdrhodT()*dhdT_constrho()-dpdrho_constT()*d2hdT2_constrho();
//		//double dAdrho_constT = d2pdrhodT()*dhdrho_constT()+dpdT_constrho()*d2hdrho2_constT()-d2pdrho2_constT()*dhdT_constrho()-dpdrho_constT()*d2hdrhodT();
//		double ddT_drhodp_h_constrho = -1.0/A()*d2hdT2_constrho()+1.0/(A()*A())*dAdT_constrho()*dhdT_constrho();
//		double ddrho_drhodp_h_constT = -1.0/A()*d2hdrhodT()      +1.0/(A()*A())*dAdrho_constT()*dhdT_constrho();
//		return ddT_drhodp_h_constrho/dhdT_constp()+ddrho_drhodp_h_constT/dhdrho_constp();
//	}
//	//virtual double AbstractState::d2rhodhdQ(void);//TODO: only two-phase, not here
//	virtual double AbstractState::d2rhodp2_constT(void){
//		return -d2pdrho2_constT()/pow(dpdrho_constT(),3);
//	}
//	//virtual double AbstractState::d2rhodpdQ(void);//TODO: only two-phase, not here
//	virtual double AbstractState::d2rhodT2_constp(void){
//		double ddrho_drhodT_p_constT = (dpdT_constrho()*d2pdrho2_constT()-dpdrho_constT()*d2pdrhodT())/pow(dpdrho_constT(),2);
//		double ddT_drhodT_p_constrho = (dpdT_constrho()*d2pdrhodT()-dpdrho_constT()*d2pdT2_constrho())/pow(dpdrho_constT(),2);
//		return ddT_drhodT_p_constrho+ddrho_drhodT_p_constT*drhodT_constp();
//	}
//	virtual double AbstractState::d2rhodTdp(void){
//		return (dpdT_constrho()*d2pdrho2_constT()-dpdrho_constT()*d2pdrhodT())/pow(dpdrho_constT(),3);
//	}
//
//	// Pressure
//	virtual double AbstractState::dpdrho_consth(void){
//		return dpdrho_constT() - dpdT_constrho()*dhdrho_constT()/dhdT_constrho();
//	}
//	virtual double AbstractState::dpdrho_constT(void){
//		return _R*_T*(1+2*delta*dalphar_dDelta()+delta*delta*d2alphar_dDelta2());
//	}
//	virtual double AbstractState::dpdT_consth(void){
//		return dpdT_constrho() - dpdrho_constT()*dhdT_constrho()/dhdrho_constT();
//	}
//	virtual double AbstractState::dpdT_constrho(void){
//		return _R*_rho*(1+delta*dalphar_dDelta()-delta*tau*d2alphar_dDelta_dTau());
//	}
//	virtual double AbstractState::d2pdrho2_constT(void){
//		return _R*_T/_rho*(2*delta*dalphar_dDelta()+4*delta*delta*d2alphar_dDelta2()+delta*delta*delta*d3alphar_dDelta3());
//	}
//	virtual double AbstractState::d2pdrhodT(void){
//		return _R*((1+2*delta*dalphar_dDelta()+delta*delta*d2alphar_dDelta2())+_T*(2*delta*d2alphar_dDelta_dTau()+delta*delta*d3alphar_dDelta2_dTau())*(-tau/_T));
//	}
//	virtual double AbstractState::d2pdT2_constrho(void){
//		return _R*_rho*delta*tau*tau/_T*d3alphar_dDelta_dTau2();
//	}
//
//	// Enthalpy
//	virtual double AbstractState::dhdp_constrho(void){
//		//	double dalphar_dDelta = pFluid->dalphar_dDelta(tau,delta);
//		//	double d2alphar_dDelta_dTau = pFluid->d2alphar_dDelta_dTau(tau,delta);
//		//	double d2alphar_dDelta2 = pFluid->d2alphar_dDelta2(tau,delta);
//		//	double dpdrho = R*T*(1+2*delta*dalphar_dDelta+delta*delta*d2alphar_dDelta2);
//		//	double dpdT = R*rho*(1+delta*dalphar_dDelta-delta*tau*d2alphar_dDelta_dTau);
//		//	double cp = -tau*tau*R*(pFluid->d2alpha0_dTau2(tau,delta)+pFluid->d2alphar_dTau2(tau,delta))+T/rho/rho*(dpdT*dpdT)/dpdrho;
//		double dpdrho = dpdrho_constT();
//		double dpdT   = dpdT_constrho();
//		double drhodT = -dpdT/dpdrho;
//		return -cp()/dpdrho/drhodT-_T*drhodT*(-1/_rho/_rho)+1/_rho;
//	}
//	virtual double AbstractState::dhdp_constT(void){
//		return dhdrho_constT()/dpdrho_constT();
//	}
//	virtual double AbstractState::dhdrho_constp(void){
//		return dhdrho_constT() - dhdT_constrho()*dpdrho_constT()/dpdT_constrho();
//	}
//	virtual double AbstractState::dhdrho_constT(void){
//		return _T*_R/_rho*(_tau*_delta*d2alphar_dDelta_dTau()+delta*dalphar_dDelta()+delta*delta*d2alphar_dDelta2());
//	}
//	virtual double AbstractState::dhdT_constp(void){
//		return dhdT_constrho() - dhdrho_constT()*dpdT_constrho()/dpdrho_constT();
//	}
//	virtual double AbstractState::dhdT_constrho(void){
//		return _R*(-tau*tau*(d2alpha0_dTau2()+d2alphar_dTau2())+1+_delta*dalphar_dDelta()-delta*tau*d2alphar_dDelta_dTau());
//	}
//	virtual double AbstractState::d2hdp2_constT(void){
//		return (d2hdrho2_constT()-dhdp_constT()*d2pdrho2_constT())/pow(dpdrho_constT(),2);
//	}
//	virtual double AbstractState::d2hdrho2_constT(void){
//		return _T*_R/_rho*(tau*delta*d3alphar_dDelta2_dTau()+tau*d2alphar_dDelta_dTau()+delta*d2alphar_dDelta2()+dalphar_dDelta()+delta*delta*d3alphar_dDelta3()+2*delta*d2alphar_dDelta2())/reducing.rho - dhdrho_constT()/_rho;
//	}
//	virtual double AbstractState::d2hdrhodT(void){
//		return _R*(-tau*tau*d3alphar_dDelta_dTau2()+delta*d2alphar_dDelta2()+dalphar_dDelta()-delta*tau*d3alphar_dDelta2_dTau()-tau*d2alphar_dDelta_dTau())/reducing.rho;
//	}
//	virtual double AbstractState::d2hdT2_constp(void){
//		double ddT_dhdT = d2hdT2_constrho()-1/pow(dpdrho_constT(),2)*(dpdrho_constT()*(dhdrho_constT()*d2pdT2_constrho()+d2hdrhodT()*dpdT_constrho())-dhdrho_constT()*dpdT_constrho()*d2pdrhodT());
//		double drho_dhdT = d2hdrhodT()-1/pow(dpdrho_constT(),2)*(dpdrho_constT()*(dhdrho_constT()*d2pdrhodT()+d2hdrho2_constT()*dpdT_constrho())-dhdrho_constT()*dpdT_constrho()*d2pdrho2_constT());
//		return ddT_dhdT-drho_dhdT*dpdT_constrho()/dpdrho_constT();
//	}
//	virtual double AbstractState::d2hdT2_constrho(void){
//		return _R*(-_tau*_tau*(d3alpha0_dTau3()+d3alphar_dTau3())-2*_tau*(d2alpha0_dTau2()+d2alphar_dTau2())-_delta*tau*d3alphar_dDelta_dTau2())*(-tau/_T);
//	}
//	virtual double AbstractState::d2hdTdp(void){
//		return 1.0/dpdrho_constT()*(d2hdrhodT()-dhdp_constT()*(drhodT_constp()*d2pdrho2_constT()+d2pdrhodT())+d2hdrho2_constT()*drhodT_constp());
//	}
//
//	// Entropy
//	virtual double AbstractState::dsdp_constT(void){
//		return dsdrho_constT()/dpdrho_constT();
//	}
//	virtual double AbstractState::dsdrho_constp(void){
//		return dsdrho_constT() - dsdT_constrho()*dpdrho_constT()/dpdT_constrho();
//	}
//	virtual double AbstractState::dsdrho_constT(void){
//		return -_R/_rho*(1+delta*dalphar_dDelta()-delta*tau*d2alphar_dDelta_dTau());
//	}
//	virtual double AbstractState::dsdT_constp(void){
//		return dsdT_constrho() - dsdrho_constT()*dpdT_constrho()/dpdrho_constT();
//	}
//	virtual double AbstractState::dsdT_constrho(void){
//		return -_R*tau*tau/_T*(d2alpha0_dTau2()+d2alphar_dTau2());
//	}
//	virtual double AbstractState::d2sdp2_constT(void){
//		return (d2sdrho2_constT()-dsdp_constT()*d2pdrho2_constT())/pow(dpdrho_constT(),2);
//	}
//	virtual double AbstractState::d2sdrho2_constT(void){
//		return -_R/_rho*(delta*d2alphar_dDelta2()+dalphar_dDelta()-tau*delta*d3alphar_dDelta2_dTau()-tau*d2alphar_dDelta_dTau())/reducing.rho+_R/_rho/_rho*(1+delta*dalphar_dDelta()-delta*tau*d2alphar_dDelta_dTau());
//	}
//	virtual double AbstractState::d2sdrhodT(void){
//		// d2alpha0_dDelta_dTau2() is zero by definition
//		return -_R*tau*tau/_T*d3alphar_dDelta_dTau2()/reducing.rho;
//	}
//	virtual double AbstractState::d2sdT2_constp(void){
//		double ddT_dsdT  = d2sdT2_constrho()-1/pow(dpdrho_constT(),2)*(dpdrho_constT()*(dsdrho_constT()*d2pdT2_constrho()+d2sdrhodT()*dpdT_constrho())-dsdrho_constT()*dpdT_constrho()*d2pdrhodT());
//		double drho_dsdT = d2sdrhodT()-1/pow(dpdrho_constT(),2)*(dpdrho_constT()*(dsdrho_constT()*d2pdrhodT()+d2sdrho2_constT()*dpdT_constrho())-dsdrho_constT()*dpdT_constrho()*d2pdrho2_constT());
//		return ddT_dsdT-drho_dsdT*dpdT_constrho()/dpdrho_constT();
//	}
//	virtual double AbstractState::d2sdT2_constrho(void){
//		return -_R/_T*(tau*tau*(d3alpha0_dTau3()+d3alphar_dTau3())+2*tau*(d2alpha0_dTau2()+d2alphar_dTau2()))*(-tau/_T)+_R*tau*tau/_T/_T*(d2alpha0_dTau2()+d2alphar_dTau2());
//	}
//	virtual double AbstractState::d2sdTdp(void){
//		return 1.0/dpdrho_constT()*(d2sdrhodT()-dsdp_constT()*(drhodT_constp()*d2pdrho2_constT()+d2pdrhodT())+d2sdrho2_constT()*drhodT_constp());
//	}
//
//	// Fundamental derivative of gas dynamics
//	virtual double AbstractState::fundamental_derivative_of_gas_dynamics(void){
//		return d2pdv2_consts()/pow(speed_sound(),2)/2.0/pow(_rho,3);
//	}
//	virtual double AbstractState::d2pdv2_consts(void){
//		double cv = this->cv();
//		// Convert each derivative in terms of volume rather than density
//		// Note drhodv = -rho^2
//		double d2pdv2_constT = _rho*_rho*_rho*_rho*d2pdrho2_constT()+2.0*_rho*_rho*_rho*dpdrho_constT();
//		double dpdT_constv = dpdT_constrho();
//		double d2pdvdT = -_rho*_rho*d2pdrhodT();
//		double d2pdT2_constv = d2pdT2_constrho();
//		double dcv_dT_constv = _R*tau/_T*(2.0*tau*(d2alpha0_dTau2()+d2alphar_dTau2())+tau*tau*(d3alpha0_dTau3()+d3alphar_dTau3()));
//		double LAMBDA1 = d2pdv2_constT;
//		double LAMBDA2 = -3.0*_T/cv*dpdT_constv*d2pdvdT;
//		double LAMBDA3 = +pow(_T/cv*dpdT_constv,2)*(3.0*d2pdT2_constv+1/_T*dpdT_constv*(1.0-_T/cv*dcv_dT_constv));
//		return LAMBDA1 + LAMBDA2 + LAMBDA3;
//	}
//
//	// Other functions and derivatives
//	virtual double AbstractState::A(void){
//		return dpdT_constrho()*dhdrho_constT()-dpdrho_constT()*dhdT_constrho();
//	}
//	virtual double AbstractState::B(void){
//		// given by B*rhoc=lim(delta --> 0) [dalphar_ddelta(tau)]
//		return 1.0/reducing.rho*dalphar_dDelta_lim();
//	}
//	virtual double AbstractState::C(void){
//		// given by C*rhoc^2=lim(delta --> 0) [d2alphar_dDelta2(tau)]
//		return 1.0/(reducing.rho*reducing.rho)*d2alphar_dDelta2_lim();
//	}
//	virtual double AbstractState::Z(void){
//		return 1+delta*dalphar_dDelta();
//	}
//	virtual double AbstractState::dAdT_constrho(void){
//		return d2pdT2_constrho()*dhdrho_constT()+dpdT_constrho()*d2hdrhodT()-d2pdrhodT()*dhdT_constrho()-dpdrho_constT()*d2hdT2_constrho();
//	}
//	virtual double AbstractState::dAdrho_constT(void){
//		return d2pdrhodT()*dhdrho_constT()+dpdT_constrho()*d2hdrho2_constT()-d2pdrho2_constT()*dhdT_constrho()-dpdrho_constT()*d2hdrhodT();
//	}
//	// TODO: Add constXX qualifier
//	virtual double AbstractState::dBdT(void){
//		return 1.0/reducing.rho*d2alphar_dDelta_dTau_lim()*-reducing.T/_T/_T;
//	}
//	virtual double AbstractState::dCdT(void){
//		return 1.0/(reducing.rho*reducing.rho)*d3alphar_dDelta2_dTau_lim()*-reducing.T/_T/_T;
//	}
//	virtual double AbstractState::dZdDelta(void){
//		return delta*d2alphar_dDelta2()+dalphar_dDelta();
//	}
//	virtual double AbstractState::dZdTau(void){
//		return delta*d2alphar_dDelta_dTau();
//	}
//
////	// ----------------------------------------
////	// Helmholtz energy and derivatives
////	// ----------------------------------------
////	virtual double AbstractState::alpha0(void);
////	virtual double AbstractState::dalpha0_dDelta(void);
////	virtual double AbstractState::dalpha0_dTau(void);
////	virtual double AbstractState::d2alpha0_dDelta2(void);
////	virtual double AbstractState::d2alpha0_dDelta_dTau(void);
////	virtual double AbstractState::d2alpha0_dTau2(void);
////	virtual double AbstractState::d3alpha0_dDelta3(void);
////	virtual double AbstractState::d3alpha0_dDelta2_dTau(void);
////	virtual double AbstractState::d3alpha0_dDelta_dTau2(void);
////	virtual double AbstractState::d3alpha0_dTau3(void);
////
////	virtual double AbstractState::alphar(void);
////	virtual double AbstractState::dalphar_dDelta(void);
////	virtual double AbstractState::dalphar_dTau(void);
////	virtual double AbstractState::d2alphar_dDelta2(void);
////	virtual double AbstractState::d2alphar_dDelta_dTau(void);
////	virtual double AbstractState::d2alphar_dTau2(void);
////	virtual double AbstractState::d3alphar_dDelta3(void);
////	virtual double AbstractState::d3alphar_dDelta2_dTau(void);
////	virtual double AbstractState::d3alphar_dDelta_dTau2(void);
////	virtual double AbstractState::d3alphar_dTau3(void);
////
////	// TODO: Add call back to calculator;
////	virtual double AbstractState::dalphar_dDelta_lim(void);
////	virtual double AbstractState::d2alphar_dDelta2_lim(void);
////	virtual double AbstractState::d2alphar_dDelta_dTau_lim(void);
////	virtual double AbstractState::d3alphar_dDelta2_dTau_lim(void);

} /* namespace CoolProp */
