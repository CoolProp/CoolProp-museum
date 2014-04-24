/*
 * AbstractBackend.h
 *
 *  Created on: 20 Dec 2013
 *      Author: jowr
 */

#ifndef REFPROPBACKEND_H_
#define REFPROPBACKEND_H_

#include "REFPROPMixtureBackend.h"

namespace CoolProp {

/**
This backend is used for pure and pseudo-pure fluids powered by 
REFPROP.  It hides all the implementation of mixture properties
and exposes just the pure fluid interface.
*/
class REFPROPBackend : public REFPROPMixtureBackend  {
private:
	/// DISABLED: Setting mole or mass fraction not allowed since pure and pseudo-pure have single component with mole fraction of 1
	void set_mole_fractions(const std::vector<double> &mole_fractions){throw NotImplementedError();};
	/// DISABLED: Setting mole or mass fraction not allowed since pure and pseudo-pure have single component with mole fraction of 1
	void set_mass_fractions(const std::vector<double> &mass_fractions){throw NotImplementedError();};
public:
	REFPROPBackend();
	REFPROPBackend(const std::string &fluid_name);
	
	virtual ~REFPROPBackend();
};

} /* namespace CoolProp */
#endif /* REFPROPBACKEND_H_ */
