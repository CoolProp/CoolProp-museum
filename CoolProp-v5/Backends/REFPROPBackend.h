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

class REFPROPBackend : public REFPROPMixtureBackend  {
public:
	REFPROPBackend();
	REFPROPBackend(const std::string &fluid_name);
	
	// Setting mole or mass fraction not allowed since pure and pseudo-pure have single component with mole fraction of 1
	void set_mole_fractions(const std::vector<double> &mole_fractions){throw NotImplementedError();};
	void set_mass_fractions(const std::vector<double> &mass_fractions){throw NotImplementedError();};

	virtual ~REFPROPBackend();
};

} /* namespace CoolProp */
#endif /* REFPROPBACKEND_H_ */
