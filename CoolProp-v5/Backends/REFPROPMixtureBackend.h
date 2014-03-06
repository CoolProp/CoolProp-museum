/*
 * AbstractBackend.h
 *
 *  Created on: 20 Dec 2013
 *      Author: jowr
 */

#ifndef REFPROPMIXTUREBACKEND_H_
#define REFPROPMIXTUREBACKEND_H_

#include "../AbstractState.h"
#include <vector>

namespace CoolProp {

class REFPROPMixtureBackend : public AbstractState  {
public:
	REFPROPMixtureBackend();
	virtual ~REFPROPMixtureBackend();

	/// Updating function for pure and pseudo-pure fluids for REFPROP
	/// @param name1 First input index for state variable
	/// @param value1 First input value
	/// @param value2 Second input value
	void update(int name1,
				double value1,
				double value2
				);

	void set_REFPROP_fluid(std::vector<std::string> fluid_names, std::vector<double> &mole_fractions);
};

} /* namespace CoolProp */
#endif /* REFPROPMIXTUREBACKEND_H_ */
