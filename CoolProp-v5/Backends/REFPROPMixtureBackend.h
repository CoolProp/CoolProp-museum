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
	/// @param input_pair Integer key corresponding to the two inputs that will be passed to the function
	/// @param value1 First input value
	/// @param value2 Second input value
	void update(int input_pair,
				double value1,
				double value2
				);

	void set_REFPROP_fluid(std::vector<std::string> fluid_names, std::vector<double> &mole_fractions);
};

} /* namespace CoolProp */
#endif /* REFPROPMIXTUREBACKEND_H_ */
