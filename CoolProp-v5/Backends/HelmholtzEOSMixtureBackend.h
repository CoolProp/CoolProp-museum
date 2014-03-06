/*
 * AbstractBackend.h
 *
 *  Created on: 20 Dec 2013
 *      Author: jowr
 */

#ifndef HELMHOLTZEOSMIXTUREBACKEND_H_
#define HELMHOLTZEOSMIXTUREBACKEND_H_

#include "../AbstractState.h"
#include <vector>

namespace CoolProp {

class HelmholtzEOSMixtureBackend : public AbstractState  {
public:
	HelmholtzEOSMixtureBackend();
	HelmholtzEOSMixtureBackend(std::vector<std::string> component_names, std::vector<double> mole_fractions);
	virtual ~HelmholtzEOSMixtureBackend();

	/// Updating function for pure and pseudo-pure fluids
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
#endif /* HELMHOLTZEOSMIXTUREBACKEND_H_ */
