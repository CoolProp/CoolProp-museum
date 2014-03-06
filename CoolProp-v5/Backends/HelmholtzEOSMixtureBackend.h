/*
 * AbstractBackend.h
 *
 *  Created on: 20 Dec 2013
 *      Author: jowr
 */

#ifndef HELMHOLTZEOSBACKEND_H_
#define HELMHOLTZEOSBACKEND_H_

#include "../AbstractState.h"
#include <vector>

namespace CoolProp {

class HelmholtzEOSMixtureBackend : public AbstractState  {
public:
	HelmholtzEOSMixtureBackend();
	HelmholtzEOSMixtureBackend(std::vector<std::string> component_names, std::vector<double> mole_fractions);
	virtual ~HelmholtzEOSMixtureBackend();

	/// Updating function for pure and pseudo-pure fluids for REFPROP
	/// @param name1 First input index for state variable
	/// @param value1 First input value
	/// @param name2 Second input index for state variable
	/// @param value2 Second input value
	void update(int name1,
				double value1,
				int name2,
				double value2
				);

	void set_REFPROP_fluid(std::vector<std::string> fluid_names, std::vector<double> &mole_fractions);
};

} /* namespace CoolProp */
#endif /* HELMHOLTZEOSBACKEND_H_ */
