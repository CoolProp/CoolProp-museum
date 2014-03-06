/*
 * AbstractBackend.h
 *
 *  Created on: 20 Dec 2013
 *      Author: jowr
 */

#ifndef HELMHOLTZEOSMIXTUREBACKEND_H_
#define HELMHOLTZEOSMIXTUREBACKEND_H_

#include "HelmholtzEOSMixtureBackend.h"
#include <vector>

namespace CoolProp {

class HelmholtzEOSBackend : public HelmholtzEOSMixtureBackend  {
public:
	HelmholtzEOSBackend();
	HelmholtzEOSBackend(std::string fluid_name);
	virtual ~HelmholtzEOSBackend();

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
#endif /* HELMHOLTZEOSMIXTUREBACKEND_H_ */