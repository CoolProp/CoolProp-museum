/*
 * AbstractBackend.h
 *
 *  Created on: 20 Dec 2013
 *      Author: jowr
 */

#ifndef HELMHOLTZEOSBACKEND_H_
#define HELMHOLTZEOSBACKEND_H_

#include "HelmholtzEOSMixtureBackend.h"
#include <vector>

namespace CoolProp {

class HelmholtzEOSBackend : public HelmholtzEOSMixtureBackend  {
public:
	HelmholtzEOSBackend();
	HelmholtzEOSBackend(std::string fluid_name);
	virtual ~HelmholtzEOSBackend();
	void set_REFPROP_fluid(std::vector<std::string> fluid_names, std::vector<double> &mole_fractions);
};

} /* namespace CoolProp */
#endif /* HELMHOLTZEOSBACKEND_H_ */
