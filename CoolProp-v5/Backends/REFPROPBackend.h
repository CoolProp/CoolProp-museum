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
	REFPROPBackend(std::string fluid_name);
	virtual ~REFPROPBackend();
};

} /* namespace CoolProp */
#endif /* REFPROPBACKEND_H_ */
