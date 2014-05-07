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
public:
    
	REFPROPBackend();
	REFPROPBackend(const std::string &fluid_name);
	
	virtual ~REFPROPBackend();
};

} /* namespace CoolProp */
#endif /* REFPROPBACKEND_H_ */
