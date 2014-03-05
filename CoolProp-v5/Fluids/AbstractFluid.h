/*
 * AbtractFluid.h
 *
 *  Created on: 20 Dec 2013
 *      Author: jowr
 */

#ifndef ABSTRACTFLUID_H_
#define ABSTRACTFLUID_H_

#include "AbstractBackend.h"

namespace CoolProp {

class AbstractFluid: public CoolProp::AbstractBackend {
public:
	AbstractFluid();
	virtual ~AbstractFluid();
};

} /* namespace CoolProp */
#endif /* ABSTRACTFLUID_H_ */
