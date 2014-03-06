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
	virtual ~REFPROPBackend();

	/// Updating function for pure and pseudo-pure fluids for REFPROP
	/// @param name1 First input index for state variable
	/// @param value1 First input value
	/// @param name2 Second input index for state variable
	/// @param value2 Second input value
	/// @param State The state class to be updated by this call
	void update(int name1,
				double value1,
				int name2,
				double value2
				);

	void set_REFPROP_fluid(std::string fluid_name);
};

} /* namespace CoolProp */
#endif /* ABSTRACTBACKEND_H_ */
