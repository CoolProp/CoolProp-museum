/*
 * AbstractBackend.cpp
 *
 *  Created on: 20 Dec 2013
 *      Author: jowr
 */

#if defined(_MSC_VER)
#define _CRTDBG_MAP_ALLOC
#define _CRT_SECURE_NO_WARNINGS
#include <crtdbg.h>
#include <sys/stat.h>
#else
#include <sys/stat.h>
#endif

#include <string>
//#include "CoolProp.h"
#include "../CoolPropTools.h"

#include "REFPROPBackend.h"

namespace CoolProp {

REFPROPBackend::REFPROPBackend(const std::string & fluid_name) {
	// Do the REFPROP instantiation for this fluid

	// Try to add this fluid to REFPROP - might want to think about making array of 
	// components and setting mole fractions if they change a lot.
	std::vector<std::string> component_names(1,fluid_name);
	REFPROPMixtureBackend::set_REFPROP_fluids(component_names);
}

REFPROPBackend::~REFPROPBackend() {
	// TODO Auto-generated destructor stub
}

} /* namespace CoolProp */
