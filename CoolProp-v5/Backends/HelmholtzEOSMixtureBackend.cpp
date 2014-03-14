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

#include "HelmholtzEOSMixtureBackend.h"

namespace CoolProp {

HelmholtzEOSMixtureBackend::HelmholtzEOSMixtureBackend(std::vector<CoolPropFluid*> components) {
	
	/// Set the components and associated flags
	set_components(components);
}

void HelmholtzEOSMixtureBackend::set_components(std::vector<CoolPropFluid*> components) {

	/// Copy the components
	this->components = components;

	if (components.size() == 1){
		is_pure_or_pseudopure = true;
	}
	else{
		is_pure_or_pseudopure = false;
	}
}

} /* namespace CoolProp */
