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

HelmholtzEOSMixtureBackend::HelmholtzEOSMixtureBackend(std::vector<std::string> component_names, std::vector<double> mole_fractions) {
	// TODO Auto-generated constructor stub
}

HelmholtzEOSMixtureBackend::~HelmholtzEOSMixtureBackend() {
	// TODO Auto-generated destructor stub
}

} /* namespace CoolProp */
