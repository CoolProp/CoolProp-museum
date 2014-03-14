/*
 * AbstractBackend.h
 *
 *  Created on: 20 Dec 2013
 *      Author: jowr
 */

#ifndef HELMHOLTZEOSMIXTUREBACKEND_H_
#define HELMHOLTZEOSMIXTUREBACKEND_H_

#include "../AbstractState.h"
#include "../Fluids/CoolPropFluid.h"
#include <vector>

namespace CoolProp {

class HelmholtzEOSMixtureBackend : public AbstractState  {
protected:
	std::vector<CoolPropFluid*> components;
	bool is_pure_or_pseudopure;
public:
	HelmholtzEOSMixtureBackend(){};
	HelmholtzEOSMixtureBackend(std::vector<CoolPropFluid*> components);
	virtual ~HelmholtzEOSMixtureBackend(){};

	/// Updating function for pure and pseudo-pure fluids
	/// @param input_pair Integer key corresponding to the two inputs that will be passed to the function
	/// @param value1 First input value
	/// @param value2 Second input value
	void update(long input_pair,
				double value1,
				double value2
				){throw std::exception();};

	/// Set the components of the mixture
	/**
	@param components The components that are to be used in this mixture
	*/
	void set_components(std::vector<CoolPropFluid*> components);

	/// Set the mole fractions
	/** 
	@param mole_fractions The vector of mole fractions of the components
	*/
	void set_mole_fractions(const std::vector<double> &mole_fractions){throw std::exception();};
	
	/// Set the mass fractions
	/** 
	@param mass_fractions The vector of mass fractions of the components
	*/
	void set_mass_fractions(const std::vector<double> &mass_fractions){throw std::exception();};
};

} /* namespace CoolProp */
#endif /* HELMHOLTZEOSMIXTUREBACKEND_H_ */
