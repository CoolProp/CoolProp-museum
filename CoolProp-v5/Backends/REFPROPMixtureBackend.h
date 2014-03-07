/*
 * AbstractBackend.h
 *
 *  Created on: 20 Dec 2013
 *      Author: jowr
 */

#ifndef REFPROPMIXTUREBACKEND_H_
#define REFPROPMIXTUREBACKEND_H_

#include "../AbstractState.h"
#include <vector>

namespace CoolProp {

class REFPROPMixtureBackend : public AbstractState  {
protected:
	bool _mole_fractions_set;
	static bool _REFPROP_supported;
	std::vector<double> mole_fractions, mass_fractions;
	std::vector<double> mole_fractions_liq, mole_fractions_vap;
public:
	REFPROPMixtureBackend(){};
	REFPROPMixtureBackend(const std::vector<std::string>& fluid_names);
	virtual ~REFPROPMixtureBackend();

	/// Updating function for pure and pseudo-pure fluids for REFPROP
	/// @param input_pair Integer key corresponding to the two inputs that will be passed to the function
	/// @param value1 First input value
	/// @param value2 Second input value
	void update(long input_pair,
				double value1,
				double value2
				);

	/// Returns true if REFPROP is supported on this platform
	bool REFPROP_supported(void);

	/// Set the fluids in REFPROP DLL
	/// @param fluid_names The vector of strings of the fluid components, without file ending
	void set_REFPROP_fluids(const std::vector<std::string> &fluid_names);

	/// Set the mole fractions
	/// @param mole_fractions The vector of mole fractions of the components
	void set_mole_fractions(const std::vector<double> &mole_fractions);
	
	/// Set the mass fractions
	/// @param mass_fractions The vector of mass fractions of the components
	void set_mass_fractions(const std::vector<double> &mass_fractions);

	/// Check if the mole fractions have been set, etc.
	void check_status();

	/// Get the viscosity [Pa-s] for given temperature and density
	double calc_viscosity(void);
	/// Get the thermal conductivity [W/m/K] for given temperature and density
	double calc_conductivity(void);
	/// Get the surface tension [N/m] for given temperature
	double calc_surface_tension(void);

};

} /* namespace CoolProp */
#endif /* REFPROPMIXTUREBACKEND_H_ */
