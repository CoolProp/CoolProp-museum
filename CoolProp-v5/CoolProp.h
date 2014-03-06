
#ifndef CoolProp_H
#define CoolProp_H

#include <iostream>
#include <map>
#include <string>
#include <vector>

class HashableObject {

	static long getHash(const std::vector<std::string> &components, const std::vector<double> &massFractions);

protected:
	long getHash(const std::string &hashstring);
	std::string getHashString(const std::vector<std::string> &components, const std::vector<double> &massFractions);

public:
	~HashableObject();
};

class StateObject : public HashableObject{
public:
	StateObject();
	~StateObject();
};

class FluidObject : public HashableObject
{
public:
	FluidObject();
	~FluidObject();
};


class CoolPropClass {

  public:
	static std::map<int, StateObject*> stateList;
	static std::map<int, FluidObject*> fluidList;

	/** Get a handle to a thermodynamic state for a pure or pseudo-pure fluid
	*/
	int request_handle(int fluid_index);

	/** Get a handle to a thermodynamic state for a mixture
	*/
	int request_handleM(std::vector<int> fluids, 
				        int molar_or_mass, 
				        std::vector<double> composition);

	long double IPropsSIM(int output,
						  int param1, 
			              double value1, 
						  int param2, 
						  double value2, 
						  std::vector<int> fluids, 
						  int molar_or_mass, 
						  std::vector<double> composition
						  );

	/// Property function that takes a handle to the state
	long double HIPropsSI(int output,
						  int param1, 
						  double value1, 
						  int param2, 
						  double value2,
						  int handle);

	/// The Props function using integers as the passed variables for speed
	long double IPropsSI(int output,
					     int param1, 
						 double value1, 
					     int param2, 
			             double value2, 
			             int fluid_index);

};

#endif /* CoolProp_H */

