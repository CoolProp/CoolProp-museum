
#ifndef CoolProp_H
#define CoolProp_H

#include <iostream>
#include <map>
#include <string>
#include <vector>

class HashableObject {

	static long getHash(const std::string components, const std::vector<double>* massFractions);

protected:
	long getHash(std::string);
	std::string getHashString(const std::string components, const std::vector<double>* massFractions);

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
	static std::map<int, StateObject* > stateList;
	static std::map<int, FluidObject* > fluidList;



};

#endif /* CoolProp_H */

