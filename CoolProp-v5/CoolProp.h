
#ifndef CoolProp_H
#define CoolProp_H

#include <iostream>
#include <map>

class HashableObject {
	static int getHash();
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

