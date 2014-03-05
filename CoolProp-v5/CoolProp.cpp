
#include "CoolProp.h"

static long HashableObject::getHash(std::string hash_string){
	int summer = 0;
	for (unsigned int i = 0; i < hash_string.length(); i++)
	{
		summer += (int)hash_string[i];
	}
	return summer;
}

StateObject::StateObject() {
	// First we need to fill names and concentrations to
	// provide input to the hashing function!
	CoolPropClass::stateList.insert(std::make_pair<int, StateObject* >(this->getHash(),this));
}

StateObject::~StateObject(){
	CoolPropClass::stateList.erase(this->getHash());
}

FluidObject::FluidObject() {
	// First we need to fill names and concentrations to
	// provide input to the hashing function!
	CoolPropClass::fluidList.insert(std::make_pair<int, FluidObject* >(this->getHash(),this));
}

FluidObject::~FluidObject(){
	CoolPropClass::fluidList.erase(this->getHash());
}

class CoolPropClass {

  public:
	static std::map<int, StateObject*> stateList;
	static std::map<int, FluidObject*> fluidList;

};
