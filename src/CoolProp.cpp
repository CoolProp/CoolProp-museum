
#include "CoolProp.h"
#include "CoolPropTools.h"
#include <assert.h>     /* assert */

//long HashableObject::getHash(const std::vector<std::string>& components, const std::vector<double>& massFractions)
//{
//	assert(components.size() == massFractions.size());
//	// Join the components by a '|'
//	std::string components_joined = strjoin(components,"|");
//	std::string massFractions_joined;
//	for (unsigned int i = 0; i < massFractions_joined.size(); i++)
//	{
//		massFractions_joined += format("%0.12e",massFractions_joined[i]);
//	}
//	return getHash(components_joined + massFractions_joined);
//}
//
//static long HashableObject::getHash(std::string hash_string){
//	int summer = 0;
//	for (unsigned int i = 0; i < hash_string.length(); i++)
//	{
//		summer += (int)hash_string[i];
//	}
//	return summer;
//}
//
//StateObject::StateObject() {
//	// First we need to fill names and concentrations to
//	// provide input to the hashing function!
//	CoolPropClass::stateList.insert(std::make_pair<int, StateObject* >(this->getHash(),this));
//}
//
//StateObject::~StateObject(){
//	CoolPropClass::stateList.erase(this->getHash());
//}
//
//FluidObject::FluidObject() {
//	// First we need to fill names and concentrations to
//	// provide input to the hashing function!
//	CoolPropClass::fluidList.insert(std::make_pair<int, FluidObject* >(this->getHash(),this));
//}
//
//FluidObject::~FluidObject(){
//	CoolPropClass::fluidList.erase(this->getHash());
//}
