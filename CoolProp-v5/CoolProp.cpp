
#include "CoolProp.h"

static int HashableObject::getHash(){
	return 0;
}

StateObject::StateObject
{
	public:
	StateObject() {
	    	// First we need to fill names and concentrations to
	    	// provide input to the hashing function!
	    	CoolPropClass::stateList.insert(std::make_pair<int, StateObject* >(this->getHash(),this));
	    }

		~StateObject(){
			CoolPropClass::stateList.erase(this->getHash());
		}
	};

class FluidObject : public HashableObject
	{
	public:
		FluidObject() {
	    	// First we need to fill names and concentrations to
	    	// provide input to the hashing function!
	    	CoolPropClass::fluidList.insert(std::make_pair<int, FluidObject* >(this->getHash(),this));
		}

		~FluidObject(){
			CoolPropClass::fluidList.erase(this->getHash());
		}
	};


class CoolPropClass {

  public:
	static std::map<int, StateObject* > stateList;
	static std::map<int, FluidObject* > fluidList;



};
