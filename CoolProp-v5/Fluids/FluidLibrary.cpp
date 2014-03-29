
#include "FluidLibrary.h"
#include "../Backends/ReducingFunctions.h"

namespace CoolProp{

static JSONFluidLibrary library;

void load()
{
	rapidjson::Document dd;
	//std::vector<BaseHelmholtzTerm*> ar;
    dd.Parse<0>(get_file_contents("../../../CoolProp/n-Propane.json").c_str());
	if (dd.HasParseError()){throw ValueError();} else{library.add_one(dd);}
    dd.Parse<0>(get_file_contents("../../../CoolProp/Water.json").c_str());
	if (dd.HasParseError()){throw ValueError();} else{library.add_one(dd);}
    dd.Parse<0>(get_file_contents("../../../CoolProp/Methane.json").c_str());
	if (dd.HasParseError()){throw ValueError();} else{library.add_one(dd);}
    dd.Parse<0>(get_file_contents("../../../CoolProp/Ethane.json").c_str());
	if (dd.HasParseError()){throw ValueError();} else{library.add_one(dd);}
}

JSONFluidLibrary & get_library(void){
	if (library.is_empty()){
		load();
	}
	return library;
}

} /* namespace CoolProp */