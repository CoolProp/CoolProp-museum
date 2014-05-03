
#include "FluidLibrary.h"
#include "../Backends/ReducingFunctions.h"

namespace CoolProp{

static JSONFluidLibrary library;

void load()
{
	rapidjson::Document dd;
    dd.Parse<0>(get_file_contents("all_fluids.json").c_str());
	if (dd.HasParseError()){throw ValueError("Unable to load all_fluids.json");} else{library.add_many(dd);}
}

JSONFluidLibrary & get_library(void){
	if (library.is_empty()){
		load();
	}
	return library;
}

} /* namespace CoolProp */