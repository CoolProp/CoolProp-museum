
#include "FluidLibrary.h"

namespace CoolProp{

static JSONFluidLibrary library;

void load()
{
	rapidjson::Document dd;
	dd.Parse<0>(get_file_contents("../../../CoolProp/Water.json").c_str());
	if (dd.HasParseError()){throw ValueError();} else{library.add_one(dd);}
}

JSONFluidLibrary & get_library(void){
	if (library.is_empty()){
		load();
	}
	return library;
}

} /* namespace CoolProp */