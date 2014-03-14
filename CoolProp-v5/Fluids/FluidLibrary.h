
#ifndef FLUIDLIBRARY_H
#define FLUIDLIBRARY_H

#include <map>

#include "../rapidjson/rapidjson_include.h"

#include "CoolPropFluid.h"

namespace CoolProp{
/// A container for the fluid parameters for the CoolProp fluids
/**
This container holds copies of all of the fluid instances for the fluids that are loaded in CoolProp. 
New fluids can be added by passing in a rapidjson::Value instance to the add_one function, or 
a rapidjson array of fluids to the add_many function.
*/
class JSONFluidLibrary
{
	/// Map from CAS code to JSON instance.  For pseudo-pure fluids, use name in place of CAS code since no CASE number is defined for mixtures
	std::map<long, CoolPropFluid> fluid_map;

	std::map<std::string, long> string_to_index_map;
	bool _is_empty;
protected:

	/// Parse the contributions to the residual Helmholtz energy
	void parse_alphar(rapidjson::Value &alphar, EquationOfState &EOS)
	{	
		for (rapidjson::Value::ValueIterator itr = alphar.Begin(); itr != alphar.End(); ++itr)
		{	
			// A reference for code cleanness
			rapidjson::Value &contribution = *itr;

			// Get the type (required!)
			std::string type = contribution["type"].GetString();

			if (!type.compare("alphar_power"))
			{
				std::vector<double> n = cpjson::get_double_array(contribution["n"]);
				std::vector<double> d = cpjson::get_double_array(contribution["d"]);
				std::vector<double> t = cpjson::get_double_array(contribution["t"]);
				std::vector<double> l = cpjson::get_double_array(contribution["l"]);
				EOS.alphar.push_back(new alphar_power());
			}
			else if (!type.compare("alphar_gaussian"))
			{
				std::vector<double> n = cpjson::get_double_array(contribution["n"]);
				std::vector<double> d = cpjson::get_double_array(contribution["d"]);
				std::vector<double> t = cpjson::get_double_array(contribution["t"]);
				std::vector<double> eta = cpjson::get_double_array(contribution["eta"]);
				std::vector<double> epsilon = cpjson::get_double_array(contribution["epsilon"]);
				std::vector<double> beta = cpjson::get_double_array(contribution["beta"]);
				std::vector<double> gamma = cpjson::get_double_array(contribution["gamma"]);
				EOS.alphar.push_back(new alphar_power());
			}
			else if (!type.compare("alphar_critical"))
			{
				std::vector<double> n = cpjson::get_double_array(contribution["n"]);
				std::vector<double> a = cpjson::get_double_array(contribution["a"]);
				std::vector<double> b = cpjson::get_double_array(contribution["b"]);
				std::vector<double> beta = cpjson::get_double_array(contribution["beta"]);
				std::vector<double> A = cpjson::get_double_array(contribution["A"]);
				std::vector<double> B = cpjson::get_double_array(contribution["B"]);
				std::vector<double> C = cpjson::get_double_array(contribution["C"]);
				std::vector<double> D = cpjson::get_double_array(contribution["D"]);
				EOS.alphar.push_back(new alphar_power());
			}
			else
			{
				throw ValueError("Unsupported alphar type");
			}
		}
	};

	/// Parse the contributions to the ideal-gas Helmholtz energy
	void parse_alpha0(rapidjson::Value &alpha0, EquationOfState &EOS)
	{
		EOS.alpha0.push_back(new alphar_power());
	};

	/// Parse the Equation of state JSON entry
	void parse_EOS(rapidjson::Value &EOS_json, CoolPropFluid &fluid)
	{
		EquationOfState EOS;

		// Universal gas constant [J/mol/K]
		EOS.R_u = cpjson::get_double(EOS_json,"R_u");
		
		parse_alphar(EOS_json["alphar"], EOS);
		parse_alpha0(EOS_json["alpha0"], EOS);
		
		// Validate the equation of state that was just created
		EOS.validate();

		// TODO: Get reducing state
		fluid.EOSVector.push_back(EOS);
	}

	/// Parse the list of possible equations of state
	void parse_EOS_listing(rapidjson::Value &EOS_array, CoolPropFluid & fluid)
	{
		for (rapidjson::Value::ValueIterator itr = EOS_array.Begin(); itr != EOS_array.End(); ++itr)
		{	
			parse_EOS(*itr,fluid);
		}
	};

	/// Parse the reducing state for the given EOS
	void parse_reducing_state(rapidjson::Value &alphar)
	{
	};

	/// Parse the critical state for the given EOS
	void parse_crit_state(rapidjson::Value &alphar)
	{
	};

	/// Validate the fluid file that was just constructed
	void validate(CoolPropFluid & fluid)
	{
		assert(fluid.EOSVector.size() > 0);
		assert(fluid.CAS.length() > 0);
		assert(fluid.name.length() > 0);
	}
public:
	// Default constructor;
	JSONFluidLibrary(){
		_is_empty = true;
	};
	bool is_empty(void){ return _is_empty;};

	/// Add all the fluid entries in the rapidjson::Value instance passed in
	void add_many(rapidjson::Value &listing)
	{
		for (rapidjson::Value::ValueIterator itr = listing.Begin(); itr != listing.End(); ++itr)
		{	
			add_one(*itr);
		}
	};
	void add_one(rapidjson::Value &fluid_json)
	{
		_is_empty = false;
		// Create an instance of the fluid
		CoolPropFluid fluid;

		// Fluid name
		fluid.name = fluid_json["NAME"].GetString();
		// CAS number
		fluid.CAS = fluid_json["CAS"].GetString();

		// Aliases
		fluid.aliases = cpjson::get_string_array(fluid_json["ALIASES"]);
		
		// EOS
		parse_EOS_listing(fluid_json["EOS"],fluid);

		// Validate the fluid
		validate(fluid);
		
		// If the fluid is ok...

		// Get the next index for this fluid
		long index = fluid_map.size();

		// Add index->fluid mapping
		fluid_map[index] = fluid;

		// Add CAS->index mapping
		string_to_index_map[fluid.CAS] = index;

		// Add name->index mapping
		string_to_index_map[fluid.name] = index;

	};
	/// Get a CoolPropFluid instance stored in this library
	/**
	@param key Either a CAS number or the name (CAS number should be preferred)
	*/
	CoolPropFluid& get(std::string key)
	{
		std::map<std::string, long>::iterator it;
		// Try to find it
		it = string_to_index_map.find(key);
		// If it is found
		if (it != string_to_index_map.end()){
			return get(it->second);
		}
		else{
			throw ValueError(format("key [%s] was not found in string_to_index_map in JSONFluidLibrary",key.c_str()));
		}
	};
	/// Get a CoolPropFluid instance stored in this library
	/**
	@param key The index of the fluid in the map
	*/
	CoolPropFluid& get(long key)
	{
		std::map<long,CoolPropFluid>::iterator it;
		// Try to find it
		it = fluid_map.find(key);
		// If it is found
		if (it != fluid_map.end()){
			return it->second;
		}
		else{
			throw ValueError(format("key [%d] was not found in JSONFluidLibrary",key));
		}
	};
};

/// Get a reference to the library instance
JSONFluidLibrary & get_library(void);

} /* namespace CoolProp */
#endif