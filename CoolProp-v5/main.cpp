
#include "Backends/REFPROPMixtureBackend.h"
#include "Backends/REFPROPBackend.h"
#include <time.h>
#include "AbstractState.h"
#include "DataStructures.h"
#include <cstdio>
using namespace CoolProp;

#include "jsoncons/json.hpp"

class JSONFluidParser
{
protected:

	/// Parse the contributions to the residual Helmholtz energy
	void parse_alphar(jsoncons::json& alphar)
	{
		for (size_t i = 0; i < alphar.size(); ++i)
		{
			// Reference to this contribution entry
			jsoncons::json& contribution = alphar[i];
			
			// Get the type (required!)
			std::string type = contribution["type"].as<std::string>();

			if (!type.compare("alphar_power"))
			{
				std::vector<double> n = contribution["n"].as<std::vector<double> >();
				std::vector<double> d = contribution["d"].as<std::vector<double> >();
				std::vector<double> t = contribution["t"].as<std::vector<double> >();
				std::vector<double> l = contribution["l"].as<std::vector<double> >();
			}
			else
			{
				throw ValueError("Unsupported alphar type");
			}
		}
	};

	/// Parse the contributions to the ideal-gas Helmholtz energy
	void parse_alpha0(const jsoncons::json& alpha0)
	{
	};

	/// Parse the Equation of state JSON
	void parse_EOS(jsoncons::json& EOS)
	{
		parse_alphar(EOS["alphar"]);
	}

	/// Parse the list of possible equations of state
	void parse_EOS_listing(jsoncons::json& EOS_array)
	{
		for (size_t i = 0; i < EOS_array.size(); ++i)
		{
			parse_EOS(EOS_array[i]);
		}
	};

	/// Parse the reducing state for the given EOS
	void parse_reducing_state(jsoncons::json& alphar)
	{
	};

	/// Parse the critical state for the given EOS
	void parse_crit_state(jsoncons::json& alphar)
	{
	};

public:
	JSONFluidParser(jsoncons::json& fluid)
	{
		std::string name                 = fluid["name"].as<std::string>();
		std::vector<std::string> aliases = fluid["aliases"].as<std::vector<std::string> >();
		
		parse_EOS_listing(fluid["EOS"]);
		double rr;
	}
};

int main()
{
	if(0)
	{
		jsoncons::json fluids = jsoncons::json::parse_file("fluids.json");

		for (size_t i = 0; i < fluids.size(); ++i)
		{
			try
			{
				JSONFluidParser JFP(fluids[i]);
				double rr = 0;
			}
			catch (const std::exception& e)
			{
				std::cerr << e.what() << std::endl;
			}
		}
        double rr = 0;
	}
	if (0)
	{
		AbstractState *State = AbstractState::factory("REFPROP-Methane|Ethane");
		std::vector<double> x(2,0.5);
		State->set_mole_fractions(x);
		State->update(DmassT_INPUTS,1,250);
		double hh = State->hmolar();
		double mu = State->viscosity();
		double sigma = State->surface_tension();
		delete State;
	}
	if (1)
	{
		time_t t1,t2;
		t1 = clock();
		long N = 100000;
		for (long ii = 0; ii < N; ii++)
		{
			AbstractState *State = AbstractState::factory("REFPROP-Methane");
			//AbstractState *State = new REFPROPBackend("Methane");
			delete State;
		}
		t2 = clock();
		double elap = ((double)(t2-t1))/CLOCKS_PER_SEC/((double)N)*1e6;
		printf("%g\n",elap);
	}

	if(0)
	{
		AbstractState *State = AbstractState::factory("REFPROP-Methane");
		State->update(DmassT_INPUTS,1,300);
		double hh = State->hmolar();
		double mu = State->viscosity();

		time_t t1,t2;
		t1 = clock();
		for (long ii = 0; ii < 1000000; ii++)
		{
			State->update(PQ_INPUTS,300000,1-ii*1e-6);
			//State->update(DmassT_INPUTS,1-ii*1e-10,180);
			//double hh1 = State->hmolar();
			//double mu2 = State->viscosity();
		}
		t2 = clock();
		double elap = ((double)(t2-t1))/CLOCKS_PER_SEC;
		printf("%g\n",elap);

		//double sigma = State->surface_tension();
		delete State;
	}

	
	
	
}