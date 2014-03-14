
#include "Backends/REFPROPMixtureBackend.h"
#include "Backends/REFPROPBackend.h"
#include <time.h>
#include "AbstractState.h"
#include "DataStructures.h"
#include <deque>
#include <cstdio>
using namespace CoolProp;


#include "rapidjson/rapidjson_include.h"
#include "Fluids\FluidLibrary.h"

int main()
{
	if(1)
	{
		JSONFluidLibrary JFL;
		
		{
		rapidjson::Document dd;
		dd.Parse<0>(get_file_contents("../../../CoolProp/Water.json").c_str());
		if (dd.HasParseError()){throw ValueError();} else{JFL.add_one(dd);}		
		}

		time_t t1,t2;
		t1 = clock();
		long N = 10000;
		for (long ii = 0; ii < N; ii++)
		{
			AbstractState *State = AbstractState::factory("CORE-Water");
			//AbstractState *State = new REFPROPBackend("Methane");
			delete State;
		}
		t2 = clock();
		double elap = ((double)(t2-t1))/CLOCKS_PER_SEC/((double)N)*1e6;
		printf("%g\n",elap);
		double eee = 0;

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