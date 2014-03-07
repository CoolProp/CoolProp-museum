
#include "Backends/REFPROPMixtureBackend.h"
#include "Backends/REFPROPBackend.h"

#include "AbstractState.h"
#include "DataStructures.h"
using namespace CoolProp;

int main()
{
	if (1)
	{
		std::vector<std::string> fluids;
		fluids.push_back("Methane");
		fluids.push_back("Ethane");

		AbstractState *State = new REFPROPMixtureBackend(fluids);
		std::vector<double> x(2,0.5);
		State->set_mole_fractions(x);
		State->update(DmassT_INPUTS,1,250);
		double hh = State->hmolar();
		double mu = State->viscosity();
		double sigma = State->surface_tension();
		delete State;
	}

	if(1)
	{
		AbstractState *State = new REFPROPBackend("Methane");
		State->update(DmassT_INPUTS,1,300);
		double hh = State->hmolar();
		double mu = State->viscosity();

		State->update(QT_INPUTS,1,180);
		double hh1 = State->hmolar();
		double mu2 = State->viscosity();

		//double sigma = State->surface_tension();
		delete State;
	}

	
	
	
}