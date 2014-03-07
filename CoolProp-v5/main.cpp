
#include "Backends/REFPROPMixtureBackend.h"
#include "AbstractState.h"
#include "DataStructures.h"
using namespace CoolProp;

int main()
{
	std::vector<std::string> fluids;
	fluids.push_back("Methane");
	fluids.push_back("Ethane");
	AbstractState *State = new REFPROPMixtureBackend(fluids);
	std::vector<double> x(2,0.5);
	State->set_mole_fractions(x);
	State->update(DmassT_INPUTS,1,500);
	double hh = State->hmolar();
	double mu = State->viscosity();
	double sigma = State->surface_tension();
	delete State;
	
}