
#include "Backends/REFPROPMixtureBackend.h"
#include "AbstractState.h"
#include "DataStructures.h"
int main()
{
	std::vector<std::string> fluids;
	fluids.push_back("Methane");
	fluids.push_back("Ethane");
	CoolProp::AbstractState *State = new CoolProp::REFPROPMixtureBackend(fluids);
	std::vector<double> x(2,0.5);
	State->set_mole_fractions(x);
	State->update(CoolProp::DmassT_INPUTS,1,500);
	double hh = State->hmolar();
	delete State;
	
}