#if defined(_MSC_VER)
#define _CRTDBG_MAP_ALLOC
#define _CRT_SECURE_NO_WARNINGS
#include <crtdbg.h>
#endif

#include "CoolProp.h"
#include "REFPROP.h"


#if defined(__ISWINDOWS__)
#include <windows.h>
#else
#ifndef DBL_EPSILON
	#include <limits>
	#define DBL_EPSILON std::numeric_limits<double>::epsilon()
#endif
#endif

#include <iostream>
#include <stdlib.h>
#include <vector>
#include <exception>
#include <stdio.h>
#include "string.h"
#include "FluidClass.h"
#include "CoolPropTools.h"
#include "CPExceptions.h"
#include "Brine.h"
#include "Solvers.h"
#include "CPState.h"
#include "IncompLiquid.h"
#include "TTSE.h"
#include "purefluids/R290.h"
#include "purefluids/R134a.h"

// Function prototypes
//double rho_TP(double T, double p);
double _Props(std::string Output,std::string Name1, double Prop1, std::string Name2, double Prop2, std::string Ref);
double _CoolProp_Fluid_Props(long iOutput, long iName1, double Value1, long iName2, double Value2, Fluid *pFluid, bool SinglePhase = false);

static std::string err_string;
static std::string warning_string;
static int debug_level=0;
static Fluid * pFluid;

// This is very hacky, but pull the git revision from the file
#include "gitrevision.h" // Contents are like "long gitrevision = "aa121435436ggregrea4t43t433";"
#include "version.h" // Contents are like "char version [] = "2.5";"

int global_Phase = -1;
bool global_SinglePhase = false;
bool global_SaturatedL = false;
bool global_SaturatedV = false;

// Default to the KSI unit system
int unit_system = UNIT_SYSTEM_KSI;

void set_warning(std::string warning){ warning_string = warning; }

int _get_standard_unit_system(){ return unit_system; }
void _set_standard_unit_system(int unit_sys){ unit_system = unit_sys; }

// This is a map of all possible strings to a unique identifier
std::pair<std::string, long> map_data[] = {
	std::make_pair(std::string("E"),iPcrit),
	std::make_pair(std::string("M"),iMM),
	std::make_pair(std::string("w"),iAccentric),
	std::make_pair(std::string("R"),iTtriple),
	std::make_pair(std::string("N"),iRhocrit),
	std::make_pair(std::string("B"),iTcrit),

	std::make_pair(std::string("pcrit"),iPcrit),
	std::make_pair(std::string("molemass"),iMM),
	std::make_pair(std::string("accentric"),iAccentric),
	std::make_pair(std::string("dipole"),iDipole),
	std::make_pair(std::string("Tmin"),iTmin),
	std::make_pair(std::string("t"),iTmin),
	std::make_pair(std::string("Ttriple"),iTtriple),
	std::make_pair(std::string("ptriple"),iPtriple),
	std::make_pair(std::string("rhocrit"),iRhocrit),
	std::make_pair(std::string("Tcrit"),iTcrit),
	std::make_pair(std::string("Treduce"),iTreduce),
	std::make_pair(std::string("rhoreduce"),iRhoreduce),
	std::make_pair(std::string("Hcrit"),iHcrit),
	std::make_pair(std::string("Scrit"),iScrit),

	std::make_pair(std::string("Q"),iQ),
	std::make_pair(std::string("T"),iT),
    std::make_pair(std::string("P"),iP),
	std::make_pair(std::string("D"),iD),
	std::make_pair(std::string("C"),iC),
	std::make_pair(std::string("C0"),iC0),
	std::make_pair(std::string("O"),iO),
	std::make_pair(std::string("U"),iU),
	std::make_pair(std::string("H"),iH),
	std::make_pair(std::string("S"),iS),
	std::make_pair(std::string("A"),iA),
	std::make_pair(std::string("G"),iG),
	std::make_pair(std::string("V"),iV),
	std::make_pair(std::string("L"),iL),
	std::make_pair(std::string("Tmax"),iTmax),
	std::make_pair(std::string("Tfreeze"),iTfreeze),
	std::make_pair(std::string("Psat"),iPsat),
	std::make_pair(std::string("I"),iI),
	std::make_pair(std::string("SurfaceTension"),iI),
	std::make_pair(std::string("Phase"),iPhase),
	std::make_pair(std::string("PHASE_LIQUID"),iPHASE_LIQUID),
	std::make_pair(std::string("PHASE_GAS"),iPHASE_GAS),
	std::make_pair(std::string("PHASE_SUPERCRITICAL"),iPHASE_SUPERCRITICAL),
	std::make_pair(std::string("PHASE_TWOPHASE"),iPHASE_TWOPHASE),
	std::make_pair(std::string("ODP"),iODP),
	std::make_pair(std::string("GWP20"),iGWP20),
	std::make_pair(std::string("GWP100"),iGWP100),
	std::make_pair(std::string("GWP500"),iGWP500),
	std::make_pair(std::string("CritSplineT"),iCritSplineT),
	// Derivatives
	std::make_pair(std::string("dhdp|rho"),iDERdh_dp__rho),
	std::make_pair(std::string("Z"),iDERZ),
	std::make_pair(std::string("dZ_dDelta"),iDERdZ_dDelta),
	std::make_pair(std::string("dZ_dTau"),iDERdZ_dTau),
	std::make_pair(std::string("VB"),iDERB),
	std::make_pair(std::string("dBdT"),iDERdB_dT),
	std::make_pair(std::string("VC"),iDERC),
	std::make_pair(std::string("dCdT"),iDERdC_dT),
	std::make_pair(std::string("phir"),iDERphir),
	std::make_pair(std::string("dphir_dTau"),iDERdphir_dTau),
	std::make_pair(std::string("dphir_dDelta"),iDERdphir_dDelta),
	std::make_pair(std::string("d2phir_dTau2"),iDERd2phir_dTau2),
	std::make_pair(std::string("d2phir_dDelta2"),iDERd2phir_dDelta2),
	std::make_pair(std::string("d2phir_dDelta_dTau"),iDERd2phir_dDelta_dTau),
	std::make_pair(std::string("d3phir_dDelta3"),iDERd3phir_dDelta3),
	std::make_pair(std::string("d3phir_dDelta2_dTau"),iDERd3phir_dDelta2_dTau),
	std::make_pair(std::string("d3phir_dDelta_dTau2"),iDERd3phir_dDelta_dTau2),
	std::make_pair(std::string("d3phir_dTau3"),iDERd3phir_dTau3),
	std::make_pair(std::string("phi0"),iDERphi0),
	std::make_pair(std::string("dphi0_dTau"),iDERdphi0_dTau),
	std::make_pair(std::string("d2phi0_dTau2"),iDERd2phi0_dTau2),
	std::make_pair(std::string("dphi0_dDelta"),iDERdphi0_dDelta),
	std::make_pair(std::string("d2phi0_dDelta2"),iDERd2phi0_dDelta2),
	std::make_pair(std::string("d2phi0_dDelta_dTau"),iDERd2phi0_dDelta_dTau),
	std::make_pair(std::string("d3phi0_dTau3"),iDERd3phi0_dTau3),
	std::make_pair(std::string("dpdT"    ),iDERdp_dT__rho),
	std::make_pair(std::string("dpdT|rho"),iDERdp_dT__rho),
	std::make_pair(std::string("dpdrho"  ),iDERdp_drho__T),
	std::make_pair(std::string("dpdrho|T"),iDERdp_drho__T),
	std::make_pair(std::string("dhdT|rho"),iDERdh_dT__rho),
	std::make_pair(std::string("dhdrho|T"),iDERdh_drho__T),
	std::make_pair(std::string("drhodT|p"),iDERdrho_dT__p),
	std::make_pair(std::string("drhodh|p"),iDERdrho_dh__p),
	std::make_pair(std::string("drhodp|h"),iDERdrho_dp__h),
	std::make_pair(std::string("rho_smoothed"),iDERrho_smoothed),
	std::make_pair(std::string("d(rho_smoothed)/dh"),iDERdrho_smoothed_dh),
	std::make_pair(std::string("d(rho_smoothed)/dp"),iDERdrho_smoothed_dp),
	std::make_pair(std::string("drhodh_constp_smoothed"),iDERdrhodh_constp_smoothed),
	std::make_pair(std::string("drhodp_consth_smoothed"),iDERdrhodp_consth_smoothed),
	std::make_pair(std::string("IsothermalCompressibility"),iDERIsothermalCompressibility)
};

//Now actually construct the map
std::map<std::string, long> param_map(map_data,
    map_data + sizeof map_data / sizeof map_data[0]);

// This is a map of all unique identifiers to std::strings with the default units
// that are used internally
std::pair<long, std::string> units_data[] = {
	std::make_pair(iPcrit, std::string("kPa")),
	std::make_pair(iMM, std::string("kg/kmol")),
	std::make_pair(iAccentric, std::string("-")),
	std::make_pair(iTtriple, std::string("K")),
	std::make_pair(iRhocrit, std::string("kg/m^3")),
	std::make_pair(iTcrit, std::string("K")),
	std::make_pair(iTmin, std::string("K")),

	std::make_pair(iQ, std::string("")),
	std::make_pair(iT, std::string("K")),
    std::make_pair(iP, std::string("kPa")),
	std::make_pair(iD, std::string("kg/m^3")),
	std::make_pair(iC, std::string("kJ/kg/K")),
	std::make_pair(iC0, std::string("kJ/kg/K")),
	std::make_pair(iO, std::string("kJ/kg/K")),
	std::make_pair(iU, std::string("kJ/kg")),
	std::make_pair(iH, std::string("kJ/kg")),
	std::make_pair(iS, std::string("kJ/kg/K")),
	std::make_pair(iA, std::string("m/s")),
	std::make_pair(iG, std::string("kJ/kg")),
	std::make_pair(iV, std::string("Pa*s")),
	std::make_pair(iL, std::string("kW/m/K")),
	std::make_pair(iTmax, std::string("K")),
	std::make_pair(iTfreeze, std::string("K")),
	std::make_pair(iPsat, std::string("kPa")),
	std::make_pair(iI, std::string("N/m")),
	std::make_pair(iDpdT, std::string("kPa/K")),
	std::make_pair(iDrhodT_p, std::string("kg/K/m^3")),
	std::make_pair(iDERdp_dT__rho, std::string("kPa/K")),
	std::make_pair(iDERdrho_dT__p, std::string("kg/K/m^3"))
};

//Now actually construct the map
std::map<long, std::string> units_map(units_data,
		units_data + sizeof units_data / sizeof units_data[0]);

FluidsContainer Fluids = FluidsContainer();

void set_err_string(std::string error_string)
{
	err_string = error_string;
}

void set_debug_level(int level){debug_level = level;}
int get_debug_level(){return debug_level;}

void set_phase(std::string Phase_str){
	if (!Phase_str.compare("Two-Phase")){
		global_SinglePhase = false;
		global_SaturatedL = false;
		global_SaturatedV = false;
		global_Phase = iTwoPhase;
	}
	else if (!Phase_str.compare("Liquid")){
		global_SinglePhase = true;
		global_SaturatedL = false;
		global_SaturatedV = false;
		global_Phase = iLiquid;
	}
	else if (!Phase_str.compare("Gas")){
		global_SinglePhase = true;
		global_SaturatedL = false;
		global_SaturatedV = false;
		global_Phase = iGas;
	}
	else if (!Phase_str.compare("Supercritical")){
		global_SinglePhase = true;
		global_SaturatedL = false;
		global_SaturatedV = false;
		global_Phase = iSupercritical;
	}
	else if (!Phase_str.compare("SaturatedL")){
		global_SinglePhase = false;
		global_SaturatedL = true;
		global_SaturatedV = false;
		global_Phase = iTwoPhase;
	}
	else if (!Phase_str.compare("SaturatedV")){
		global_SinglePhase = false;
		global_SaturatedL = false;
		global_SaturatedV = true;
		global_Phase = iTwoPhase;
	}
}

// Returns a pointer to the fluid class
Fluid* get_fluid(long iFluid){
	return Fluids.get_fluid(iFluid);
}
long get_Fluid_index(std::string FluidName)
{
	// Try to get the fluid from Fluids by name
	pFluid = Fluids.get_fluid(FluidName);
	// If NULL, didn't find it (or its alias)
	if (pFluid!=NULL)
	{
		// Find the fluid index
		return Fluids.get_fluid_index(pFluid);
	}
	else
		return -1;
}

bool add_REFPROP_fluid(std::string FluidName)
{
	double x[100];
	// Starts with REFPROP- keep going
	if (!(FluidName.find("REFPROP-") == 0)) return false;
	// Stop here if there is no REFPROP support
	if (!REFPROPFluidClass::refpropSupported()) return false;
	// Try to load this fluid, index >= 0 if already added
	long iFluid = get_Fluid_index(FluidName);
	// If not added yet, and a valid fluid, then continue
	if (iFluid < 0 && set_REFPROP_fluid(FluidName, x)) // If you can set the fluid, it's a valid fluid
	{
		Fluids.add_REFPROP_fluid(FluidName,std::vector<double>(1,1));
		return true;
		std::cout << format("Added the fluid %s\n",FluidName.c_str()).c_str();
	}
	return true;
}
std::string get_index_units(long index)
{
	std::map<long, std::string>::iterator it;
	// Try to find using the map
	it = units_map.find(index);
	// If it is found the iterator will not be equal to end
	if (it != units_map.end() )
	{
		// Return the index of the parameter
		return (*it).second;
	}
	else
	{
		return std::string("Didn't match parameter");
	}
}

std::string get_ASHRAE34(std::string fluid)
{
	long iFluid = get_Fluid_index(fluid);
	if (iFluid > -1)
	{
		Fluid *pFluid = get_fluid(iFluid);
		return pFluid->environment.ASHRAE34;
	}
	else
	{
		return "Fluid name invalid";
	}
}
long get_param_index(std::string param)
{
	std::map<std::string,long>::iterator it;
	// Try to find using the map
	it = param_map.find(param);
	// If it is found the iterator will not be equal to end
	if (it != param_map.end() )
	{
		// Return the index of the parameter
		return (*it).second;
	}
	else
	{
		return -1;
	}
}

static int IsCoolPropFluid(std::string FluidName)
{
	// Try to get the fluid from Fluids by name
	try
	{
		pFluid = Fluids.get_fluid(FluidName);
	}
	catch (NotImplementedError &)
	{
		return false;
	}
	// If NULL, didn't find it (or its alias)
	if (pFluid!=NULL)
	{
		return true;
	}
	else
		return false;
}

static int IsBrine(char* Ref)
{
	// First check whether it is one of the Brines that does
	// not have a pure-fluid equivalent in CoolProp
    if (
        strcmp(Ref,"HC-10")==0 || 
        strncmp(Ref,"PG-",3)==0 ||
		strncmp(Ref,"EG-",3)==0 ||
		strncmp(Ref,"EA-",3)==0 ||
		strncmp(Ref,"MA-",3)==0 ||
		strncmp(Ref,"Glycerol-",9)==0 ||
		strncmp(Ref,"K2CO3-",6)==0 ||
		strncmp(Ref,"CaCl2-",6)==0 ||
		strncmp(Ref,"MgCl2-",6)==0 ||
		strncmp(Ref,"NaCl-",5)==0 ||
		strncmp(Ref,"KAC-",4)==0 ||
		strncmp(Ref,"KFO-",4)==0 ||
		strncmp(Ref,"LiCl-",4)==0 ||
        strncmp(Ref,"NH3/H2O-",8)==0
       )
    {
        return 1;
    }
	// Then check for diluants that are also pure fluids in CoolProp
	else if ( (strncmp(Ref,"Methanol-",9)==0 && Ref[8] == '-') ||
		      (strncmp(Ref,"Ethanol-",8)==0 && Ref[7] == '-') ||
			  (strncmp(Ref,"NH3-",4)==0 && Ref[3] == '-')
		)
	{
		return 1;
	}
    else
    {
        return 0;
    }
}
static int IsREFPROP(std::string Ref)
{
    if (!Ref.compare(0,8,"REFPROP-"))
        return 1;
    else
        return 0;
}

long getFluidType(std::string FluidName){
	if (IsREFPROP(FluidName)) { return FLUID_TYPE_REFPROP;}
	else if(IsIncompressibleLiquid(FluidName)){ return FLUID_TYPE_INCOMPRESSIBLE_LIQUID;}
	// TODO SOLUTION: Check if working
	else if(IsIncompressibleSolution(FluidName)){ return FLUID_TYPE_INCOMPRESSIBLE_SOLUTION; }
	else {
		// Try to get the index of the fluid
		long iFluid = get_Fluid_index(FluidName);
		// If iFluid is greater than -1, it is a CoolProp Fluid, otherwise not
		if (iFluid > -1) {
			// Get a pointer to the fluid object
			pFluid = get_fluid(iFluid);
			if (pFluid->pure())	{ return FLUID_TYPE_PURE;}
			else { return FLUID_TYPE_PSEUDOPURE; }
		} else {
			throw ValueError(format("Bad Fluid name [%s] - not a CoolProp fluid",FluidName.c_str()));
		}
	}
	return -1;
}

EXPORT_CODE int CONVENTION IsFluidType(char *Ref, char *Type)
{
	pFluid = Fluids.get_fluid(Ref);

	if (IsBrine(Ref)){ // TODO Solution: Remove this part
		if (!strcmp(Type,"Brine")){
			return 1;
		}
		else{
			return 0;
		}
	}
	else if (IsIncompressibleSolution(Ref)){
		if (!strcmp(Type,"Solution")){
			return 1;
		}
		else{
			return 0;
		}
	}
	else if (IsIncompressibleLiquid(Ref)){
		if (!strcmp(Type,"Liquid")){
			return 1;
		}
		else{
			return 0;
		}
	}
	else if (IsREFPROP(Ref)){
		if (!strcmp(Type,"PureFluid")){
			return 1;
		}
		else{
			return 0;
		}
	}
	else if (!pFluid->pure()){
		if (!strcmp(Type,"PseudoPure") || !strcmp(Type,"PseudoPureFluid")){
			return 1;
		}
		else{
			return 0;
		}
	}
	else if (pFluid->pure()){
		if (!strcmp(Type,"PureFluid")){
			return 1;
		}
		else{
			return 0;
		}
	}
    else
    {
        return 0;
    }
}

double _T_hp_secant(std::string Ref, double h, double p, double T_guess)
{
    double x1=0,x2=0,x3=0,y1=0,y2=0,eps=1e-8,change=999,f=999,T=300;
    int iter=1;

    while ((iter<=3 || fabs(f)>eps) && iter<100)
    {
        if (iter==1){x1=T_guess; T=x1;}
        if (iter==2){x2=T_guess+0.1; T=x2;}
        if (iter>2) {T=x2;}
            f=Props("H",'T',T,'P',p,Ref)-h;
        if (iter==1){y1=f;}
        if (iter>1)
        {
            y2=f;
            x3=x2-y2/(y2-y1)*(x2-x1);
            change=fabs(y2/(y2-y1)*(x2-x1));
            y1=y2; x1=x2; x2=x3;
        }
        iter=iter+1;
        if (iter>100)
        {
			throw SolutionError(format("iter %d: T_hp not converging with inputs(%s,%g,%g,%g) value: %0.12g\n",iter,(char*)Ref.c_str(),h,p,T_guess,f));
        }
    }
    return T;
}

std::string Phase_Trho(std::string Fluid, double T, double rho)
{
	try{
		// Try to load the CoolProp Fluid
		pFluid = Fluids.get_fluid(Fluid);
		double pL,pV,rhoL,rhoV;
		return pFluid->phase_Trho(T,rho, &pL, &pV, &rhoL, &rhoV);
	}
	catch(NotImplementedError &){
		return std::string("");
	}
	return std::string("");
}

std::string Phase(std::string Fluid, double T, double p)
{
	try{
		// Try to load the CoolProp Fluid
		pFluid = Fluids.get_fluid(Fluid);
		double pL,pV,rhoL,rhoV;
		return pFluid->phase_Tp(T, p, &pL, &pV, &rhoL, &rhoV);
	}
	catch(NotImplementedError &){
		return std::string("");
	}
	return std::string("");
}

std::string Phase_Tp(std::string Fluid, double T, double p)
{
	return Phase(Fluid,T,p);
}


/*
 * Start with the internal functions to handle different inputs
 * First we handle the constants: Props1
 */
// Internal one to do the actual calculations
double _Props1(std::string FluidName, std::string Output)
{
	// Try to load the CoolProp Fluid
	pFluid = Fluids.get_fluid(FluidName);
	if (pFluid != NULL)
	{
		// It's a CoolProp fluid
		// Convert the parameter to integer
		long iOutput = get_param_index(Output);
		if (iOutput < 0){
			throw ValueError(format("Your output key [%s] is not valid. (names are case sensitive)",Output.c_str()));
		}
		// Get the output using the conventional function
		return _CoolProp_Fluid_Props(iOutput,0,0,0,0,pFluid);
	}
	else if (IsREFPROP(FluidName))
	{
		// REFPROP fluid, or something else that is invalid
		try{
			long iOutput = get_param_index(Output);
			switch (iOutput)
			{
				case iTtriple:
				case iTcrit:
				case iPcrit:
				case iTmin:
				case iMM:
				case iRhocrit:
				case iAccentric:
					return Props(Output,'T',0,'P',0,FluidName);
					break;
				default:
					throw ValueError(format("Output parameter \"%s\" is invalid for REFPROP fluid",Output.c_str()));
					break;
			}
		}
		// Catch any error that subclasses the std::exception
		catch(std::exception &e){
			err_string = std::string("CoolProp error: ").append(e.what());
			return _HUGE;
		}
	}
	else
	{
		throw ValueError(format("Fluid \"%s\" is an invalid fluid",FluidName.c_str()));
	}
	return -_HUGE;
}
// Define the functions from the header file
double Props(char *FluidName, char *Output){
	return _Props1(std::string(FluidName), std::string(Output));
}
double Props(std::string FluidName,std::string Output){
	return _Props1(FluidName, Output);
}
double Props1(std::string FluidName,std::string Output){
	return _Props1(FluidName, Output);
}

/*
 * Now we need an internal functions to handle different
 * inputs for non-constant values: Props
 */
// Internal one to do the actual calculations, make this a wrapped function so
// that error bubbling can be done properly
double _Props(std::string Output, std::string Name1, double Prop1, std::string Name2, double Prop2, std::string Ref)
{
	if (get_debug_level()>5){
		std::cout<<__FILE__<<": "<<Output.c_str()<<","<<Name1.c_str()<<","<<Prop1<<","<<Name2.c_str()<<","<<Prop2<<","<<Ref.c_str()<<std::endl;
	}
	/* 
    If the fluid name is not actually a refrigerant name, but a string beginning with "REFPROP-",
    then REFPROP is used to calculate the desired property.
    */
    if (IsREFPROP(Ref))  // First eight characters match "REFPROP-"
    {
    	if (get_debug_level()>7) std::cout<<__FILE__<<": Identified Refprop fluid - "<<Ref.c_str()<<std::endl;
        // Stop here if there is no REFPROP support
    	if (REFPROPFluidClass::refpropSupported()) {
			return REFPROP(Output,Name1,Prop1,Name2,Prop2,Ref);
    	} else {
    		throw AttributeError(format("Your refrigerant [%s] is from REFPROP, but CoolProp does not support REFPROP on this platform, yet.",Ref.c_str()));
    		return -_HUGE;
    	}
    }
	else if (IsCoolPropFluid(Ref))
	{
		if (get_debug_level()>7) std::cout<<__FILE__<<": Identified CoolProp fluid - "<<Ref.c_str()<<std::endl;
		pFluid = Fluids.get_fluid(Ref);
		// Convert all the parameters to integers
		long iOutput = get_param_index(Output);
		if (iOutput<0) 
			throw ValueError(format("Your output key [%s] is not valid. (names are case sensitive)",Output.c_str()));
		long iName1 = get_param_index(std::string(Name1));  
		if (iName1<0) 
			throw ValueError(format("Your input key #1 [%s] is not valid. (names are case sensitive)",Name1.c_str()));
		long iName2 = get_param_index(std::string(Name2));  
		if (iName2<0) 
			throw ValueError(format("Your input key #2 [%s] is not valid. (names are case sensitive)",Name2.c_str()));
		// Call the internal method that uses the parameters converted to longs
		return _CoolProp_Fluid_Props(iOutput,iName1,Prop1,iName2,Prop2,pFluid);
	}
    // It's a brine, call the brine routine // TODO Solutions: remove this part
	else if (IsBrine((char*)Ref.c_str()))
    {
		if (get_debug_level()>7) std::cout<<__FILE__<<": Identified brine - "<<Ref.c_str()<<std::endl;
		//Enthalpy and pressure are the inputs
		if ((Name1.c_str()[0]=='H' && Name2.c_str()[0]=='P') || (Name2.c_str()[0]=='H' && Name1.c_str()[0]=='P'))
        {
			if (Name2.c_str()[0]=='H' && Name1.c_str()[0]=='P')
			{
				std::swap(Prop1,Prop2);
				std::swap(Name1,Name2);
			}
			// Start with a guess of 10 K below max temp of fluid
			double Tguess = SecFluids('M',Prop1,Prop2,(char*)Ref.c_str())-10;
			// Solve for the temperature
			double T =_T_hp_secant(Ref,Prop1,Prop2,Tguess);
			// Return whatever property is desired
			return SecFluids(Output[0],T,Prop2,(char*)Ref.c_str());
		}
		else if ((Name1.c_str()[0] == 'T' && Name2.c_str()[0] =='P') || (Name1.c_str()[0] == 'P' && Name2.c_str()[0] == 'T'))
        {
			if (Name1.c_str()[0] =='P' && Name2.c_str()[0] =='T'){
				std::swap(Prop1,Prop2);
			}
			return SecFluids(Output[0],Prop1,Prop2,(char*)Ref.c_str());
		}
		else
		{
			throw ValueError("For brine, inputs must be (order does not matter) 'T' and 'P', or 'H' and 'P'");
		}
    }
	// It's an incompressible liquid, call the routine
	else if (IsIncompressibleLiquid(Ref))
    {
		if (get_debug_level()>7) std::cout<<__FILE__<<": Identified incompressible liquid - "<<Ref.c_str()<<std::endl;
		//Enthalpy and pressure are the inputs
		if ((Name1.c_str()[0]=='H' && Name2.c_str()[0]=='P') || (Name2.c_str()[0]=='H' && Name1.c_str()[0]=='P'))
        {
			if (Name2.c_str()[0]=='H' && Name1.c_str()[0]=='P')
			{
				std::swap(Prop1,Prop2);
				std::swap(Name1,Name2);
			}
			
			// Solve for the temperature
			double Tma     = IncompLiquid(get_param_index(std::string("Tmax")),0.0,0.0,Ref);
			double T_guess = Tma - 10.0 ;
			double T =_T_hp_secant(Ref,Prop1,Prop2,T_guess);
			// Return whatever property is desired
			return IncompLiquid(get_param_index(Output),T,Prop2,Ref);
		}
		else if ((Name1.c_str()[0] == 'T' && Name2.c_str()[0] =='P') || (Name1.c_str()[0] == 'P' && Name2.c_str()[0] == 'T'))
        {
			if (Name1.c_str()[0] =='P' && Name2.c_str()[0] =='T'){
				std::swap(Prop1,Prop2);
			}
			return IncompLiquid(get_param_index(Output),Prop1,Prop2,Ref);
		}
		else
		{
			throw ValueError("For incompressible fluids, inputs must be (order does not matter) 'T' and 'P', or 'H' and 'P'");
		}
    }
    // It's an incompressible solution, call the routine
	else if (IsIncompressibleSolution(Ref))
	{
		if (get_debug_level()>7) std::cout<<__FILE__<<": Identified incompressible solution - "<<Ref.c_str()<<std::endl;
		//Enthalpy and pressure are the inputs
		if ((Name1.c_str()[0]=='H' && Name2.c_str()[0]=='P') || (Name2.c_str()[0]=='H' && Name1.c_str()[0]=='P'))
		{
			if (Name2.c_str()[0]=='H' && Name1.c_str()[0]=='P')
			{
				std::swap(Prop1,Prop2);
				std::swap(Name1,Name2);
			}

			// Solve for the temperature
			double Tma     = IncompSolution(get_param_index(std::string("Tmax")),0.0,0.0,Ref);
			double Tmi     = IncompSolution(get_param_index(std::string("Tmin")),0.0,0.0,Ref);
			double T_guess = (Tma+Tmi)/2.0 ;
			double T =_T_hp_secant(Ref,Prop1,Prop2,T_guess);
			// Return whatever property is desired
			return IncompSolution(get_param_index(Output),T,Prop2,Ref);
		}
		else if ((Name1.c_str()[0] == 'T' && Name2.c_str()[0] =='P') || (Name1.c_str()[0] == 'P' && Name2.c_str()[0] == 'T'))
		{
			if (Name1.c_str()[0] =='P' && Name2.c_str()[0] =='T'){
				std::swap(Prop1,Prop2);
			}
			return IncompSolution(get_param_index(Output),Prop1,Prop2,Ref);
		}
		else
		{
			throw ValueError("For incompressible solutions, inputs must be (order does not matter) 'T' and 'P', or 'H' and 'P'");
		}
	}
	else
	{
		throw ValueError(format("Your fluid name [%s] is not a CoolProp fluid, a REFPROP fluid, a brine or a liquid",Ref.c_str()));
	}
}
double _CoolProp_Fluid_Props(long iOutput, long iName1, double Prop1, long iName2, double Prop2, Fluid *pFluid, bool SinglePhase)
{
	double val = _HUGE;
	// This private method uses the indices directly for speed

	if (get_debug_level()>3){
		std::cout<<__FILE__<<" _CoolProp_Fluid_Props: "<<iOutput<<","<<iName1<<","<<Prop1<<","<<iName2<<","<<Prop2<<","<<pFluid->get_name().c_str()<<std::endl;
	}

	// Generate a State instance wrapped around the Fluid instance
	CoolPropStateClass CPS(pFluid);

	// Check if it is an output that doesn't require a state input
	// Deal with it and return

	switch (iOutput)
	{
		case iMM:
		case iPcrit:
		case iTcrit:
		case iTtriple:
		case iPtriple:
		case iRhocrit:
		case iTmin:
		case iAccentric:
		case iPHASE_LIQUID:
		case iPHASE_GAS:
		case iPHASE_SUPERCRITICAL:
		case iPHASE_TWOPHASE:
		case iGWP20:
		case iGWP100:
		case iGWP500:
		case iODP:
		case iCritSplineT:
		case iScrit:
		case iHcrit:
		case iTreduce:
		case iRhoreduce:
			return CPS.keyed_output(iOutput);
	}

	// Update the class
	CPS.update(iName1,Prop1,iName2,Prop2);
	// Get the output
	val = CPS.keyed_output(iOutput);
	
	if (get_debug_level()>5){
		std::cout<<__FILE__<<" _CoolProp_Fluid_Props return: "<<val<<std::endl;
	}
	// Return the value
	return val;
}
EXPORT_CODE double CONVENTION IProps(long iOutput, long iName1, double Prop1, long iName2, double Prop2, long iFluid)
{
	pFluid = Fluids.get_fluid(iFluid);
	// Didn't work
	if (pFluid == NULL){
		err_string=std::string("CoolProp error: ").append(format("%d is an invalid fluid index to IProps",iFluid));
		return _HUGE;
	}
	else{
		// In this function the error catching happens;
		try{
			// This is already converted to the right units since it comes from the CoolPropStateClass which
			// does the unit handling
			return _CoolProp_Fluid_Props(iOutput,iName1,Prop1,iName2,Prop2,pFluid);
		}
		catch(std::exception &e){
			err_string=std::string("CoolProp error: ").append(e.what());
			return _HUGE;
		}
		catch(...){
			err_string=std::string("CoolProp error: Indeterminate error");
			return _HUGE;
		}
	}
}
double Props(char Output,char Name1, double Prop1, char Name2, double Prop2, char* Ref)
{
	// Go to the std::string, std::string version
	return Props(std::string(1,Output),Name1,Prop1,Name2,Prop2,std::string(Ref));
}
//double Props(char* Output,char Name1, double Prop1, char Name2, double Prop2, char* Ref)
//{
//	// Go to the std::string, std::string version
//	return Props(std::string(Output),Name1,Prop1,Name2,Prop2,std::string(Ref));
//}
double Props(std::string Output,char Name1, double Prop1, char Name2, double Prop2, std::string Ref)
{
	// In this function the error catching happens;
	try{
		return _Props(Output,std::string(1,Name1),Prop1,std::string(1,Name2),Prop2,Ref);
	}
	catch(const CoolPropBaseError& e){
		err_string = std::string("CoolProp error: ").append(e.what());
		return _HUGE;
	}
	catch(const std::exception& e){
			err_string = std::string("CoolProp error: ").append(e.what());
			return _HUGE;
		}
	catch(...){
		err_string = std::string("CoolProp error: Indeterminate error");
		return _HUGE;
	}
	return _HUGE;
}


/// Calculate some interesting derivatives
double _CoolProp_Deriv_Terms(long iTerm, double T, double rho, Fluid * pFluid)
{
	double val = _HUGE;
	// This private method uses the indices directly for speed

	if (get_debug_level()>3){
		std::cout<<__FILE__<<" _CoolProp_Deriv_Terms return: "<<val<<std::endl;
	}

	// Generate a State instance wrapped around the Fluid instance
	CoolPropStateClass CPS(pFluid);

	// Update the class
	CPS.update(iT,T,iD,rho);

	switch (iTerm) {
	case iDERdh_dp__rho:
		case iDERdh_dp__v:
		case iDERZ:
		case iDERdZ_dDelta:
		case iDERdZ_dTau:
		case iDERB:
		case iDERdB_dT:
		case iDERC:
		case iDERdC_dT:
		case iDERphir:
		case iDERdphir_dTau:
		case iDERdphir_dDelta:
		case iDERd2phir_dTau2:
		case iDERd2phir_dDelta2:
		case iDERd2phir_dDelta_dTau:
		case iDERd3phir_dDelta3:
		case iDERd3phir_dDelta2_dTau:
		case iDERd3phir_dDelta_dTau2:
		case iDERd3phir_dTau3:
		case iDERphi0:
		case iDERdphi0_dTau:
		case iDERd2phi0_dTau2:
		case iDERdphi0_dDelta:
		case iDERd2phi0_dDelta2:
		case iDERd2phi0_dDelta_dTau:
		case iDERd3phi0_dTau3:
		case iDERdp_dT__rho:
		case iDERdp_drho__T:
		case iDERdh_dT__rho:
		case iDERdh_drho__T:
		case iDERdrho_dT__p:
		case iDERdrho_dh__p:
		case iDERdrho_dp__h:
		case iDERrho_smoothed:
		case iDERdrho_smoothed_dh:
		case iDERdrho_smoothed_dp:
		case iDERdrhodh_constp_smoothed:
		case iDERdrhodp_consth_smoothed:
		case iDERIsothermalCompressibility:
			val = CPS.keyed_output(iTerm);
			break;
		default:
			throw ValueError(format("Sorry DerivTerms is a work in progress, your derivative term [%d] is not available!",iTerm));
	}

	if (get_debug_level()>5){
		std::cout<<__FILE__<<" _CoolProp_Deriv_Terms return: "<<val<<std::endl;
	}
	// Return the value
	return val;
}

// Define the functions from the header
double DerivTerms(long iTerm, double T, double rho, Fluid * pFluid){
	return _CoolProp_Deriv_Terms(iTerm,T,rho,pFluid);
}
double DerivTerms(std::string Term, double T, double rho, std::string Fluidname){
	if (get_debug_level()>5){
			std::cout<<__FILE__<<": "<<Term.c_str()<<",T="<<T<<",rho="<<rho<<","<<Fluidname.c_str()<<std::endl;
		}
		/*
	    Derivatives are only supported for CoolProp fluids
	    */
	    if (IsCoolPropFluid(Fluidname))
		{
			pFluid = Fluids.get_fluid(Fluidname);
			// for compatibility, replace B and C with VB and VC
			if ((!Term.compare("B")) || (!Term.compare("C"))) {
				Term = std::string("V").append(Term);
			}
			// Convert all the parameters to integers
			long iOutput = get_param_index(Term);
			if (iOutput<0)
				throw ValueError(format("Your output key [%s] is not valid. (names are case sensitive)",Term.c_str()));

			if (T<=0)
				throw ValueError(format("Your input temperature [%f] is not valid.",T));

			if (rho<=0)
				throw ValueError(format("Your input density [%f] is not valid.",rho));
			// Call the internal method that uses the parameters converted to longs
			return _CoolProp_Deriv_Terms(iOutput,T,rho,pFluid);
		}
		else
		{
			throw ValueError(format("Your fluid name [%s] is not a CoolProp fluid.",Fluidname.c_str()));
		}
}
//double DerivTerms(char *Term, double T, double rho, char * FluidName)
//{
//	return DerivTerms(std::string(Term),T,rho,std::string(FluidName));
//}



int set_reference_stateS(std::string Ref, std::string reference_state)
{
	Fluid *pFluid=Fluids.get_fluid(Ref);
	if (pFluid!=NULL)
	{
		return set_reference_stateP(pFluid, reference_state);
	}
	else{
		return -1;
	}
}

int set_reference_stateP(Fluid *pFluid, std::string reference_state)
{
	CoolPropStateClassSI CPS(pFluid);
	if (!reference_state.compare("IIR"))
	{
		CoolPropStateClassSI CPS(pFluid);
		CPS.update(iT,273.15,iQ,0);
		// Get current values for the enthalpy and entropy
		double h1 = CPS.h();
		double s1 = CPS.s();
		double deltah = h1-200000; // offset from 200 kJ/kg enthalpy
		double deltas = s1-1000; // offset from 1 kJ/kg/K entropy
		double delta_a1 = deltas/((8314.472/pFluid->params.molemass));
		double delta_a2 = -deltah/((8314.472/pFluid->params.molemass)*pFluid->reduce.T);
		pFluid->phi0list.push_back(new phi0_enthalpy_entropy_offset(delta_a1, delta_a2));
		return 0;
	}
	else if (!reference_state.compare("ASHRAE"))
	{
		CoolPropStateClassSI CPS(pFluid);
		CPS.update(iT,233.15,iQ,0);
		// Get current values for the enthalpy and entropy
		double h1 = CPS.h();
		double s1 = CPS.s();
		double deltah = h1-0; // offset from 0 kJ/kg enthalpy
		double deltas = s1-0; // offset from 0 kJ/kg/K entropy
		double delta_a1 = deltas/((8314.472/pFluid->params.molemass));
		double delta_a2 = -deltah/((8314.472/pFluid->params.molemass)*pFluid->reduce.T);
		pFluid->phi0list.push_back(new phi0_enthalpy_entropy_offset(delta_a1, delta_a2));
		return 0;
	}
	else if (!reference_state.compare("NBP"))
	{
		CoolPropStateClassSI CPS(pFluid);
		CPS.update(iP,101325.0,iQ,0); // Saturated boiling point at 1 atmosphere
		// Get current values for the enthalpy and entropy
		double h1 = CPS.h();
		double s1 = CPS.s();
		double deltah = h1-0; // offset from 0 kJ/kg enthalpy
		double deltas = s1-0; // offset from 0 kJ/kg/K entropy
		double delta_a1 = deltas/((8314.472/pFluid->params.molemass));
		double delta_a2 = -deltah/((8314.472/pFluid->params.molemass)*pFluid->reduce.T);
		pFluid->phi0list.push_back(new phi0_enthalpy_entropy_offset(delta_a1, delta_a2));
		return 0;
	}
	else
	{ 
		return -1;
	}

}
int set_reference_stateD(std::string Ref, double T, double rho, double h0, double s0)
{
	pFluid=Fluids.get_fluid(Ref);
	if (pFluid!=NULL)
	{
		CoolPropStateClassSI CPS(pFluid);
		CPS.update(iT,T,iD,rho);
		// Get current values for the enthalpy and entropy
		double h1 = CPS.h();
		double s1 = CPS.s();
		double deltah = h1-h0; // offset from given enthalpy in SI units
		double deltas = s1-s0; // offset from given enthalpy in SI units
		double delta_a1 = deltas/((8314.472/pFluid->params.molemass));
		double delta_a2 = -deltah/((8314.472/pFluid->params.molemass)*pFluid->reduce.T);
		pFluid->phi0list.push_back(new phi0_enthalpy_entropy_offset(delta_a1, delta_a2));
		return 0;
	}
	else{
		return -1;
	}
}

std::string get_BibTeXKey(std::string Ref, std::string item)
{
	pFluid=Fluids.get_fluid(Ref);
	if (pFluid!=NULL)
	{
		
		if (!item.compare("EOS")){ return pFluid->BibTeXKeys.EOS; }
		else if (!item.compare("CP0")){ return pFluid->BibTeXKeys.CP0; }
		else if (!item.compare("VISCOSITY")){ return pFluid->BibTeXKeys.VISCOSITY; }
		else if (!item.compare("CONDUCTIVITY")){ return pFluid->BibTeXKeys.CONDUCTIVITY; }
		else if (!item.compare("ECS_LENNARD_JONES")){ return pFluid->BibTeXKeys.ECS_LENNARD_JONES; }
		else if (!item.compare("ECS_FITS")){ return pFluid->BibTeXKeys.ECS_FITS; }
		else if (!item.compare("SURFACE_TENSION")){ return pFluid->BibTeXKeys.SURFACE_TENSION; }
		else{ return "Bad key";}
	}
	else{
		return std::string("");
	}
}

std::string get_global_param_string(std::string ParamName)
{
	if (!ParamName.compare("version"))
	{
		return std::string(version);
	}
	else if (!ParamName.compare("errstring"))
	{
		std::string temp = err_string;
		err_string = std::string("");
		return temp;
	}
	else if (!ParamName.compare("warnstring"))
	{
		std::string temp = warning_string;
		warning_string = std::string("");
		return temp;
	}
	else if (!ParamName.compare("gitrevision"))
	{
		return gitrevision;
	}
	else if (!ParamName.compare("FluidsList") || !ParamName.compare("fluids_list"))
	{
		return Fluids.FluidList();
	}
	else
	{
		return format("Input value [%s] is invalid",ParamName.c_str()).c_str();
	}
};
std::string get_fluid_param_string(std::string FluidName, std::string ParamName)
{
	try{
		pFluid = Fluids.get_fluid(FluidName);
		// Didn't work
		if (pFluid == NULL){
			err_string=std::string("CoolProp error: ").append(format("%s is an invalid fluid for get_fluid_param_string",FluidName.c_str()));
			return format("%s is an invalid fluid for get_fluid_param_string",FluidName.c_str()).c_str();
		}
		else{
			if (!ParamName.compare("aliases"))
			{
				std::vector<std::string> v = pFluid->get_aliases();
				return strjoin(v,", ");
			}
			else if (!ParamName.compare("CAS") || !ParamName.compare("CAS_number"))
			{
				return pFluid->params.CAS;
			}
			else if (!ParamName.compare("ASHRAE34"))
			{
				return pFluid->environment.ASHRAE34;
			}
			else if (!ParamName.compare("REFPROPName") || !ParamName.compare("REFPROP_name"))
			{
				return pFluid->get_REFPROPname();
			}
			else if (!ParamName.compare("TTSE_mode"))
			{
				int mode = pFluid->TTSESinglePhase.get_mode();
				switch (mode)
				{
				case TTSE_MODE_TTSE:
					return "TTSE";
				case TTSE_MODE_BICUBIC:
					return "BICUBIC";
				default:
					throw ValueError("TTSE mode is invalid");
				}
			}
			else
			{
				return format("Input value [%s] is invalid for Fluid [%s]",ParamName.c_str(),FluidName.c_str()).c_str();
			}
		}
	}
	catch(std::exception &e)
	{
		return(std::string("CoolProp error: ").append(e.what()));
	}
	catch(...){
		return(std::string("CoolProp error: Indeterminate error"));
	}
}



