#include "CoolProp.h"
#include "AbstractState.h"

#if defined(_MSC_VER)
#define _CRTDBG_MAP_ALLOC
#define _CRT_SECURE_NO_WARNINGS
#include <crtdbg.h>
#endif

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
#include <string>
#include "CoolPropTools.h"
#include "Solvers.h"

namespace CoolProp
{
///// The lower-level methods that can throw exceptions
//double _CoolProp_Fluid_PropsSI(long iOutput, long iName1, double Value1, long iName2, double Value2, Fluid *pFluid);
double _PropsSI(std::string &Output,std::string &Name1, double Prop1, std::string &Name2, double Prop2, std::string &Ref);
//double _Props1SI(std::string FluidName, std::string Output);

static std::string err_string;
static std::string warning_string;

//// This is very hacky, but pull the git revision from the file
//#include "gitrevision.h" // Contents are like "long gitrevision = "aa121435436ggregrea4t43t433";"
//#include "version.h" // Contents are like "char version [] = "2.5";"

//int global_Phase = -1;
//bool global_SinglePhase = false;
//bool global_SaturatedL = false;
//bool global_SaturatedV = false;

//void set_warning(std::string warning){ warning_string = warning; }
//
//FluidsContainer Fluids = FluidsContainer();
//
//void set_err_string(std::string error_string)
//{
//	err_string = error_string;
//}
//
//void set_phase(std::string Phase_str){
//	if (!Phase_str.compare("Two-Phase")){
//		global_SinglePhase = false;
//		global_SaturatedL = false;
//		global_SaturatedV = false;
//		global_Phase = iTwoPhase;
//	}
//	else if (!Phase_str.compare("Liquid")){
//		global_SinglePhase = true;
//		global_SaturatedL = false;
//		global_SaturatedV = false;
//		global_Phase = iLiquid;
//	}
//	else if (!Phase_str.compare("Gas")){
//		global_SinglePhase = true;
//		global_SaturatedL = false;
//		global_SaturatedV = false;
//		global_Phase = iGas;
//	}
//	else if (!Phase_str.compare("Supercritical")){
//		global_SinglePhase = true;
//		global_SaturatedL = false;
//		global_SaturatedV = false;
//		global_Phase = iSupercritical;
//	}
//	else if (!Phase_str.compare("SaturatedL")){
//		global_SinglePhase = false;
//		global_SaturatedL = true;
//		global_SaturatedV = false;
//		global_Phase = iTwoPhase;
//	}
//	else if (!Phase_str.compare("SaturatedV")){
//		global_SinglePhase = false;
//		global_SaturatedL = false;
//		global_SaturatedV = true;
//		global_Phase = iTwoPhase;
//	}
//}
//
//long get_Fluid_index(std::string FluidName)
//{
//	// Try to get the fluid from Fluids by name
//	pFluid = Fluids.get_fluid(FluidName);
//	// If NULL, didn't find it (or its alias)
//	if (pFluid!=NULL)
//	{
//		// Find the fluid index
//		return Fluids.get_fluid_index(pFluid);
//	}
//	else
//		return -1;
//}
//
//long get_param_index(std::string param)
//{
//	std::map<std::string,long>::iterator it;
//	// Try to find using the map
//	it = param_map.find(param);
//	// If it is found the iterator will not be equal to end
//	if (it != param_map.end() )
//	{
//		// Return the index of the parameter
//		return (*it).second;
//	}
//	else
//	{
//		return -1;
//	}
//}
//
//static int IsCoolPropFluid(std::string FluidName)
//{
//	// Try to get the fluid from Fluids by name
//	try
//	{
//		pFluid = Fluids.get_fluid(FluidName);
//	}
//	catch (NotImplementedError &)
//	{
//		return false;
//	}
//	// If NULL, didn't find it (or its alias)
//	if (pFluid!=NULL)
//	{
//		return true;
//	}
//	else
//		return false;
//}
//
//static int IsBrine(const char* Ref)
//{
//	// First check whether it is one of the Brines that does
//	// not have a pure-fluid equivalent in CoolProp
//    if (
//        strcmp(Ref,"HC-10")==0 || 
//        strncmp(Ref,"PG-",3)==0 ||
//		strncmp(Ref,"EG-",3)==0 ||
//		strncmp(Ref,"EA-",3)==0 ||
//		strncmp(Ref,"MA-",3)==0 ||
//		strncmp(Ref,"Glycerol-",9)==0 ||
//		strncmp(Ref,"K2CO3-",6)==0 ||
//		strncmp(Ref,"CaCl2-",6)==0 ||
//		strncmp(Ref,"MgCl2-",6)==0 ||
//		strncmp(Ref,"NaCl-",5)==0 ||
//		strncmp(Ref,"KAC-",4)==0 ||
//		strncmp(Ref,"KFO-",4)==0 ||
//		strncmp(Ref,"LiCl-",4)==0 ||
//        strncmp(Ref,"NH3/H2O-",8)==0
//       )
//    {
//        return 1;
//    }
//	// Then check for diluants that are also pure fluids in CoolProp
//	else if ( (strncmp(Ref,"Methanol-",9)==0 && Ref[8] == '-') ||
//		      (strncmp(Ref,"Ethanol-",8)==0 && Ref[7] == '-') ||
//			  (strncmp(Ref,"NH3-",4)==0 && Ref[3] == '-')
//		)
//	{
//		return 1;
//	}
//    else
//    {
//        return 0;
//    }
//}
//
//long getFluidType(std::string FluidName){
//	if (IsREFPROP(FluidName)) { return FLUID_TYPE_REFPROP;}
//	else if(IsIncompressibleLiquid(FluidName)){ return FLUID_TYPE_INCOMPRESSIBLE_LIQUID;}
//	// TODO SOLUTION: Check if working
//	else if(IsIncompressibleSolution(FluidName)){ return FLUID_TYPE_INCOMPRESSIBLE_SOLUTION; }
//	else {
//		// Try to get the index of the fluid
//		long iFluid = get_Fluid_index(FluidName);
//		// If iFluid is greater than -1, it is a CoolProp Fluid, otherwise not
//		if (iFluid > -1) {
//			// Get a pointer to the fluid object
//			pFluid = get_fluid(iFluid);
//			if (pFluid->pure())	{ return FLUID_TYPE_PURE;}
//			else { return FLUID_TYPE_PSEUDOPURE; }
//		} else {
//			throw ValueError(format("Bad Fluid name [%s] - not a CoolProp fluid",FluidName.c_str()));
//		}
//	}
//	return -1;
//}
//
//EXPORT_CODE int CONVENTION IsFluidType(const char *Ref, const char *Type)
//{
//	pFluid = Fluids.get_fluid(Ref);
//
//	if (IsBrine(Ref)){ // TODO Solution: Remove this part
//		if (!strcmp(Type,"Brine")){
//			return 1;
//		}
//		else{
//			return 0;
//		}
//	}
//	else if (IsIncompressibleSolution(Ref)){
//		if (!strcmp(Type,"Solution")){
//			return 1;
//		}
//		else{
//			return 0;
//		}
//	}
//	else if (IsIncompressibleLiquid(Ref)){
//		if (!strcmp(Type,"Liquid")){
//			return 1;
//		}
//		else{
//			return 0;
//		}
//	}
//	else if (IsREFPROP(Ref)){
//		if (!strcmp(Type,"PureFluid")){
//			return 1;
//		}
//		else{
//			return 0;
//		}
//	}
//	else if (!pFluid->pure()){
//		if (!strcmp(Type,"PseudoPure") || !strcmp(Type,"PseudoPureFluid")){
//			return 1;
//		}
//		else{
//			return 0;
//		}
//	}
//	else if (pFluid->pure()){
//		if (!strcmp(Type,"PureFluid")){
//			return 1;
//		}
//		else{
//			return 0;
//		}
//	}
//    else
//    {
//        return 0;
//    }
//}
//
//
//std::string Phase_Trho(std::string Fluid, double T, double rho)
//{
//	try{
//		// Try to load the CoolProp Fluid
//		pFluid = Fluids.get_fluid(Fluid);
//		double pL,pV,rhoL,rhoV;
//		return pFluid->phase_Trho(T,rho, pL, pV, rhoL, rhoV);
//	}
//	catch(NotImplementedError &){
//		return std::string("");
//	}
//	return std::string("");
//}
//
//std::string Phase(std::string Fluid, double T, double p)
//{
//	try{
//		// Try to load the CoolProp Fluid
//		pFluid = Fluids.get_fluid(Fluid);
//		double pL,pV,rhoL,rhoV;
//		return pFluid->phase_Tp(T, p, pL, pV, rhoL, rhoV);
//	}
//	catch(NotImplementedError &){
//		return std::string("");
//	}
//	return std::string("");
//}
//
//std::string Phase_Tp(std::string Fluid, double T, double p)
//{
//	return Phase(Fluid,T,p);
//}
//
///*
// * Start with the internal functions to handle different inputs
// * First we handle the constants: Props1
// */
//// Internal one to do the actual calculations
//double _Props1SI(std::string FluidName, std::string Output)
//{
//    double out = _HUGE;
//	// Try to load the CoolProp Fluid
//	pFluid = Fluids.get_fluid(FluidName);
//	if (pFluid != NULL)
//	{
//		// It's a CoolProp fluid
//		// Convert the parameter to integer
//		long iOutput = get_param_index(Output);
//		if (iOutput < 0){
//			throw ValueError(format("Your output key [%s] is not valid. (names are case sensitive)",Output.c_str()));
//		}
//		// Get the output using the conventional function
//		return _CoolProp_Fluid_PropsSI(iOutput,0,0,0,0,pFluid);
//	}
//	else if (IsREFPROP(FluidName))
//	{
//		// REFPROP fluid
//		long iOutput = get_param_index(Output);
//		switch (iOutput)
//		{
//			case iTtriple:
//			case iTcrit:
//			case iPcrit:
//			case iTmin:
//			case iMM:
//			case iRhocrit:
//			case iAccentric:
//				out = Props(Output,'T',0,'P',0,FluidName);
//				break;
//			default:
//				throw ValueError(format("Output parameter \"%s\" is invalid for REFPROP fluid",Output.c_str()));
//		}
//        return convert_from_SI_to_unit_system(iOutput,out,get_standard_unit_system());
//	}
//	else
//	{
//		throw ValueError(format("Fluid \"%s\" is an invalid fluid",FluidName.c_str()));
//	}
//	return -_HUGE;
//}
//double _Props1(std::string FluidName, std::string Output)
//{
//    double val = _Props1SI(FluidName, Output);
//    if (ValidNumber(val))
//    {
//        long iOutput = get_param_index(Output);
//        return convert_from_SI_to_unit_system(iOutput,val,get_standard_unit_system());
//    }
//    else
//    {
//        return _HUGE;
//    }
//}
//// Define the functions from the header file
//double Props1SI(std::string FluidName,std::string Output){
//    // Redirect to the Props() function that takes const char *
//	// In this function the error catching happens;
//	try{
//		return _Props1SI(FluidName, Output);
//	}
//	catch(const std::exception& e){
//			err_string = std::string("CoolProp error: ").append(e.what());
//			return _HUGE;
//		}
//	catch(...){
//		err_string = std::string("CoolProp error: Indeterminate error");
//		return _HUGE;
//	}
//	return _HUGE;
//}
//double Props1(std::string FluidName, std::string Output){
//
//    // Redirect to the Props() function that takes const char *
//	// In this function the error catching happens;
//	try{
//		return _Props1(FluidName, Output);
//	}
//	catch(const std::exception& e){
//			err_string = std::string("CoolProp error: ").append(e.what());
//			return _HUGE;
//		}
//	catch(...){
//		err_string = std::string("CoolProp error: Indeterminate error");
//		return _HUGE;
//	}
//	return _HUGE;
//}

double PropsSI(std::string Output, std::string Name1, double Prop1, std::string Name2, double Prop2, std::string Ref)
{
    // In this function the error catching happens;
	try{
		return _PropsSI(Output,Name1,Prop1,Name2,Prop2,Ref);
	}
	catch(const std::exception& e){
        std::cout << e.what() << std::endl;
        //set_err_string(e.what());
		return _HUGE;
	}
	catch(...){
		//set_err_string(std::string("CoolProp error: Indeterminate error"));
		return _HUGE;
	}
	return _HUGE;
}

// Internal function to do the actual calculations, make this a wrapped function so
// that error bubbling can be done properly
double _PropsSI(std::string &Output, std::string &Name1, double Prop1, std::string &Name2, double Prop2, std::string &Ref)
{
    double x1, x2;
	/*if (get_debug_level()>5){
		std::cout << format("%s:%d: _Props(%s,%s,%g,%s,%g,%s)\n",__FILE__,__LINE__,Output.c_str(),Name1.c_str(),Prop1,Name2.c_str(),Prop2,Ref.c_str()).c_str();
	}*/

    // Convert all the input and output parameters to integers
	long iOutput = get_parameter_index(Output);
	long iName1 = get_parameter_index(Name1);
	long iName2 = get_parameter_index(Name2);

    // The state we are going to use
    AbstractState *State = NULL;
    try
    {
        // Generate the State class pointer using the factory function
        // TODO: need to parse ref string to obtain the backend, or use "?" if unknown
        State = AbstractState::factory("HEOS", Ref);

        // Obtain the input pair
        long pair = generate_update_pair(iName1, Prop1, iName2, Prop2, x1, x2);

        // First check if it is a trivial input (critical/max parameters for instance)
        // TODO: check for trivial inputs that do not require the use of the eos
        /*if (State->is_trivial_output(iOutput))
        { 
            double val = State->trivial_keyed_output(iOutput);
            delete(State);
            return val;
        };*/

        // Update the state
        State->update(pair, x1, x2);

        // Return the desired output
        return State->keyed_output(iOutput);
    }
    catch(...){	
        delete(State); throw;
	}
}

 ///* 
 //   If the fluid name is not actually a refrigerant name, but a string beginning with "REFPROP-",
 //   then REFPROP is used to calculate the desired property.
 //   */
 //   if (IsREFPROP(Ref))  // First eight characters match "REFPROP-"
 //   {
 //   	if (get_debug_level()>7) std::cout<<__FILE__<<": Identified Refprop fluid - "<<Ref.c_str()<<std::endl;
 //       // Stop here if there is no REFPROP support
 //   	if (REFPROPFluidClass::refpropSupported()) {
	//		return REFPROPSI(iOutput,iName1,Prop1,iName2,Prop2,Ref);
 //   	} else {
 //   		throw AttributeError(format("Your fluid [%s] is from REFPROP, but CoolProp does not support REFPROP on this platform, yet.",Ref.c_str()));
 //   		return -_HUGE;
 //   	}
 //   }
	//else if (IsCoolPropFluid(Ref))
	//{
	//	if (get_debug_level()>7) std::cout << format("%s:%d: Identified CoolProp fluid - %s\n",__FILE__,__LINE__,Ref.c_str()).c_str();
	//	pFluid = Fluids.get_fluid(Ref);
	//	// Call the internal method that uses the parameters converted to longs
	//	return _CoolProp_Fluid_PropsSI(iOutput,iName1,Prop1,iName2,Prop2,pFluid);
	//}

 //   // It's a brine, call the brine routine // TODO Solutions: remove this part
	//else if (IsBrine(Ref.c_str()))
 //   {
	//	if (get_debug_level()>7) std::cout<<__FILE__<<": Identified brine - "<<Ref.c_str()<<std::endl;
	//	//Enthalpy and pressure are the inputs
	//	if ((iName1 == iH && iName2 == iP) || (iName2 == iH && iName1 == iP))
 //       {
	//		if (iName2 == iH && iName1 == iP)
	//		{
	//			std::swap(Prop1,Prop2);
	//			std::swap(Name1,Name2);
	//		}
	//		// Start with a guess of 10 K below max temp of fluid
	//		double Tguess = SecFluids("Tmax",Prop1,Prop2,Ref)-10;
	//		// Solve for the temperature
	//		double T =_T_hp_secant(Ref,Prop1,Prop2,Tguess);
	//		// Return whatever property is desired
	//		return SecFluidsSI(Output,T,Prop2,Ref);
	//	}
	//	else if ((iName1 == iT && iName2 == iP) || (iName1 == iP && iName2 == iT))
 //       {
	//		if (iName1 == iP && iName2 == iT){
	//			std::swap(Prop1,Prop2);
	//		}
	//		return SecFluidsSI(Output,Prop1,Prop2,Ref);
	//	}
	//	else
	//	{
	//		throw ValueError("For brine, inputs must be (order does not matter) 'T' and 'P', or 'H' and 'P'");
	//	}
 //   }

	//// It's an incompressible liquid, call the routine
	//else if (IsIncompressibleLiquid(Ref))
 //   {
	//	if (get_debug_level()>7) std::cout<<__FILE__<<": Identified incompressible liquid - "<<Ref.c_str()<<std::endl;
	//	//Enthalpy and pressure are the inputs
	//	if ((iName1 == iH && iName2 == iP) || (iName2 == iH && iName1 == iP))
 //       {
	//		if (iName2 == iH && iName1 == iP)
	//		{
	//			std::swap(Prop1,Prop2);
	//			std::swap(Name1,Name2);
	//		}
	//		
	//		// Solve for the temperature
	//		double Tma     = IncompLiquidSI(get_param_index("Tmax"),0.0,0.0,Ref);
	//		double T_guess = Tma - 10.0 ;
	//		double T =_T_hp_secant(Ref,Prop1,Prop2,T_guess);
	//		// Return whatever property is desired
	//		return IncompLiquidSI(iOutput,T,Prop2,Ref);
	//	}
	//	else if ((iName1 == iT && iName2 == iP) || (iName1 == iP && iName2 == iT))
 //       {
	//		if (iName1 == iP && iName2 == iT){
	//			std::swap(Prop1, Prop2);
	//		}
	//		return IncompLiquidSI(iOutput,Prop1,Prop2,Ref);
	//	}
	//	else
	//	{
	//		throw ValueError("For incompressible fluids, inputs must be (order does not matter) 'T' and 'P', or 'H' and 'P'");
	//	}
 //   }
 //   // It's an incompressible solution, call the routine
	//else if (IsIncompressibleSolution(Ref))
	//{
	//	if (get_debug_level()>7) std::cout<<__FILE__<<": Identified incompressible solution - "<<Ref.c_str()<<std::endl;
	//	//Enthalpy and pressure are the inputs
	//	if ((iName1 == iH && iName2 == iP) || (iName2 == iH && iName1 == iP))
	//	{
	//		if (iName2 == iH && iName1 == iP)
	//		{
	//			std::swap(Prop1,Prop2);
	//			std::swap(Name1,Name2);
	//		}

	//		// Solve for the temperature
	//		double Tma     = IncompSolutionSI(get_param_index("Tmax"),0.0,0.0,Ref);
	//		double Tmi     = IncompSolutionSI(get_param_index("Tmin"),0.0,0.0,Ref);
	//		double T_guess = (Tma+Tmi)/2.0 ;
	//		double T =_T_hp_secant(Ref,Prop1,Prop2,T_guess);
	//		// Return whatever property is desired
	//		return IncompSolutionSI(iOutput,T,Prop2,Ref);
	//	}
	//	else if ((iName1 == iT && iName2 == iP) || (iName1 == iP && iName2 == iT))
	//	{
	//		if (iName1 == iP && iName2 == iT){
	//			std::swap(Prop1,Prop2);
	//		}
	//		return IncompSolutionSI(iOutput,Prop1,Prop2,Ref);
	//	}
	//	else
	//	{
	//		throw ValueError("For incompressible solutions, inputs must be (order does not matter) 'T' and 'P', or 'H' and 'P'");
	//	}
	//}
	//else
	//{
	//	throw ValueError(format("Your fluid name [%s] is not a CoolProp fluid, a REFPROP fluid, a brine or a liquid",Ref.c_str()));
	//}
//}
//EXPORT_CODE double CONVENTION IProps(long iOutput, long iName1, double Prop1, long iName2, double Prop2, long iFluid)
//{
//    Prop1 = convert_from_unit_system_to_SI(iName1, Prop1, get_standard_unit_system());
//    Prop2 = convert_from_unit_system_to_SI(iName2, Prop2, get_standard_unit_system());
//	double out = IPropsSI(iOutput,iName1,Prop1,iName2,Prop2,iFluid);
//    return convert_from_SI_to_unit_system(iOutput,out,get_standard_unit_system());
//}
//double _CoolProp_Fluid_PropsSI(long iOutput, long iName1, double Prop1, long iName2, double Prop2, Fluid *pFluid)
//{
//	double val = _HUGE, T = _HUGE;
//	// This private method uses the indices directly for speed
//
//	if (get_debug_level()>3){
//		std::cout << format("%s:%d: _CoolProp_Fluid_PropsSI(%d,%d,%g,%d,%g,%s)\n",__FILE__,__LINE__,iOutput,iName1, Prop1, iName2, Prop2, pFluid->get_name().c_str()).c_str();
//	}
//    if (iName1 == iT){ 
//        T = Prop1;} 
//    else if (iName2 == iT){
//        T = Prop2;} 
//
//	// Generate a State instance wrapped around the Fluid instance
//	CoolPropStateClassSI CPS(pFluid);
//
//	// Check if it is an output that doesn't require a state input
//	// Deal with it and return
//	switch (iOutput)
//	{
//        case iI:
//            {
//            if (!ValidNumber(T)){throw ValueError(format("T must be provided as an input to use this output").c_str());}
//            return CPS.pFluid->surface_tension_T(T);
//            }
//        case iRhosatLanc:
//            {
//            if (!ValidNumber(T)){throw ValueError(format("T must be provided as an input to use this output").c_str());}
//            return CPS.pFluid->rhosatL(T);
//            }
//        case iRhosatVanc:
//            {
//            if (!ValidNumber(T)){throw ValueError(format("T must be provided as an input to use this output").c_str());}
//            return CPS.pFluid->rhosatV(T);
//            }
//        case iPsatLanc:
//            {
//            if (!ValidNumber(T)){throw ValueError(format("T must be provided as an input to use this output").c_str());}
//            if (CPS.pFluid->pure()){
//                return CPS.pFluid->psat(T);
//            }
//            else{
//                return CPS.pFluid->psatL(T);
//            }
//            }
//        case iPsatVanc:
//            {
//            if (!ValidNumber(T)){throw ValueError(format("T must be provided as an input to use this output").c_str());}
//            if (CPS.pFluid->pure()){
//                return CPS.pFluid->psat(T);
//            }
//            else{
//                return CPS.pFluid->psatV(T);
//            }
//            }
//		case iMM:
//		case iPcrit:
//		case iTcrit:
//		case iTtriple:
//		case iPtriple:
//        case iPmax:
//        case iTmax:
//		case iRhocrit:
//		case iTmin:
//		case iAccentric:
//		case iPHASE_LIQUID:
//		case iPHASE_GAS:
//		case iPHASE_SUPERCRITICAL:
//		case iPHASE_TWOPHASE:
//		case iGWP20:
//		case iGWP100:
//		case iGWP500:
//		case iODP:
//		case iCritSplineT:
//		case iScrit:
//		case iHcrit:
//		case iTreduce:
//		case iRhoreduce:
//			return CPS.keyed_output(iOutput);
//	}
//
//	// Update the class
//	CPS.update(iName1,Prop1,iName2,Prop2);
//
//	// Debug
//	if (get_debug_level()>9){std::cout << format("%s:%d: State update successful\n",__FILE__,__LINE__).c_str();}
//
//	// Get the output
//	val = CPS.keyed_output(iOutput);
//
//	// Debug
//	if (get_debug_level()>5){std::cout << format("%s:%d: _CoolProp_Fluid_PropsSI returns: %g\n",__FILE__,__LINE__,val).c_str();}
//
//	// Return the value
//	return val;
//}
//EXPORT_CODE double CONVENTION IPropsSI(long iOutput, long iName1, double Prop1, long iName2, double Prop2, long iFluid)
//{
//	pFluid = Fluids.get_fluid(iFluid);
//	// Didn't work
//	if (pFluid == NULL){
//		err_string=std::string("CoolProp error: ").append(format("%d is an invalid fluid index to IProps",iFluid));
//		return _HUGE;
//	}
//	else{
//		// In this function the error catching happens;
//		try{
//			// This is already converted to the right units since we take in SI units
//			return _CoolProp_Fluid_PropsSI(iOutput,iName1,Prop1,iName2,Prop2,pFluid);
//		}
//		catch(std::exception &e){
//			err_string=std::string("CoolProp error: ").append(e.what());
//			return _HUGE;
//		}
//		catch(...){
//			err_string=std::string("CoolProp error: Indeterminate error");
//			return _HUGE;
//		}
//	}
//}


//
//double Props(std::string Output, std::string Name1, double Prop1, std::string Name2, double Prop2, std::string Ref)
//{
//	// Go to the std::string version
//    return PropsS(Output.c_str(),Name1.c_str(),Prop1,Name2.c_str(),Prop2,Ref.c_str());
//}
//double Props(std::string Output,char Name1, double Prop1, char Name2, double Prop2, std::string Ref)
//{
//    return Props(Output, std::string(1,Name1), Prop1, std::string(1,Name2), Prop2, Ref);
//}

///// Calculate some interesting derivatives
//double _CoolProp_Deriv_Terms(long iTerm, double T, double rho, Fluid * pFluid)
//{
//	double val = _HUGE;
//	// This private method uses the indices directly for speed
//
//	if (get_debug_level()>3){
//		std::cout<<__FILE__<<" _CoolProp_Deriv_Terms return: "<<val<<std::endl;
//	}
//
//	switch (iTerm) {
//		case iDERdh_dp__rho:
//		case iDERdh_dp__v:
//		case iDERZ:
//		case iDERdZ_dDelta:
//		case iDERdZ_dTau:
//		case iDERB:
//		case iDERdB_dT:
//		case iDERC:
//		case iDERdC_dT:
//		case iDERphir:
//		case iDERdphir_dTau:
//		case iDERdphir_dDelta:
//		case iDERd2phir_dTau2:
//		case iDERd2phir_dDelta2:
//		case iDERd2phir_dDelta_dTau:
//		case iDERd3phir_dDelta3:
//		case iDERd3phir_dDelta2_dTau:
//		case iDERd3phir_dDelta_dTau2:
//		case iDERd3phir_dTau3:
//		case iDERphi0:
//		case iDERdphi0_dTau:
//		case iDERd2phi0_dTau2:
//		case iDERdphi0_dDelta:
//		case iDERd2phi0_dDelta2:
//		case iDERd2phi0_dDelta_dTau:
//		case iDERd3phi0_dTau3:
//		case iDERdp_dT__rho:
//		case iDERdp_drho__T:
//		case iDERdh_dT__rho:
//		case iDERdh_drho__T:
//		case iDERdrho_dT__p:
//		case iDERdrho_dh__p:
//		case iDERdrho_dp__h:
//			{
//			// Generate a State instance wrapped around the Fluid instance
//			CoolPropStateClass CPS(pFluid);
//
//			// Force the update to consider the inputs as single-phase inputs
//			CPS.flag_SinglePhase =  true;
//
//			// Update the class
//			CPS.update(iT,T,iD,rho);
//			
//			// Get the output value
//			val = CPS.keyed_output(iTerm);
//			break;
//			}
//
//		case iDERrho_smoothed:
//		case iDERdrho_smoothed_dh:
//		case iDERdrho_smoothed_dp:
//		case iDERdrhodh_constp_smoothed:
//		case iDERdrhodp_consth_smoothed:
//		case iDERIsothermalCompressibility:
//			{
//			// Generate a State instance wrapped around the Fluid instance
//			CoolPropStateClass CPS(pFluid);
//
//			// Update the class
//			CPS.update(iT,T,iD,rho);
//
//			// Get the output value
//			val = CPS.keyed_output(iTerm);
//			break;
//			}
//			
//		default:
//			throw ValueError(format("Sorry DerivTerms is a work in progress, your derivative term [%d] is not available!",iTerm));
//	}
//
//	if (get_debug_level()>5){
//		std::cout<<__FILE__<<" _CoolProp_Deriv_Terms return: "<<val<<std::endl;
//	}
//	// Return the value
//	return val;
//}
//
//// Define the functions from the header
//double DerivTerms(long iTerm, double T, double rho, Fluid * pFluid){
//	return _CoolProp_Deriv_Terms(iTerm,T,rho,pFluid);
//}
//double DerivTerms(std::string Term, double T, double rho, std::string Fluidname){
//	if (get_debug_level()>5){
//			std::cout<<__FILE__<<": "<<Term.c_str()<<",T="<<T<<",rho="<<rho<<","<<Fluidname.c_str()<<std::endl;
//		}
//		/*
//	    Derivatives are only supported for CoolProp fluids
//	    */
//	    if (IsCoolPropFluid(Fluidname))
//		{
//			pFluid = Fluids.get_fluid(Fluidname);
//			// for compatibility, replace B and C with VB and VC
//			if ((!Term.compare("B")) || (!Term.compare("C"))) {
//				Term = std::string("V").append(Term);
//			}
//			// Convert all the parameters to integers
//			long iOutput = get_param_index(Term);
//			if (iOutput<0)
//				throw ValueError(format("Your output key [%s] is not valid. (names are case sensitive)",Term.c_str()));
//
//			if (T<=0)
//				throw ValueError(format("Your input temperature [%f] is not valid.",T));
//
//			if (rho<=0)
//				throw ValueError(format("Your input density [%f] is not valid.",rho));
//			// Call the internal method that uses the parameters converted to longs
//			return _CoolProp_Deriv_Terms(iOutput,T,rho,pFluid);
//		}
//		else
//		{
//			throw ValueError(format("Your fluid name [%s] is not a CoolProp fluid.",Fluidname.c_str()));
//		}
//}
//
//int set_reference_stateS(std::string Ref, std::string reference_state)
//{
//	Fluid *pFluid=Fluids.get_fluid(Ref);
//	if (pFluid!=NULL)
//	{
//		return set_reference_stateP(pFluid, reference_state);
//	}
//	else{
//		return -1;
//	}
//}
//
//int set_reference_stateP(Fluid *pFluid, std::string reference_state)
//{
//	CoolPropStateClassSI CPS(pFluid);
//	if (!reference_state.compare("IIR"))
//	{
//		CoolPropStateClassSI CPS(pFluid);
//		CPS.update(iT,273.15,iQ,0);
//		// Get current values for the enthalpy and entropy
//		double h1 = CPS.h();
//		double s1 = CPS.s();
//		double deltah = h1-200000; // offset from 200 kJ/kg enthalpy
//		double deltas = s1-1000; // offset from 1 kJ/kg/K entropy
//		double delta_a1 = deltas/((8314.472/pFluid->params.molemass));
//		double delta_a2 = -deltah/((8314.472/pFluid->params.molemass)*pFluid->reduce.T);
//		pFluid->phi0list.push_back(new phi0_enthalpy_entropy_offset(delta_a1, delta_a2));
//		return 0;
//	}
//	else if (!reference_state.compare("ASHRAE"))
//	{
//		CoolPropStateClassSI CPS(pFluid);
//		CPS.update(iT,233.15,iQ,0);
//		// Get current values for the enthalpy and entropy
//		double h1 = CPS.h();
//		double s1 = CPS.s();
//		double deltah = h1-0; // offset from 0 kJ/kg enthalpy
//		double deltas = s1-0; // offset from 0 kJ/kg/K entropy
//		double delta_a1 = deltas/((8314.472/pFluid->params.molemass));
//		double delta_a2 = -deltah/((8314.472/pFluid->params.molemass)*pFluid->reduce.T);
//		pFluid->phi0list.push_back(new phi0_enthalpy_entropy_offset(delta_a1, delta_a2));
//		return 0;
//	}
//	else if (!reference_state.compare("NBP"))
//	{
//		CoolPropStateClassSI CPS(pFluid);
//		CPS.update(iP,101325.0,iQ,0); // Saturated boiling point at 1 atmosphere
//		// Get current values for the enthalpy and entropy
//		double h1 = CPS.h();
//		double s1 = CPS.s();
//		double deltah = h1-0; // offset from 0 kJ/kg enthalpy
//		double deltas = s1-0; // offset from 0 kJ/kg/K entropy
//		double delta_a1 = deltas/((8314.472/pFluid->params.molemass));
//		double delta_a2 = -deltah/((8314.472/pFluid->params.molemass)*pFluid->reduce.T);
//		pFluid->phi0list.push_back(new phi0_enthalpy_entropy_offset(delta_a1, delta_a2));
//		return 0;
//	}
//	else
//	{ 
//		return -1;
//	}
//
//}
//int set_reference_stateD(std::string Ref, double T, double rho, double h0, double s0)
//{
//	pFluid=Fluids.get_fluid(Ref);
//	if (pFluid!=NULL)
//	{
//		CoolPropStateClassSI CPS(pFluid);
//		CPS.update(iT,T,iD,rho);
//		// Get current values for the enthalpy and entropy
//		double h1 = CPS.h();
//		double s1 = CPS.s();
//		double deltah = h1-h0; // offset from given enthalpy in SI units
//		double deltas = s1-s0; // offset from given enthalpy in SI units
//		double delta_a1 = deltas/((8314.472/pFluid->params.molemass));
//		double delta_a2 = -deltah/((8314.472/pFluid->params.molemass)*pFluid->reduce.T);
//		pFluid->phi0list.push_back(new phi0_enthalpy_entropy_offset(delta_a1, delta_a2));
//		return 0;
//	}
//	else{
//		return -1;
//	}
//}

//std::string get_BibTeXKey(std::string Ref, std::string item)
//{
//	pFluid=Fluids.get_fluid(Ref);
//	if (pFluid!=NULL)
//	{
//		
//		if (!item.compare("EOS")){ return pFluid->BibTeXKeys.EOS; }
//		else if (!item.compare("CP0")){ return pFluid->BibTeXKeys.CP0; }
//		else if (!item.compare("VISCOSITY")){ return pFluid->BibTeXKeys.VISCOSITY; }
//		else if (!item.compare("CONDUCTIVITY")){ return pFluid->BibTeXKeys.CONDUCTIVITY; }
//		else if (!item.compare("ECS_LENNARD_JONES")){ return pFluid->BibTeXKeys.ECS_LENNARD_JONES; }
//		else if (!item.compare("ECS_FITS")){ return pFluid->BibTeXKeys.ECS_FITS; }
//		else if (!item.compare("SURFACE_TENSION")){ return pFluid->BibTeXKeys.SURFACE_TENSION; }
//		else{ return "Bad key";}
//	}
//	else{
//		return std::string("");
//	}
//}
//std::string get_global_param_string(std::string ParamName)
//{
//	if (!ParamName.compare("version"))
//	{
//		return std::string(version);
//	}
//	else if (!ParamName.compare("errstring"))
//	{
//		std::string temp = err_string;
//		err_string = std::string("");
//		return temp;
//	}
//	else if (!ParamName.compare("warnstring"))
//	{
//		std::string temp = warning_string;
//		warning_string = std::string("");
//		return temp;
//	}
//	else if (!ParamName.compare("gitrevision"))
//	{
//		return gitrevision;
//	}
//	else if (!ParamName.compare("FluidsList") || !ParamName.compare("fluids_list"))
//	{
//		return Fluids.FluidList();
//	}
//	else
//	{
//		return format("Input value [%s] is invalid",ParamName.c_str()).c_str();
//	}
//};
//std::string get_fluid_param_string(std::string FluidName, std::string ParamName)
//{
//	try{
//		pFluid = Fluids.get_fluid(FluidName);
//		// Didn't work
//		if (pFluid == NULL){
//			err_string=std::string("CoolProp error: ").append(format("%s is an invalid fluid for get_fluid_param_string",FluidName.c_str()));
//			return format("%s is an invalid fluid for get_fluid_param_string",FluidName.c_str()).c_str();
//		}
//		else{
//			if (!ParamName.compare("aliases"))
//			{
//				std::vector<std::string> v = pFluid->get_aliases();
//				return strjoin(v,", ");
//			}
//			else if (!ParamName.compare("CAS") || !ParamName.compare("CAS_number"))
//			{
//				return pFluid->params.CAS;
//			}
//			else if (!ParamName.compare("ASHRAE34"))
//			{
//				return pFluid->environment.ASHRAE34;
//			}
//			else if (!ParamName.compare("REFPROPName") || !ParamName.compare("REFPROP_name") || !ParamName.compare("REFPROPname"))
//			{
//				return pFluid->get_REFPROPname();
//			}
//			else if (!ParamName.compare("TTSE_mode"))
//			{
//				int mode = pFluid->TTSESinglePhase.get_mode();
//				switch (mode)
//				{
//				case TTSE_MODE_TTSE:
//					return "TTSE";
//				case TTSE_MODE_BICUBIC:
//					return "BICUBIC";
//				default:
//					throw ValueError("TTSE mode is invalid");
//				}
//			}
//			else
//			{
//				return format("Input value [%s] is invalid for Fluid [%s]",ParamName.c_str(),FluidName.c_str()).c_str();
//			}
//		}
//	}
//	catch(std::exception &e)
//	{
//		return(std::string("CoolProp error: ").append(e.what()));
//	}
//	catch(...){
//		return(std::string("CoolProp error: Indeterminate error"));
//	}
//}

} /* namespace CoolProp */