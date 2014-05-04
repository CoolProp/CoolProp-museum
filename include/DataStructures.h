/*
 * DataStructures.h
 *
 *  Created on: 21 Dec 2013
 *      Author: jowr
 */

#ifndef DATASTRUCTURES_H_
#define DATASTRUCTURES_H_

#include <map>
namespace CoolProp {

struct SimpleState
{
    double rhomolar, T, p, hmolar, smolar, umolar;
};

/// --------------------------------------------------
/// Define some constants that will be used throughout
/// --------------------------------------------------
/// These are constants for the input and output parameters
/// The structure is taken directly from the AbstractState class.
enum parameters{
    // Bulk properties
    iT, iDmolar, iP, iQ, iHmolar, iSmolar, iCpmolar, iCvmolar, iUmolar, ispeed_sound, iisothermal_compressibility, iisobaric_expansion_coefficient,

    // Smoothing functions for density
    idrhodh_constp_smoothed, idrhodp_consth_smoothed, irho_smoothed,

    // Transport properties
    iviscosity, iconductivity, isurface_tension,

    // Derivatives of properties
    idvdp_constT, idvdT_constp,
    // Density
    idrhodh_constp, idrhodp_consth, idrhodp_constT, idrhodT_constp, id2rhodh2_constp,
    id2rhodhdp, id2rhodhdQ, id2rhodp2_constT, id2rhodpdQ, id2rhodT2_constp, id2rhodTdp,
    // Pressure
    idpdrho_consth, idpdrho_constT, idpdT_consth, idpdT_constrho, id2pdrho2_constT,
    id2pdrhodT, id2pdT2_constrho,
    // Enthalpy
    idhdp_constrho, idhdp_constT, idhdrho_constp, idhdrho_constT, idhdT_constp,
    idhdT_constrho, id2hdp2_constT, id2hdrho2_constT, id2hdrhodT, id2hdT2_constp,
    id2hdT2_constrho, id2hdTdp,
    // Entropy
    idsdp_constT, idsdrho_constp, idsdrho_constT, idsdT_constp, idsdT_constrho,
    id2sdp2_constT,	id2sdrho2_constT, id2sdrhodT, id2sdT2_constp, id2sdT2_constrho,
    id2sdTdp,

    // Fundamental derivative of gas dynamics
    ifundamental_derivative_of_gas_dynamics, id2pdv2_consts,

    // Other functions and derivatives
    iB, iC, iZ, idBdT, idCdT, idZdDelta, idZdTau
};

/// These are constants for the phases of the fluid
enum phases {iphase_liquid, iphase_supercritical, iphase_gas, iphase_twophase, iphase_unknown};

/// These are unit types for the fluid
enum fluid_types{FLUID_TYPE_PURE, FLUID_TYPE_PSEUDOPURE, FLUID_TYPE_REFPROP, FLUID_TYPE_INCOMPRESSIBLE_LIQUID, FLUID_TYPE_INCOMPRESSIBLE_SOLUTION, FLUID_TYPE_UNDEFINED};

// !! If you add a parameter, update the map in the corresponding CPP file !!
/// These are input pairs that can be used (in each pair, input keys are sorted alphabetically)
enum input_pairs{
    QT_INPUTS, ///< Molar quality, Temperature in K
    PQ_INPUTS, ///< Pressure in Pa, Molar quality
                 
    PT_INPUTS, ///< Pressure in Pa, Temperature in K

    DmassT_INPUTS, ///< Mass density in kg/m^3, Temperature in K
    DmolarT_INPUTS, ///< Molar density in mol/m^3, Temperature in K
    HmolarT_INPUTS, ///< Enthalpy in J/mol, Temperature in K
    HmassT_INPUTS, ///< Enthalpy in J/kg, Temperature in K
    SmolarT_INPUTS, ///< Entropy in J/mol/K, Temperature in K
    SmassT_INPUTS, ///< Entropy in J/kg/K, Temperature in K
    TUmolar_INPUTS, ///< Temperature in K, Internal energy in J/mol
    TUmass_INPUTS, ///< Temperature in K, Internal energy in J/kg
                 
    DmassP_INPUTS, ///< Mass density in kg/m^3, Pressure in Pa
    DmolarP_INPUTS, ///< Molar density in mol/m^3, Pressure in Pa
    HmassP_INPUTS, ///< Enthalpy in J/kg, Pressure in Pa
    HmolarP_INPUTS, ///< Enthalpy in J/mol, Pressure in Pa
    PSmass_INPUTS, ///< Pressure in Pa, Entropy in J/kg/K
    PSmolar_INPUTS, ///< Pressure in Pa, Entropy in J/mol/K 
    PUmass_INPUTS, ///< Pressure in Pa, Internal energy in J/kg
    PUmolar_INPUTS, ///< Pressure in Pa, Internal energy in J/mol
                 
    HmassSmass_INPUTS, ///< Enthalpy in J/kg, Entropy in J/kg/K
    HmolarSmolar_INPUTS, ///< Enthalpy in J/mol, Entropy in J/mol/K
    SmassUmass_INPUTS, ///< Entropy in J/kg/K, Internal energy in J/kg
    SmolarUmolar_INPUTS, ///< Entropy in J/mol/K, Internal energy in J/mol
                 
    DmassHmass_INPUTS, ///< Mass density in kg/m^3, Enthalpy in J/kg
    DmolarHmolar_INPUTS, ///< Molar density in mol/m^3, Enthalpy in J/mol
    DmassSmass_INPUTS, ///< Mass density in kg/m^3, Entropy in J/kg/K
    DmolarSmolar_INPUTS, ///< Molar density in mol/m^3, Entropy in J/mol/K
    DmassUmass_INPUTS, ///< Mass density in kg/m^3, Internal energy in J/kg
    DmolarUmolar_INPUTS, ///< Molar density in mol/m^3, Internal energy in J/mol
};
// !! If you add or remove a parameter, update the map in the corresponding CPP file !!

/// Return the short description of an input pair key ("DmolarT_INPUTS" for instance)
std::string get_input_pair_short_desc(int pair);
/// Return the long description of an input pair key ("Molar density in mol/m^3, Temperature in K" for instance)
std::string get_input_pair_long_desc(int pair);

} /* namespace CoolProp */
#endif /* DATASTRUCTURES_H_ */
