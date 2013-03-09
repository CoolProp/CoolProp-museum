
#cython: embedsignature = True
from __future__ import division
#
#
# This file provides wrapper functions of all the CoolProp functions
#
#
# Each of the functions from the CoolProp header are renamed in cython code to
# an underscored name so that the same name can be used in the exposed functions below
    
#Check for the existence of quantities
cdef bint _quantities_supported
try:
    import quantities as pq
    _quantities_supported = True
except ImportError:
    _quantities_supported = False

import cython
import math

# Default string in Python 3.x is a unicode string (type str)
# Default string in Python 2.x is a byte string(type bytes) 
#
# Create a fused type that allows for either unicode string or bytestring
# We encode unicode strings using the ASCII encoding since we know they are all
# ASCII strings 
ctypedef fused bytes_or_str:
    cython.bytes
    cython.str

include "CyState.pyx"

from param_constants import *
from param_constants_header cimport *

from phase_constants import *
from phase_constants_header cimport *

include "HumidAirProp.pyx"

cdef double _convert_to_desired_units(double value, bytes_or_str parameter_type, bytes_or_str desired_units) except *:
    """
    A convenience function to convert a double value back to the units that
    are desired by the user
    """
    # Convert to a byte-string
    cdef bytes _parameter_type = parameter_type if bytes_or_str is bytes else parameter_type.encode('ascii')
    cdef bytes _desired_units = desired_units if bytes_or_str is bytes else desired_units.encode('ascii')
    # Convert parameter string to index
    cdef long index = _get_param_index(parameter_type)
    # Get the units string by the index
    cdef bytes default_units = _get_index_units(index)
    # Create the Quantity instance
    old = pq.Quantity(value, default_units)
    # Convert the units
    old.units = _desired_units
    #Return the scaled units
    return old.magnitude

cdef _convert_to_default_units(bytes_or_str parameter_type, object parameter):
    """
    A convenience function to convert a quantities instance to the default 
    units required for CoolProp
    """
    cdef bytes _parameter_type = parameter_type if bytes_or_str is bytes else parameter_type.encode('ascii')
    #Convert parameter string to index
    cdef long index = _get_param_index(_parameter_type)
    #Get the units string by the index
    cdef bytes default_units = _get_index_units(index)
    #Rescale the units of the parameter to the default units
    parameter.units = default_units
    #Return the scaled units
    return parameter
    
cpdef double Props(str in1, str in2, in3 = None, in4 = None, in5 = None, in6 = None, in7 = None) except *:
    """
    Call Type #1::

        Props(Fluid,PropName) --> float

    Where ``Fluid`` is a string with a valid CoolProp fluid name, and ``PropName`` is one of the following strings:
    
    =============  ============================
    ``Tcrit``      Critical temperature [K]
    ``pcrit``      Critical pressure [kPa]
    ``rhocrit``    Critical density [kg/m3]
    ``molemass``   Molecular mass [kg/kmol]
    ``Ttriple``    Triple-point temperature [K]
    ``Tmin``       Minimum temperature [K]
    ``ptriple``    Triple-point pressure [kPa]
    ``accentric``  Accentric factor [-]
    =============  ============================
   
    This type of call is used to get fluid-specific parameters that are not 
    dependent on the state 
     
    Call Type #2:
    
    Alternatively, Props can be called in the form::
    
        Props(OutputName,InputName1,InputProp1,InputName2,InputProp2,Fluid) --> float
    
    where ``Fluid`` is a string with a valid CoolProp fluid name.  The value 
    ``OutputName`` is either a single-character or a string alias.  This list 
    shows the possible values
    
    ==========================  ======================================================
    ``OutputName``              Description
    ==========================  ======================================================
    ``Q``                       Quality [-]
    ``T``                       Temperature [K]
    ``P``                       Pressure [kPa]
    ``D``                       Density [kg/m3]
    ``C0``                      Ideal-gas specific heat at constant pressure [kJ/kg]
    ``C``                       Specific heat at constant pressure [kJ/kg]
    ``O``                       Specific heat at constant volume [kJ/kg]
    ``U``                       Internal energy [kJ/kg]
    ``H``                       Enthalpy [kJ/kg]
    ``S``                       Entropy [kJ/kg/K]
    ``A``                       Speed of sound [m/s]
    ``G``                       Gibbs function [kJ/kg]
    ``V``                       Viscosity [Pa-s]
    ``L``                       Thermal conductivity [kW/m/K]
    ``I`` or `SurfaceTension`   Surface Tension [N/m]
    ``w`` or `accentric`        Accentric Factor [-]
    ==========================  ======================================================
    
    The following sets of input values are valid (order doesn't matter):
    
    =========================  ======================================
    ``InputName1``             ``InputName2``
    =========================  ======================================
    ``T``                      ``P``
    ``T``                      ``D``
    ``T``                      ``Q``
    ``P``                      ``Q``
    ``H``                      ``P``
    ``S``                      ``P``
    =========================  ======================================
    
    
    If `InputName1` is `T` and `OutputName` is ``I`` or ``SurfaceTension``, the second input is neglected
    since surface tension is only a function of temperature
    
    Call Type #3:
    New in 2.2
    If you provide InputName1 or InputName2 as a derived class of Quantity, the value will be internally
    converted to the required units as long as it is dimensionally correct.  Otherwise a ValueError will 
    be raised by the conversion
    """
    cdef bytes errs,_in1,_in2,_in3,_in4,_in5,_in6,_in7
    cdef char in2_char, in4_char
        
    if (in4 is None and in6 is None and in7 is None):
        val = _Props1(in1.encode('ascii'), in2.encode('ascii'))
        if math.isinf(val) or math.isnan(val):
            err_string = _get_errstring()
            if not len(err_string) == 0:
                raise ValueError("{err:s} :: inputs were :\"{in1:s}\",\"{in2:s}\"".format(err= err_string,in1=in1,in2=in2))
            else:
                raise ValueError("Props failed ungracefully with inputs:\"{in1:s}\",\"{in2:s}\"; please file a ticket at https://sourceforge.net/p/coolprop/tickets/".format(in1=in1,in2=in2))
        else:
            return val
    else:
        if _quantities_supported:
            if isinstance(in3,pq.Quantity):
                in3 = _convert_to_default_units(in2,in3).magnitude
            if isinstance(in5,pq.Quantity):
                in5 = _convert_to_default_units(<bytes?>in4,in5).magnitude

        in2_char = (<bytes>(in2.encode('ascii')))[0]
        in4_char = (<bytes>(in4.encode('ascii')))[0]
        val = _Props(in1.encode('ascii'), in2_char, in3, in4_char, in5, in6.encode('ascii'))
        if math.isinf(val) or math.isnan(val):
            err_string = _get_errstring()
            if not len(err_string) == 0:
                raise ValueError("{err:s} :: inputs were:\"{in1:s}\",\'{in2:s}\',{in3:0.16e},\'{in4:s}\',{in5:0.16e},\"{in6:s}\"".format(err=err_string,in1=in1,in2=in2,in3=in3,in4=in4,in5=in5,in6=in6))
            else:
                raise ValueError("Props failed ungracefully with inputs:\"{in1:s}\",\'{in2:s}\',{in3:0.16e},\'{in4:s}\',{in5:0.16e},\"{in6:s}\"; please file a ticket at https://sourceforge.net/p/coolprop/tickets/".format(in1=in1,in2=in2,in3=in3,in4=in4,in5=in5,in6=in6))
            
        if not _quantities_supported and in7 is not None:
            raise ValueError("Cannot use output units because quantities package is not installed")
        elif _quantities_supported and in7 is not None: #Then in7 contains a string representation of the units
            #Convert the units to the units given by in7
            return _convert_to_desired_units(val,in1,in7)
        else:
            return val #Error raised by Props2 on failure
    
cpdef double DerivTerms(bytes_or_str Output, double T, double rho, bytes_or_str Fluid):
    """

    .. |cubed| replace:: \ :sup:`3`\ 
    .. |squared| replace:: \ :sup:`2`\ 
    .. |IC| replace:: ``IsothermalCompressibility``
    
    Call signature::
    
        DerivTerms(OutputName, T, rho, Fluid) --> float
    
    where ``Fluid`` is a string with a valid CoolProp fluid name, and ``T`` and ``rho`` are the temperature in K and density in kg/m |cubed| .  The value 
    ``OutputName`` is one of the strings in the table below:
    
    ========================  =====================================================================================================================================
    OutputName                Description
    ========================  =====================================================================================================================================
    ``dpdT``                  Derivative of pressure with respect to temperature at constant density [kPa/K]
    ``dpdrho``                Derivative of pressure with respect to density at constant temperature [kPa/(kg/m\ |cubed|\ )]
    ``Z``                     Compressibility factor [-]
    ``dZ_dDelta``             Derivative of Z with respect to reduced density [-]
    ``dZ_dTau``               Derivative of Z with respect to inverse reduced temperature [-]
    ``B``                     Second virial coefficient [m\ |cubed|\ /kg]
    ``dBdT``                  Derivative of second virial coefficient with respect to temperature [m\ |cubed|\ /kg/K]
    ``C``                     Third virial coefficient [m\ :sup:`6`\ /kg\ |squared|\ ]
    ``dCdT``                  Derivative of third virial coefficient with respect to temperature [m\ :sup:`6`\ /kg\ |squared|\ /K]
    ``phir``                  Residual non-dimensionalized Helmholtz energy [-]
    ``dphir_dTau``            Partial of residual non-dimensionalized Helmholtz energy with respect to inverse reduced temperature [-]
    ``d2phir_dTau2``          Second partial of residual non-dimensionalized Helmholtz energy with respect to inverse reduced temperature [-]
    ``dphir_dDelta``          Partial of residual non-dimensionalized Helmholtz energy with respect to reduced density [-]
    ``d2phir_dDelta2``        Second partial of residual non-dimensionalized Helmholtz energy with respect to reduced density [-]
    ``d2phir_dDelta_dTau``    First cross-partial of residual non-dimensionalized Helmholtz energy [-]
    ``d3phir_dDelta2_dTau``   Second cross-partial of residual non-dimensionalized Helmholtz energy [-]
    ``phi0``                  Ideal-gas non-dimensionalized Helmholtz energy [-]
    ``dphi0_dTau``            Partial of ideal-gas non-dimensionalized Helmholtz energy with respect to inverse reduced temperature [-]
    ``d2phi0_dTau2``          Second partial of ideal-gas non-dimensionalized Helmholtz energy with respect to inverse reduced temperature [-]
    ``dphi0_dDelta``          Partial of ideal-gas non-dimensionalized Helmholtz energy with respect to reduced density [-]
    ``d2phi0_dDelta2``        Second partial of ideal-gas non-dimensionalized Helmholtz energy with respect to reduced density [-]
    |IC|                      Isothermal compressibility [1/kPa]
    ========================  =====================================================================================================================================
    """
    cdef bytes _Fluid = Fluid if bytes_or_str is bytes else Fluid.encode('ascii')
    cdef bytes _Output = Output if bytes_or_str is bytes else Output.encode('ascii')
    return _DerivTerms(_Output,T,rho,_Fluid)

cpdef string Phase_Tp(bytes_or_str Fluid, double T, double p):
    cdef bytes _Fluid = Fluid if bytes_or_str is bytes else Fluid.encode('ascii')
    return _Phase_Tp(_Fluid,T,p)
    
cpdef string Phase_Trho(bytes_or_str Fluid, double T, double rho):
    cdef bytes _Fluid = Fluid if bytes_or_str is bytes else Fluid.encode('ascii')
    return _Phase_Trho(_Fluid,T,rho)
    
cpdef string Phase(bytes_or_str Fluid, double T, double p):
    """
    Given a set of temperature and pressure, returns one of the following strings

    * Gas
    * Liquid
    * Supercritical
    * Two-Phase
    
    Phase diagram::
    
            |         |     
            |         |    Supercritical
            |         |
        p   | Liquid (b)------------
            |        /
            |       / 
            |      /       Gas
            |     / 
            |   (a)
            |  /
            |------------------------
    
                       T
    
           a: triple point
           b: critical point
           a-b: Saturation line
    """
    cdef bytes _Fluid = Fluid if bytes_or_str is bytes else Fluid.encode('ascii')
    return _Phase(_Fluid,T,p)

cpdef F2K(double T_F):
    """
    Convert temperature in degrees Fahrenheit to Kelvin
    """
    return _F2K(T_F)

cpdef K2F(double T_K):
    """
    Convert temperature in Kelvin to degrees Fahrenheit
    """
    return _K2F(T_K)

cpdef list FluidsList():
    """
    Return a list of strings of all fluid names
    
    Returns
    -------
    FluidsList : list of strings of fluid names
        All the fluids that are included in CoolProp
    
    Notes
    -----
    
    Here is an example::
        
       In [0]: from CoolProp.CoolProp import FluidsList
    
       In [1]: FluidsList()
       
    """ 
    cdef string FL = _FluidsList()
    return [F.encode('ascii') for F in FL.decode('ascii').split(',')]

cpdef get_aliases(bytes_or_str Fluid):
    """
    Return a comma separated string of aliases for the given fluid
    """
    cdef bytes _Fluid = Fluid if bytes_or_str is bytes else Fluid.encode('ascii')
    return [F.encode('ascii') for F in (<str>_get_aliases(_Fluid)).decode('ascii').split(',')]
    
cpdef string get_REFPROPname(bytes_or_str Fluid):
    """
    Return the REFPROP compatible name for the fluid (only useful on windows)
    
    Some fluids do not use the REFPROP name.  For instance, 
    ammonia is R717, and propane is R290.  You can still can still call CoolProp
    using the name ammonia or R717, but REFPROP requires that you use a limited
    subset of names.  Therefore, this function that returns the REFPROP compatible
    name.  To then use this to call REFPROP, you would do something like::
    
       In [0]: from CoolProp.CoolProp import get_REFPROPname, Props
    
       In [1]: Fluid = 'REFPROP-' + get_REFPROPname('R290')
       
       In [2]: Props('D', 'T', 300, 'P', 300, Fluid)
    """
    cdef bytes _Fluid = Fluid if bytes_or_str is bytes else Fluid.encode('ascii')
    return _get_REFPROPname(_Fluid)

cpdef string get_errstr():
    """
    Return the current error string
    """
    return _get_errstring()

cpdef set_debug(int level):
    """
    Set the current debug level as integer in the range [0,10]
    
    Parameters
    ----------
    level : int
        If level is 0, no output will be written to screen, if >0, 
        some output will be written to screen.  The larger level is, 
        the more verbose the output will be
    """
    _debug(level)

cpdef get_debug():
    """
    Return the current debug level as integer
    """
    return _get_debug()
    
cpdef string get_EOSReference(bytes_or_str Fluid):
    """
    Return a string with the reference for the equation of state
    """
    cdef bytes _Fluid = Fluid if bytes_or_str is bytes else Fluid.encode('ascii')
    return _get_EOSReference(_Fluid)

cpdef string get_TransportReference(bytes_or_str Fluid):
    """
    Return a string with the reference for the transport properties (thermal conductivity, viscosity, surface tension)
    """
    cdef bytes _Fluid = Fluid if bytes_or_str is bytes else Fluid.encode('ascii')
    return _get_TransportReference(_Fluid)
    
cpdef bint IsFluidType(bytes_or_str Fluid, bytes_or_str Type):
    """
    Check if a fluid is of a given type
    
    Valid types are:
    - Brine
    - PseudoPure (or equivalently PseudoPureFluid)
    - PureFluid
    """
    cdef bytes _Fluid = Fluid if bytes_or_str is bytes else Fluid.encode('ascii')
    cdef bytes _Type = Type if bytes_or_str is bytes else Type.encode('ascii')
    if _IsFluidType(_Fluid, _Type):
        return True
    else:
        return False
    
#: Enable the TTSE
cpdef bint enable_TTSE_LUT(char *FluidName): return _enable_TTSE_LUT(FluidName)
#: Check if TTSE is enabled
cpdef bint isenabled_TTSE_LUT(char *FluidName): return _isenabled_TTSE_LUT(FluidName)
#: Disable the TTSE
cpdef bint disable_TTSE_LUT(char *FluidName): return _disable_TTSE_LUT(FluidName)

#: Enable the TTSE
cpdef bint enable_TTSE_LUT_writing(char *FluidName): return _enable_TTSE_LUT_writing(FluidName)
#: Check if TTSE is enabled
cpdef bint isenabled_TTSE_LUT_writing(char *FluidName): return _isenabled_TTSE_LUT_writing(FluidName)
#: Disable the TTSE
cpdef bint disable_TTSE_LUT_writing(char *FluidName): return _disable_TTSE_LUT_writing(FluidName)

#: Over-ride the default size of both of the saturation LUT
cpdef bint set_TTSESat_LUT_size(char *FluidName, int N): return _set_TTSESat_LUT_size(FluidName, N)
#: Over-ride the default size of the single-phase LUT
cpdef bint set_TTSESinglePhase_LUT_size(char *FluidName, int Np, int Nh): return _set_TTSESinglePhase_LUT_size(FluidName, Np, Nh)
#: Over-ride the default range of the single-phase LUT
cpdef bint set_TTSESinglePhase_LUT_range(char *FluidName, double hmin, double hmax, double pmin, double pmax): return _set_TTSESinglePhase_LUT_range(FluidName, hmin, hmax, pmin, pmax)

cpdef tuple get_TTSESinglePhase_LUT_range(char *FluidName):
    """
    Get the current range of the single-phase LUT
    
    Returns
    -------
    tuple of hmin,hmax,pmin,pmax
    """
    cdef double hmin = 0, hmax = 0, pmin = 0, pmax = 0
    #In cython, hmin[0] to get pointer rather than &hmin
    cdef bint valsok = _get_TTSESinglePhase_LUT_range(FluidName, &hmin, &hmax, &pmin, &pmax)
    if valsok:
        return (hmin, hmax, pmin, pmax)
    else:
        raise ValueError("Either your FluidName was invalid or LUT bounds not available since no call has been made to tables")
    
cpdef viscosity_dilute(bytes_or_str Fluid, double T, double rho, double e_k, double sigma):
    cdef bytes _Fluid = Fluid if bytes_or_str is bytes else Fluid.encode('ascii')
    return _viscosity_dilute(_Fluid,T,rho,e_k,sigma)

cpdef tuple conformal_Trho(bytes_or_str Fluid, bytes_or_str ReferenceFluid, double T, double rho):
    cdef double T0 = 0,rho0= 0
    cdef bytes _Fluid = Fluid if bytes_or_str is bytes else Fluid.encode('ascii')
    cdef bytes _ReferenceFluid = ReferenceFluid if bytes_or_str is bytes else ReferenceFluid.encode('ascii')
    _conformal_Trho(_Fluid,_ReferenceFluid,T,rho,&T0,&rho0)
    return T0,rho0

cpdef rhosatL_anc(bytes_or_str Fluid, double T):
    cdef bytes _Fluid = Fluid if bytes_or_str is bytes else Fluid.encode('ascii')
    return _rhosatL_anc(_Fluid,T)

cpdef rhosatV_anc(bytes_or_str Fluid, double T):
    cdef bytes _Fluid = Fluid if bytes_or_str is bytes else Fluid.encode('ascii')
    return _rhosatV_anc(_Fluid,T)

cpdef psatL_anc(bytes_or_str Fluid, double T):
    cdef bytes _Fluid = Fluid if bytes_or_str is bytes else Fluid.encode('ascii')
    return _psatL_anc(_Fluid,T)

cpdef psatV_anc(bytes_or_str Fluid, double T):
    cdef bytes _Fluid = Fluid if bytes_or_str is bytes else Fluid.encode('ascii')
    return _psatV_anc(_Fluid,T)

cpdef int debug(int level):
    """
    Sets the debug level
    
    Parameters
    ----------
    level : int
        Flag indicating how verbose the debugging should be.
            0 : no debugging output
            ...
            ...
            10 : very annoying debugging output - every function call debugged
    """
    _debug(level)
        
from math import pow as pow_

#A dictionary mapping parameter index to string for use with non-CoolProp fluids
cdef dict paras = {iD : 'D',
                   iQ : 'Q',
                   iMM : 'M',
                   iT : 'T',
                   iH : 'H',
                   iP : 'P',
                   iC : 'C',
                   iC0 : 'C0',
                   iO : 'O',
                   iV : 'V',
                   iL : 'L',
                   iS : 'S',
                   iU : 'U',
                   iDpdT : 'dpdT'}

cdef dict paras_inverse = {v:k for k,v in paras.iteritems()}

cdef class State:
    """
    A class that contains all the code that represents a thermodynamic state
    """
        
    def __init__(self, bytes Fluid, dict StateDict, double xL=-1.0, object Liquid = None, object phase = None):
        """
        Parameters
        ----------
        Fluid, string
        StateDict, dictionary
            The state of the fluid - passed to the update function
        phase, string
            The phase of the fluid, it it is known.  One of 'Gas','Liquid','Supercritical','TwoPhase'
        xL, float
            Liquid mass fraction (not currently supported)
        Liquid, string
            The name of the liquid (not currently supported)
        
            
        """
        cdef bytes _Fluid = Fluid
        cdef bytes _Liquid
        
        if Liquid is None:
            _Liquid = b''
        elif isinstance(Liquid,str):
            _Liquid = Liquid.encode('ascii')
        elif isinstance(Liquid,bytes):
            _Liquid = Liquid
        else:
            raise TypeError()
            
        if phase is None:
            _phase = b''
        elif isinstance(phase,str):
            _phase = phase.encode('ascii')
        elif isinstance(phase,bytes):
            _phase = phase
        else:
            raise TypeError()
        
        self.Fluid = _Fluid
        self.iFluid = _get_Fluid_index(_Fluid)
        #Try to get the fluid from CoolProp
        if self.iFluid >= 0:
            #It is a CoolProp Fluid so we can use the faster integer passing function
            self.is_CPFluid = True
            self.PFC = PureFluidClass(Fluid)
        else:
            self.is_CPFluid = False
        self.xL = xL
        self.Liquid = _Liquid
        self.phase = _phase
        #Parse the inputs provided
        self.update(StateDict)
        #Set the phase flag
        if self.phase == str('Gas') or self.phase == str('Liquid') or self.phase == str('Supercritical'):
            if self.is_CPFluid and (self.phase == str('Gas') or self.phase == str('Liquid')):
                self.PFC.CPS.flag_SinglePhase = True
            elif not self.is_CPFluid and self.phase is not None:
                _set_phase(self.phase)
            
    def __reduce__(self):
        d={}
        d['xL']=self.xL
        d['Liquid']=self.Liquid
        d['Fluid']=self.Fluid
        d['T']=self.T_
        d['rho']=self.rho_
        d['phase'] = self.phase
        return rebuildState,(d,)
          
    cpdef update_ph(self, double p, double h):
        """
        Use the pressure and enthalpy directly
        
        Parameters
        ----------
        p: float
            Pressure (absolute) [kPa]
        h: float
            Enthalpy [kJ/kg]
        
        """
        self.p_ = p
        cdef double T
        
        if self.is_CPFluid:
            self.PFC.update(iP, p, iH, h)
            self.T_ = self.PFC.T()
            self.rho_ = self.PFC.rho()
        else:
            T = _Props('T','P',p,'H',h,self.Fluid)
            if abs(T)<1e90:
                self.T_=T
            else:
                errstr = _get_errstring()
                raise ValueError(errstr)
            self.rho_ = _Props('D','P',p,'H',h,self.Fluid)
            
    cpdef update_Trho(self, double T, double rho):
        """
        Just use the temperature and density directly for speed
        
        Parameters
        ----------
        T: float
            Temperature [K]
        rho: float
            Density [kg/m^3]
            
        """
        cdef double p
        self.T_ = T
        self.rho_ = rho
        
        if self.is_CPFluid:
            self.PFC.update(iT,T,iD,rho)
            p = self.PFC.p()
        else:
            p = _Props('P','T',T,'D',rho,self.Fluid)
        
        if abs(p)<1e90:
            self.p_=p
        else:
            errstr = _get_errstring()
            raise ValueError(errstr)
        
    cpdef update(self, dict params, double xL=-1.0):
        """
        Parameters
        params, dictionary 
            A dictionary of terms to be updated, with keys equal to single-char inputs to the Props function,
            for instance ``dict(T=298, P = 101.325)`` would be one standard atmosphere
        """
            
        cdef double p
        cdef long iInput1, iInput2
        cdef bytes errstr
        
        # If no value for xL is provided, it will have a value of -1 which is 
        # impossible, so don't update xL
        if xL > 0:
            #There is liquid
            self.xL=xL
            self.hasLiquid=True
        else:
            #There's no liquid
            self.xL=0.0
            self.hasLiquid=False
        
        #Given temperature and pressure, determine density of gas 
        # (or gas and oil if xL is provided)
        if abs(self.xL)<=1e-15:
            
            if self.is_CPFluid:
                items = params.items()
                iInput1 = paras_inverse[items[0][0]]
                iInput2 = paras_inverse[items[1][0]]
                try:
                    self.PFC.update(iInput1, items[0][1], iInput2, items[1][1])
                except:
                    raise
                self.T_ = self.PFC.T()
                self.p_ = self.PFC.p()
                self.rho_ = self.PFC.rho()
                return
            
            #Get the density if T,P provided, or pressure if T,rho provided
            if 'P' in params:
                self.p_=params['P']
                rho = _Props('D','T',self.T_,'P',self.p_,self.Fluid)
                
                if abs(rho) < 1e90:
                    self.rho_=rho
                else:
                    errstr = _get_errstring()
                    raise ValueError(errstr)
            elif 'D' in params:
                self.rho_=params['D']
                p = _Props('P','T',self.T_,'D',self.rho_,self.Fluid)
                
                if abs(p)<1e90:
                    self.p_=p
                else:
                    errstr = _get_errstring()
                    raise ValueError(errstr+str(params))
            elif 'Q' in params:
                p = _Props('P','T',self.T_,'Q',params['Q'],self.Fluid)
                self.rho_ = _Props('D','T',self.T_,'Q',params['Q'],self.Fluid)
                
                if abs(self.rho_)<1e90:
                    pass
                else:
                    errstr = _get_errstring()
                    raise ValueError(errstr+str(params))
            else:
                raise KeyError("Dictionary must contain the key 'T' and one of 'P' or 'D'")
            
        elif self.xL>0 and self.xL<=1:
            raise ValueError('xL is out of range - value for xL is [0,1]')
        else:
            raise ValueError('xL must be between 0 and 1')
        
    cpdef long Phase(self) except *:
        """
        Returns an integer flag for the phase of the fluid, where the flag value
        is one of iLiquid, iSupercritical, iGas, iTwoPhase
        
        These constants are defined in the phase_constants module, and are imported
        into this module
        """
        
        if self.is_CPFluid:
            return self.PFC.phase()
        else:
            raise NotImplementedError("Phase not defined for fluids other than CoolProp fluids")
        
    cpdef double Props(self, long iOutput) except *: 
        if iOutput<0:
            raise ValueError('Your output is invalid') 
        
        if self.is_CPFluid:
            return self.PFC.keyed_output(iOutput)
        else:
            return _Props(paras[iOutput],'T',self.T_,'D',self.rho_,self.Fluid)
            
    cpdef double get_Q(self) except *:
        """ Get the quality [-] """
        return self.Props(iQ)
    property Q:
        """ The quality [-] """
        def __get__(self):
            return self.get_Q()
    
    cpdef double get_MM(self) except *:
        """ Get the mole mass [kg/kmol] or [g/mol] """
        return _Props1(self.Fluid,'molemass')
    
    cpdef double get_rho(self) except *:
        """ Get the density [kg/m^3] """ 
        return self.rho_
    property rho:
        """ The density [kg/m^3] """
        def __get__(self):
            return self.rho_
            
    cpdef double get_p(self) except *:
        """ Get the pressure [kPa] """ 
        return self.p_
    property p:
        """ The pressure [K] """
        def __get__(self):
            return self.p_
    
    cpdef double get_T(self) except *: 
        """ Get the temperature [K] """
        return self.T_
    property T:
        """ The temperature [K] """
        def __get__(self):
            return self.T_
    
    cpdef double get_h(self) except *: 
        """ Get the specific enthalpy [kJ/kg] """
        return self.Props(iH)
    property h:
        """ The specific enthalpy [kJ/kg] """
        def __get__(self):
            return self.get_h()
          
    cpdef double get_u(self) except *: 
        """ Get the specific internal energy [kJ/kg] """
        return self.Props(iU)
    property u:
        """ The internal energy [kJ/kg] """
        def __get__(self):
            return self.get_u()
            
    cpdef double get_s(self) except *: 
        """ Get the specific enthalpy [kJ/kg/K] """
        return self.Props(iS)
    property s:
        """ The specific enthalpy [kJ/kg/K] """
        def __get__(self):
            return self.get_s()
    
    cpdef double get_cp0(self) except *:
        """ Get the specific heat at constant pressure for the ideal gas [kJ/kg] """
        return self.Props(iC0)
    
    cpdef double get_cp(self) except *: 
        """ Get the specific heat at constant pressure  [kJ/kg] """   
        return self.Props(iC)
    property cp:
        """ The specific heat at constant pressure  [kJ/kg] """
        def __get__(self):
            return self.get_cp()
            
    cpdef double get_cv(self) except *: 
        """ Get the specific heat at constant volume  [kJ/kg] """
        return self.Props(iO)
    property cv:
        """ The specific heat at constant volume  [kJ/kg] """
        def __get__(self):
            return self.get_cv()
            
    cpdef double get_visc(self) except *:
        """ Get the viscosity, in [Pa-s]"""
        return self.Props(iV)
    property visc:
        """ The viscosity, in [Pa-s]"""
        def __get__(self):
            return self.get_visc()

    cpdef double get_cond(self) except *:
        """ Get the thermal conductivity, in [kW/m/K]"""
        return self.Props(iL)
    property k:
        """ The thermal conductivity, in [kW/m/K]"""
        def __get__(self):
            return self.get_cond()
            
    property Prandtl:
        """ The Prandtl number (cp*mu/k) [-] """
        def __get__(self):
            return self.cp * self.visc / self.k
            
    cpdef double get_dpdT(self) except *:
        return self.Props(iDpdT)
    property dpdT:
        def __get__(self):
            return self.get_dpdT()
        
    cpdef speed_test(self, int N):
        from time import clock
        cdef int i
        cdef char * k
        cdef string Fluid = self.Fluid
        cdef long IT = 'T'
        cdef long ID = 'D'
        import CoolProp as CP
        
        print 'Call to the Python call layer'
        print "'M' involves basically no computational effort and is a good measure of the function call overhead"
        keys = ['H','P','S','U','C','O','V','L','M','C0','dpdT']
        for key in keys:
            t1=clock()
            for i in range(N):
                CP.Props(key,'T',self.T_,'D',self.rho_,Fluid)
            t2=clock()
            print 'Elapsed time for {0:d} calls for "{1:s}" at {2:g} us/call'.format(N,key,(t2-t1)/N*1e6)
            
        print 'Direct c++ call to CoolProp without the Python call layer'
        print "'M' involves basically no computational effort and is a good measure of the function call overhead"
        keys = ['H','P','S','U','C','O','V','L','M','C0','dpdT']
        for key in keys:
            t1=clock()
            for i in range(N):
                _Props(key,'T',self.T_,'D',self.rho_,Fluid)
            t2=clock()
            print 'Elapsed time for {0:d} calls for "{1:s}" at {2:g} us/call'.format(N,key,(t2-t1)/N*1e6)
        
        print 'Call to the c++ layer using integers'
        keys = [iH,iP,iS,iU,iC,iO,iV,iL,iMM,iC0,iDpdT]
        for key in keys:
            t1=clock()
            for i in range(N):
                _IProps(key,iT,self.T_,iD,self.rho_,self.iFluid)
            t2=clock()
            print 'Elapsed time for {0:d} calls for "{1:s}" at {2:g} us/call'.format(N,paras[key],(t2-t1)/N*1e6)
            
        print 'Call to the c++ layer using integers'
        keys = [iH,iP,iS,iU,iC,iO,iV,iL,iMM,iC0,iDpdT]
        for key in keys:
            t1=clock()
            self.PFC.update(iT,self.T_,iD,self.rho_)
            for i in range(N):
                self.PFC.keyed_output(key)
            t2=clock()
            print 'Elapsed time for {0:d} calls for "{1:s}" at {2:g} us/call'.format(N,paras[key],(t2-t1)/N*1e6)
        
        print 'Using CoolPropStateClass with T,rho'
        keys = [iH,iP,iC,iO,iDpdT]
        t1=clock()
        for i in range(N):
            self.PFC.update(iT,self.T_,iD,self.rho_)
            for key in keys:
                self.PFC.keyed_output(key)
        t2=clock()
        print 'Elapsed time for {0:d} calls of iH,iP,iC,iO,iDpdT takes {1:g} us/call'.format(N,(t2-t1)/N*1e6)
            
        
        keys = [iH,iP,iS,iU,iC,iO,iV,iL,iMM,iC0,iDpdT]
        isenabled = _isenabled_TTSE_LUT(<bytes>Fluid)
        _enable_TTSE_LUT(<bytes>Fluid)
        _IProps(iH,iT,self.T_,iD,self.rho_,self.iFluid)
        
        print 'Call using TTSE with T,rho'
        print "'M' involves basically no computational effort and is a good measure of the function call overhead"
        for key in keys:
            t1=clock()
            self.PFC.update(iT,self.T_,iD,self.rho_)
            for i in range(N):
                self.PFC.keyed_output(key)
            t2=clock()
            print 'Elapsed time for {0:d} calls for "{1:s}" at {2:g} us/call'.format(N,paras[key],(t2-t1)/N*1e6)
            
        print 'Call using TTSE with p,h'
        print "'M' involves basically no computational effort and is a good measure of the function call overhead"
        cdef double hh = self.h
        for key in keys:
            t1=clock()
            self.PFC.update(iP,self.p_,iH,hh)
            for i in range(N):
                self.PFC.keyed_output(key)
            t2=clock()
            print 'Elapsed time for {0:d} calls for "{1:s}" at {2:g} us/call'.format(N,paras[key],(t2-t1)/N*1e6)
        
        print 'Using CoolPropStateClass with T,rho with LUT'
        keys = [iH,iP,iC,iO,iDpdT]
        t1=clock()
        for i in range(N):
            self.PFC.update(iT,self.T_,iD,self.rho_)
            for key in keys:
                self.PFC.keyed_output(key)
        t2=clock()
        print 'Elapsed time for {0:d} calls of iH,iP,iC,iO,iDpdT takes {1:g} us/call'.format(N,(t2-t1)/N*1e6)
        
        if not isenabled:
            _disable_TTSE_LUT(<bytes>Fluid)
    
    def __str__(self):
        """
        Return a string representation of the state
        """
        units={'T': 'K', 
               'p': 'kPa', 
               'rho': 'kg/m^3',
               'Q':'kg/kg',
               'h':'kJ/kg',
               'u':'kJ/kg',
               's':'kJ/kg/K',
               'visc':'Pa-s',
               'k':'kW/m/K',
               'cp':'kJ/kg/K',
               'cv':'kJ/kg/K',
               'dpdT':'kPa/K'}
        s='phase = '+self.phase+'\n'
        for k in ['T','p','rho','Q','h','u','s','visc','k','cp','cv','dpdT','Prandtl']:
            if k in units:
                s+=k+' = '+str(getattr(self,k))+' '+units[k]+'\n'
            else:
                s+=k+' = '+str(getattr(self,k))+' NO UNITS'+'\n'
        return s.rstrip()
        
    cpdef copy(self):
        cdef double T = self.T_
        cdef double rho = self.rho_
        ST = State(self.Fluid,{'T':T,'D':rho},phase = self.phase)
        return ST
    
def rebuildState(d):
    S=State(d['Fluid'],{'T':d['T'],'D':d['rho']},phase=d['phase'])
    S.xL = d['xL']
    S.Liquid=d['Liquid']
    return S
    
