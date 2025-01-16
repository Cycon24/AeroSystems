# -*- coding: utf-8 -*-
"""
Created on Tue Nov 14 10:25:46 2023

@author: cycon
"""
import EngineErrors as EE


def mdot_1(V, A, rho=None, T=None, P=None, R=287):
    '''
    Calculates the mass flow rate through a single area/nozzle

    Parameters
    ----------
    V : Float [m/s]
        Mean air velocity through Area A.
    A : Float [m^2]
        Area at point of interest.
    rho : Float [kg/m^3], optional
        Static Air Density, if not included T and P is necessary. The default is None.
    T : Float [K], optional
        Static Air Temperatire. Not needed if rho is included. The default is None.
    P : Float [Pa], optional
        Static Air Pressure. Not needed if rho is included.. The default is None.
    R : Float/Int [J/kg*K], optional
        The gas constant of the fluid (assuming air). Not needed of rho is included. The default is 287.

    Returns
    -------
    mdot : Float
        Mass flow rate in kg/s through area A.

    '''

    if rho == None:
        if T == None or P == None:
            raise EE.IncompleteInputs(
                'Rho or Pressure and Temperature', 'mdot')
        else:
            rho = P/(R*T)

    mdot = rho*V*A
    return mdot


def mdot_2(F, BPR, C9, C19, Ca, f=0, P9=None, P19=None, Pa=None, A9=None, A19=None):
    '''

    Parameters
    ----------
    Calculates the mass flow rate through an engine with a bypass/two nozzles based on the known thrust value
    
    Parameters
    ----------
    F : Float [N]
        Thrust of the engine.
    BPR : Float, >= 1
        Bypass ratio of the engine.
    C9 : Float [m/s]
        Velocity of the core exhaust nozzle (rel to engine).
    C19 : Float [m/s]
        Velocty of the bypass exhaust nozzle (rel to engine).
    Ca : Float [m/s]
        Airspeed of the inlet flow.
    f : Float, optional
        Fuel to air ratio, when inputted the airmass flow is adjusted to account for the momentum added by fuel mass flow. The default is 0.
    P9 : Float [Pa], optional
        Static exit pressure of core nozzle. The default is None.
    P19 : Float [Pa], optional
        Static exit pressure of bypass nozzle. The default is None.
    Pa : Float [Pa], optional
        Static atmospheric pressure. The default is None.
    A9 : Float [m^2], optional
        Area at the exit of the core nozzle. The default is None.
    A19 : Float [m^2], optional
        Area at the exit of the bypass nozzle. The default is None.
    
    Returns
    -------
    mdot: Float [kg/s]
        Mass flow rate through the entire engine.

    '''
    press_thrust = 0 
    if A9 != None and P9 != None and Pa != None:
        press_thrust += A9*(P9-Pa)
    if A19 != None and P19 != None and Pa != None:
        press_thrust += A19*(P19-Pa) 
        
    den = ((1+f)/(BPR+1))*(C9-Ca) + (BPR/(BPR+1))*(C19-Ca)
    mdot = (F-press_thrust)/den
    return


def Thrust_1(mdot, BPR, C9, C19, Ca, f=0, P9=None, P19=None, Pa=None, A9=None, A19=None):
    '''
    Calculates the mass flow rate through an engine with a bypass/two nozzles based on the known thrust value

    Parameters
    ---------
    mdot: Float [kg/s]
        Mass flow rate through the entire engine.
    BPR : Float, >= 1
        Bypass ratio of the engine.
    C9 : Float [m/s]
        Velocity of the core exhaust nozzle (rel to engine).
    C19 : Float [m/s]
        Velocty of the bypass exhaust nozzle (rel to engine).
    Ca : Float [m/s]
        Airspeed of the inlet flow.
    f : Float, optional
        Fuel to air ratio, when inputted the airmass flow is adjusted to account for the momentum added by fuel mass flow. The default is 0.
    P9 : Float [Pa], optional
        Static exit pressure of core nozzle. The default is None.
    P19 : Float [Pa], optional
        Static exit pressure of bypass nozzle. The default is None.
    Pa : Float [Pa], optional
        Static atmospheric pressure. The default is None.
    A9 : Float [m^2], optional
        Area at the exit of the core nozzle. The default is None.
    A19 : Float [m^2], optional
        Area at the exit of the bypass nozzle. The default is None.
    Returns
    -------
    F : Float [N]
        Thrust of the engine.

    '''
    press_thrust = 0 
    if A9 != None and P9 != None and Pa != None:
        press_thrust += A9*(P9-Pa)
    if A19 != None and P19 != None and Pa != None:
        press_thrust += A19*(P19-Pa) 
        
    F = mdot*(((1+f)/(BPR+1))*(C9-Ca) + (BPR/(BPR+1))*(C19-Ca)) + press_thrust
    return F

def f_1(T03, T04, nb=1, Qf=43500, cpg=1.148, cpa=1.005):
    '''
    Calculates the real Fuel to Air ratio through combustion. If nb = 1 then it is for ideal combustion.

    Parameters
    ----------
    T03 : Flaot [K]
        Stagnation Temperature into the combustor.
    T04 : Float [K]
        Stagnation Temperature out of the combustor.
    nb : Float [0 to 1]
        Combustion efficiency.
    Qf : Float, [kJ/kg] optional
        Fuel energy. The default is 43500 kJ/kg.
    cpg : Float, [kJ/kg*K] optional
        Cp of the combustion gas. The default is 1.148 kJ/kg*K.
    cpa : Float, [kJ/kg*K]
        Cp of the air. The default is 1.005 kJ/kg*K.

    Returns
    -------
    f : Float
        Real fuel to air ratio.

    '''
    num = cpg*T04 - cpa*T03
    den = nb*(Qf - cpg*T04)
    f = num/den
    return f


def TSFC_1(f, BPR, F, mdot):
    '''
    Returns the TSFC for a turbofan engine. Equation given my Professor Cuppoletti

    Parameters
    ----------
    f : Float
        Fuel to air ratio.
    BPR : Float/Int
        Bypass Ratio.
    F : Float/Int
        Total thrust of the engine.
    mdot : Float
        Mass flow rate of air through engine.

    Returns
    -------
    TSFC : Float
        Turbofan Specific Fuel consumption in kg fuel / hr per N of thrust.

    '''
    TSFC = f*3600/((1+BPR)*F/mdot)
    return TSFC


def nT_1(mdot, BPR, f, Ca, C9, C19, Qf=43500):
    '''
    Calculates the Thermal efficiency of a turbofan engine

    Parameters
    ----------
    mdot : Float [kg/s]
        Air mass flow into the engine.
    BPR : Float/Int
        Bypass Ratio.
    f : Float
        Fuel to air ratio.
    Ca : Float [m/s]
        Rel intake velocity.
    C9 : Float [m/s]
        Rel exhaust velocity of the core flow.
    C19 : Float [m/s]
        Rel exhaust velocity of the bypass flow.
    Qf : Float, optional
        Fuel energy, [kJ/kg]. The default is 43500.

    Returns
    -------
    nT : Float [0 - 1]
        Thermal efficiency of a turbofan engine.

    '''
    mdot_h = mdot/(1 + BPR)  # Hot  (core-without fuel)
    mdot_c = mdot - mdot_h   # Cold  (bypass)
    mdot_f = mdot_h*f        # Fuel flow kg/s
    
    num = 0.5*(mdot_h*C9**2 + mdot_c*C19**2 - mdot*Ca**2)
    den = mdot_f*(Qf*1e3)
    nT = num/den
    return nT

def nP_1(mdot, BPR, Ca, C9, C19):
    '''
    Calculates the propulsive efficiency of a turbofan engine

    Parameters
    ----------
    mdot : Float [kg/s]
        Air mass flow into the engine.
    BPR : Float/Int
        Bypass Ratio.
    Ca : Float [m/s]
        Rel intake velocity.
    C9 : Float [m/s]
        Rel exhaust velocity of the core flow.
    C19 : Float [m/s]
        Rel exhaust velocity of the bypass flow.

    Returns
    -------
    nP : Float [0 - 1]
        Propulsive efficiency of a turbofan engine.

    '''
    mdot_h = mdot/(1 + BPR)  # Hot  (core-without fuel)
    mdot_c = mdot - mdot_h   # Cold  (bypass)
    
    num = Ca * (mdot_c*(C19 - Ca) + mdot_h*(C9 - Ca))
    den = 0.5 * (mdot_h*C9**2 + mdot_c*C19**2 - mdot*Ca**2)
    nP = num/den
    return nP

def nO(nT, nP):
    '''
    Calculates the overall efficincy of an engine

    Parameters
    ----------
    nT : Float [0 - 1]
        Thermal efficiency.
    nP : Float [0 - 1]
        Propulsive efficiency.

    Returns
    -------
    nO : Float [0 - 1]
        Overall efficiency.

    '''
    nO = nT*nP
    return nO
    
