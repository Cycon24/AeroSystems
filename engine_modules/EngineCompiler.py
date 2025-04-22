# -*- coding: utf-8 -*-
"""
Created on Sun Jun 30 07:37:41 2024

@author: cycon

To-Do List:
    - Define the EngineCreator init fnct
        - Determine how to handle the parameters.
            - There should be two types of params, overall/intra params (engine/between component params)
            and specific/inter params (within a component).
            - Overall params are set within EngineCreator module
        - 
    - Restructure calculate to handle n number of airstreams and m number components
    - Create some sort of packet of data that can be transfered to the database/returned
        - Look at Mareks code regarding these packets
        
        
20240630:
    Created from previously named file EngineModule.py, which is renamed to SimpleStages.py
20240714:
    
"""

import engine_modules.EnginePerformanceFunctions as EPF
import engine_modules.SimpleStages as SS
import numpy as np


# =============================================================================
# Main Class that can generate any engine 
# =============================================================================

class EngineCreator():
    def __init__(self, **kwargs):
        print('bla bla')
        
        self.MAIN_CALCULATION_UPDATED = False 
        
    
    def calculate(self, printVals=True):
        '''
        Calculates the properties of the air throughout the engine.

        Parameters
        ----------
        printVals : Bool, optional
            When true, it will print out each components exit conditions. The default is True.

        Returns
        -------
        None.

        '''
        for i in range(0,len(self.AllStages)):
            # Calculate each row and print outputs
            self.AllStages[i][0].calculate()
            if printVals: self.AllStages[i][0].printOutputs()
            # Check if current stage has a parallel (ie, prev stage passes air to 2 stages)
            if self.AllStages[i][1] != None:
                self.AllStages[i][1].calculate()
                if printVals: self.AllStages[i][1].printOutputs()
                
            # Move forward/propogate
            if i != len(self.AllStages)-1: # It is not at the end, so forward
                if self.AllStages[i+1][1] != None: 
                    # Means that this stage delivers to two stages: fan -> HPC & BP Noz
                    self.AllStages[i][0].forward(self.AllStages[i+1][0],self.AllStages[i+1][1])
                else:
                    # Stage delivers to one stage
                    self.AllStages[i][0].forward(self.AllStages[i+1][0])
                    
        self.MAIN_CALCULATION_UPDATED = True 
                    
    def engine_performance(self):
        
        # Uninstalled thrust at the moment
        # Installed thrust will include inlet and outlet losses
        if not self.MAIN_CALCULATION_UPDATED:
            self.calculate(False) 
            
        # Pull out types of stages, 
        
        
        
        
            
# =============================================================================
# Example Classes
# =============================================================================
class Turbofan_SingleSpool():
    def __init__(self, **kwargs):
        '''
        A signle spool turbofan which has one turbine to power the
        fan and compressor. It has a cold-air bypass which is after the fan
        and goes straight to a nozzle. The core-flow goes to a second compressor
        and then to the combustor, followed by a single turbine and lastly the
        core nozzle. 
        Parameters
        ----------
        **kwargs : Dictionary
            Contains all needed and optional parameters with the keys listed below.
            Required:
            'Ta': Atmospheric static temperature
            'Pa': Atmospheric static pressure
            'rfan': Fan Pressure Ratio
            'rc':   Compressor Pressure Ratio
            'BPR':  Bypass Ratio (If not passed, assumed 1 so no bypass occurs)
            'T_turb_in': Turbine inlet temp
            
            Optional
            'Vinf': None, # Or Minf
            'Minf': None, # Or Vinf, if none its assumed stationary
            'mdot_a': # Mass flow rate of air into engine (kg/s)
            'Q_fuel': Heat energy of fuel, pass if real to calculate f and include fuel flow in power/velecity calcs
            'F': Thrust of the engine produced
            Efficiencies (Assumed to be 1 if not passed)
            'ni': Inlet Isentropic Efficiency
            'nj': Nozzle Isentropic Efficiency
            'nf': Fan Isentropic Efficiency
            'nc': Compressor Isentropic Efficiency
            'nt': Turbine Isentropic Efficiency
            'nb': Cobustor Efficincy
            'nm': Mechanical Efficiency
            'npf': Fan Polytropic Efficiency (overrides isentropic)
            'npc': Compressor Polytropic Efficiency (overrides isentropic)
            'npt': Turbine Polytropic Efficiency (overrides isentropic)
            'dP_combustor': Decimal pressure drop in combustor (ex: 0.05 for 5% P loss, 0 for ideal)
            
        Returns
        -------
        None.

        '''
        # Stages
        # Atm moving
        # Inlet
        # Fan  (is a compressor)
        # Bypass Nzzle
        # LP Compressor
        # HP Compressor
        # Combustor
        # HP Turbine
        # LP Turbine
        # Nozzle
        self.inputs = kwargs.copy()
        gen_kwargs = {
            'Ta': kwargs.get('Ta'),
            'Pa': kwargs.get('Pa'),
            'Vinf': kwargs.get('Vinf'),
            'Minf': kwargs.get('Minf'),
            'BPR': kwargs.get('BPR',1), # Need this here on case mdot=None
            'Q_fuel':  kwargs.get('Q_fuel')}# kJ/kg
        # Efficiencies
        ni = kwargs.get('ni',1) # Inlet
        nj = kwargs.get('nj',1) # Nozzle
        nf = kwargs.get('nf',1) # Compressor - Isentropic
        nc = kwargs.get('nc',1) # Compressor - Isentropic
        nt = kwargs.get('nt',1) # Turbine - Isentropic
        nb = kwargs.get('nb',1) # Cobustor
        nm = kwargs.get('nm',1) # Mechanical
        npf = kwargs.get('npf') # Fan - Polytropic
        npc = kwargs.get('npc') # Compressor - Polytropic
        npt = kwargs.get('npt') # Turbine - Polytropic
        # Pressure Ratios/Relations
        dP_b = kwargs.get('dP_combustor') # Decimal pressure drop in combustor
        rfan = kwargs.get('rfan') # Fan PR
        rc   = kwargs.get('rc')   # Compressor PR
        # Turbine Inlet
        To_ti = kwargs.get('T_turb_in') # K - Turbine inlet temp
        # Air Mass flow
        mdot = kwargs.get('mdot_a') # kg/s
        
        self.F = kwargs.get('F')
        
        # Define each stage and pass in parameters
        self.inlet     = SS.Intake(**gen_kwargs,ni=ni,m_dot=mdot)
        self.fan       = SS.Compressor(**gen_kwargs, rc=rfan, np=npf, ni=nf)
        self.BP_nozzle = SS.Nozzle('cold',**gen_kwargs, ni=nj)
        self.HP_comp   = SS.Compressor(**gen_kwargs, rc=rc, np=npc, ni=nc)
        self.combustor = SS.Combustor(**gen_kwargs, Toe=To_ti, dPb_dec=dP_b, ni=nb)
        self.HP_turb   = SS.Turbine([self.fan, self.HP_comp], **gen_kwargs, nm=nm, ni=nt, np=npt)
        self.nozzle    = SS.Nozzle(**gen_kwargs, ni=nj) # Nozzle/Exhaust?
        
        # Set names for easier readout checks
        self.fan.StageName = 'Fan'
        self.BP_nozzle.StageName = 'Cold Nozzle'
        self.nozzle.StageName = 'Hot Nozzle'
        
        # Define all stages in engine to iterate through
        # Two dimensional since there is a bypass, ie one stage
        # passes params to two different stages
        self.AllStages = [[self.inlet, None ],
                          [self.fan, None], 
                          [self.HP_comp,  self.BP_nozzle],
                          [self.combustor,None],
                          [self.HP_turb, None],
                          [self.nozzle, None]]
        
    def calculate(self, printVals=True):
        '''
        Calculates the properties of the air throughout the engine.

        Parameters
        ----------
        printVals : Bool, optional
            When true, it will print out each components exit conditions. The default is True.

        Returns
        -------
        None.

        '''
        for i in range(0,len(self.AllStages)):
            # Calculate each row and print outputs
            self.AllStages[i][0].calculate()
            if printVals: self.AllStages[i][0].printOutputs()
            # Check if current stage has a parallel (ie, prev stage passes air to 2 stages)
            if self.AllStages[i][1] != None:
                self.AllStages[i][1].calculate()
                if printVals: self.AllStages[i][1].printOutputs()
                
            # Move forward/propogate
            if i != len(self.AllStages)-1: # It is not at the end, so forward
                if self.AllStages[i+1][1] != None: 
                    # Means that this stage delivers to two stages: fan -> HPC & BP Noz
                    self.AllStages[i][0].forward(self.AllStages[i+1][0],self.AllStages[i+1][1])
                else:
                    # Stage delivers to one stage
                    self.AllStages[i][0].forward(self.AllStages[i+1][0])
                    
    def getOutputs(self):
        '''
        Returns a dictionary containing important values fromwithin the engine
        using typical engine notaion. Must be run after calculate in order to contain values.

        Returns
        -------
        outs : Dictionary
            Conatins the following items:
                'mdot_c' - bypass flow [kg/s]
                'mdot_h1' - core flow  [kg/s]
                'mdot_h2' - core flow + fuel [kg/s]
                'mdot' - air flow into the engine [kg/s]
                'Ca'  - Initial vel of air (rel to engine) [m/s]
                'C9'  - Exhast vel of core    [m/s]
                'C19' - Exhuast vel of bypass [m/s]
                'Pa'   - Static atmospheric pressure [Pa]
                'P9'   - Exit pressure of core noz   [Pa]
                'P19'  - Exit pressure of bypass noz [Pa]
                'To3': - Combustor inlet stagnation temp [K]
                'To4': - Combustor outlet stagnation temp [K]
                'nb':  - combustor efficiency
                'dH':  - fuel energy  [kJ/kg]
                'f':   - air to fuel ratio

        '''
        outs = {
            'mdot_c': self.BP_nozzle.m_dot, # bypass flow
            'mdot_h1': self.HP_comp.m_dot, # core flow
            'mdot_h2': self.nozzle.m_dot, # core flow + fuel
            'mdot': self.inlet.m_dot, # air flow in
            'Ca': self.inlet.Vi, # Initial vel of air
            'C9': self.nozzle.Ve, # Exhast vel of core
            'C19': self.BP_nozzle.Ve, # Exhuast vel of bypass
            'Pa': self.inputs.get('Pa'),
            'P9': self.nozzle.Pe,    # Exit pressure of core noz
            'P19': self.BP_nozzle.Pe,# Exit pressure of bypass noz
            'To3': self.combustor.Toi, # combustor inlet temp
            'To4': self.combustor.Toe, # combustor outlet temp
            'nb': self.combustor.ni, # combustor efficiency
            'dH': self.combustor.Q,  # fuel energy
            'f': self.combustor.f    # air to fuel ratio
            }
        return outs
    
    def printInputs(self):
        '''
        Prints out all of the kwargs entered on intialization

        Returns
        -------
        None.

        '''
        print('Inputs')
        for key,val in self.inputs.items():
            print('\t {}  =  {}'.format(key,val))
            
    def calculatePerformanceParams(self):
        outs = self.getOutputs()
        
        # Calculate thrust or mdot if one is given and not the other
        if self.F == None:
            if outs['mdot'] != None:
                # Calcualte thrust
                self.F = EPF.Thrust_1(outs['mdot'], self.inputs['BPR'], outs['C9'], outs['C19'], outs['Ca'])
        else:
            if outs['mdot'] == None:
                # calculate mdot
                mdot = EPF.mdot_2(self.F, self.inputs('BPR'), outs['C9'], outs['C19'], outs['Ca'])
                
                
                



class Turbofan_DoubleSpool():
    def __init__(self, **kwargs):
        '''
        A double spool turbofan which has two turbines, one to power the
        fan and another for the compressor. It has a cold-air bypass which is after the fan
        and goes straight to a nozzle. The core-flow goes to a second compressor
        and then to the combustor, followed by a single turbine and lastly the
        core nozzle. 
        Parameters
        ----------
        **kwargs : Dictionary
            Contains all needed and optional parameters with the keys listed below.
            Required:
            'Ta': Atmospheric static temperature
            'Pa': Atmospheric static pressure
            'rfan': Fan Pressure Ratio
            'rc':   Compressor Pressure Ratio
            'BPR':  Bypass Ratio (If not passed, assumed 1 so no bypass occurs)
            'T_turb_in': Turbine inlet temp (HP Turbine)
            
            Optional
            'Vinf':  Or Minf
            'Minf':  Or Vinf, if none its assumed stationary
            'mdot_a': Mass flow rate of air into engine (kg/s)
            'Q_fuel': Heat energy of fuel, pass if real to calculate f and include fuel flow in power/velecity calcs
            'F': Thrust of the engine produced
            Efficiencies (Assumed to be 1 if not passed)
            'ni': Inlet Isentropic Efficiency
            'nj': Nozzle Isentropic Efficiency
            'nf': Fan Isentropic Efficiency
            'nc': Compressor Isentropic Efficiency
            'nt': Turbine Isentropic Efficiency
            'nb': Cobustor Efficincy
            'nm': Mechanical Efficiency
            'npf': Fan Polytropic Efficiency (overrides isentropic)
            'npc': Compressor Polytropic Efficiency (overrides isentropic)
            'npt': Turbine Polytropic Efficiency (overrides isentropic)
            'npt_lp': Low Pressure Turbine Polytropic Efficiency (overrides isentropic)
            'dP_combustor': Decimal pressure drop in combustor (ex: 0.05 for 5% P loss, 0 for ideal)

        Returns
        -------
        None.

        '''
        # Stages
        # Atm moving
        # Inlet
        # Fan  (is a compressor)
        # Bypass Nzzle
        # LP Compressor
        # HP Compressor
        # Combustor
        # HP Turbine
        # LP Turbine
        # Nozzle
        self.inputs = kwargs.copy()
        gen_kwargs = {
            'Ta': kwargs.get('Ta'),
            'Pa': kwargs.get('Pa'),
            'Vinf': kwargs.get('Vinf'),
            'Minf': kwargs.get('Minf'),
            'BPR': kwargs.get('BPR',1), # Need this here on case mdot=None
            'Q_fuel':  kwargs.get('Q_fuel')}# kJ/kg
        # Efficiencies
        ni = kwargs.get('ni',1) # Inlet
        nj = kwargs.get('nj',1) # Nozzle
        nf = kwargs.get('nf',1) # Compressor - Isentropic
        nc = kwargs.get('nc',1) # Compressor - Isentropic
        nt = kwargs.get('nt',1) # Turbine - Isentropic
        nt_lp = kwargs.get('nt_lp',1) # LP Turbine - Isentropic
        nb = kwargs.get('nb',1) # Cobustor
        nm = kwargs.get('nm',1) # Mechanical
        npf = kwargs.get('npf') # Fan - Polytropic
        npc = kwargs.get('npc') # Compressor - Polytropic
        npt = kwargs.get('npt') # Turbine - Polytropic
        npt_lp = kwargs.get('npt_lp') # LP Turbine - Polytropic
        # Pressure Ratios/Relations
        dP_b = kwargs.get('dP_combustor') # Decimal pressure drop in combustor
        rfan = kwargs.get('rfan') # Fan PR
        rc   = kwargs.get('rc')   # Compressor PR
        # Turbine Inlet
        To_ti = kwargs.get('T_turb_in') # K - Turbine inlet temp
        # Air Mass flow
        mdot = kwargs.get('mdot_a') # kg/s
        
        self.F = kwargs.get('F')
        
        # Define each stage and pass in parameters
        self.inlet     = SS.Intake(**gen_kwargs,ni=ni,m_dot=mdot)
        self.fan       = SS.Compressor(**gen_kwargs, rc=rfan, np=npf, ni=nf)
        self.BP_nozzle = SS.Nozzle('cold',**gen_kwargs, ni=nj)
        self.HP_comp   = SS.Compressor(**gen_kwargs, rc=rc, np=npc, ni=nc)
        self.combustor = SS.Combustor(**gen_kwargs, Toe=To_ti, dPb_dec=dP_b, ni=nb)
        self.HP_turb   = SS.Turbine(self.HP_comp, **gen_kwargs, nm=nm, ni=nt, np=npt)
        self.LP_turb   = SS.Turbine(self.fan, **gen_kwargs, nm=nm, ni=nt_lp, np=npt_lp)
        self.nozzle    = SS.Nozzle(**gen_kwargs, ni=nj) # Nozzle/Exhaust?
        
        # Set names for easier readout checks
        self.fan.StageName = 'Fan'
        self.BP_nozzle.StageName = 'Cold Nozzle'
        self.nozzle.StageName = 'Hot Nozzle'
        self.HP_turb.StageName = 'HP Turbine'
        self.LP_turb.StageName = 'LP Turbine'
        
        # Define all stages in engine to iterate through
        # Two dimensional since there is a bypass, ie one stage
        # passes params to two different stages
        self.AllStages = [[self.inlet, None ],
                          [self.fan, None], 
                          [self.HP_comp,  self.BP_nozzle],
                          [self.combustor,None],
                          [self.HP_turb, None],
                          [self.LP_turb, None],
                          [self.nozzle, None]]
        
    def calculate(self, printVals=True):
        '''
        Calculates the properties of the air throughout the engine.

        Parameters
        ----------
        printVals : Bool, optional
            When true, it will print out each components exit conditions. The default is True.

        Returns
        -------
        None.

        '''
        for i in range(0,len(self.AllStages)):
            # Calculate each row and print outputs
            self.AllStages[i][0].calculate()
            if printVals: self.AllStages[i][0].printOutputs()
            # Check if current stage has a parallel (ie, prev stage passes air to 2 stages)
            if self.AllStages[i][1] != None:
                self.AllStages[i][1].calculate()
                if printVals: self.AllStages[i][1].printOutputs()
                
            # Move forward/propogate
            if i != len(self.AllStages)-1: # It is not at the end, so forward
                if self.AllStages[i+1][1] != None: 
                    # Means that this stage delivers to two stages: fan -> HPC & BP Noz
                    self.AllStages[i][0].forward(self.AllStages[i+1][0],self.AllStages[i+1][1])
                else:
                    # Stage delivers to one stage
                    self.AllStages[i][0].forward(self.AllStages[i+1][0])
                    
    def getOutputs(self):
        '''
        Returns a dictionary containing important values fromwithin the engine
        using typical engine notaion. Must be run after calculate in order to contain values.

        Returns
        -------
        outs : Dictionary
            Conatins the following items:
                'mdot_c' - bypass flow [kg/s]
                'mdot_h1' - core flow  [kg/s]
                'mdot_h2' - core flow + fuel [kg/s]
                'mdot' - air flow into the engine [kg/s]
                'Ca'  - Initial vel of air (rel to engine) [m/s]
                'C9'  - Exhast vel of core    [m/s]
                'C19' - Exhuast vel of bypass [m/s]
                'Pa'   - Static atmospheric pressure [Pa]
                'P9'   - Exit pressure of core noz   [Pa]
                'P19'  - Exit pressure of bypass noz [Pa]
                'To3': - Combustor inlet stagnation temp [K]
                'To4': - Combustor outlet stagnation temp [K]
                'nb':  - combustor efficiency
                'dH':  - fuel energy  [kJ/kg]
                'f':   - air to fuel ratio

        '''
        outs = {
            'mdot_c': self.BP_nozzle.m_dot, # bypass flow
            'mdot_h1': self.HP_comp.m_dot, # core flow
            'mdot_h2': self.nozzle.m_dot, # core flow + fuel
            'mdot': self.inlet.m_dot, # air flow in
            'Ca': self.inlet.Vi, # Initial vel of air
            'C9': self.nozzle.Ve, # Exhast vel of core
            'C19': self.BP_nozzle.Ve, # Exhuast vel of bypass
            'Pa': self.inputs.get('Pa'),
            'P9': self.nozzle.Pe,    # Exit pressure of core noz
            'P19': self.BP_nozzle.Pe,# Exit pressure of bypass noz
            'To3': self.combustor.Toi, # combustor inlet temp
            'To4': self.combustor.Toe, # combustor outlet temp
            'nb': self.combustor.ni, # combustor efficiency
            'dH': self.combustor.Q,  # fuel energy
            'f': self.combustor.f    # air to fuel ratio
            }
        return outs
    
    def printInputs(self):
        '''
        Prints out all of the kwargs entered on intialization

        Returns
        -------
        None.

        '''
        print('Turbofan Engine Inputs')
        for key,val in self.inputs.items():
            print('\t {}  =  {}'.format(key,val))
            
    def calculatePerformanceParams(self):
        outs = self.getOutputs()
        
        # Calculate thrust or mdot if one is given and not the other
        if self.F == None:
            if outs['mdot'] != None:
                # Calcualte thrust
                self.F = EPF.Thrust_1(outs['mdot'], self.inputs['BPR'], outs['C9'], outs['C19'], outs['Ca'])
        else:
            if outs['mdot'] == None:
                # calculate mdot
                mdot = EPF.mdot_2(self.F, self.inputs('BPR'), outs['C9'], outs['C19'], outs['Ca'])
                
                
                
                

class Turbojet_AfterBurner():
    def __init__(self, **kwargs):
        '''
        A single spool turbojet engine with a fan with bypass, HPC, burner, HPT,
        mixer, afterburner and nozzle. 
        Parameters
        ----------
        **kwargs : Dictionary
            Contains all needed and optional parameters with the keys listed below.
            Required:
            'Ta': Atmospheric static temperature
            'Pa': Atmospheric static pressure
            'pi_f' Fan Pressure Ratio
            'pi_c':   Compressor Pressure Ratio
            'T_turb_in': Turbine inlet temp (HP Turbine)
            
            Optional
            'Vinf':  Or Minf
            'Minf':  Or Vinf, if none its assumed stationary
            'mdot_a': Mass flow rate of air into engine (kg/s) (lbm/s)
            'Q_fuel': Heat energy of fuel, pass if real to calculate f and include fuel flow in power/velecity calcs
            'F': Thrust of the engine produced
            Efficiencies (Assumed to be 1 if not passed)
            'ni': Inlet Isentropic Efficiency
            'nj': Nozzle Isentropic Efficiency
            'nf': Fan Isentropic Efficiency
            'nc': Compressor Isentropic Efficiency
            'nt': Turbine Isentropic Efficiency
            'nb': Cobustor Efficincy
            'nm': Mechanical Efficiency
            'npf': Fan Polytropic Efficiency (overrides isentropic)
            'npc': Compressor Polytropic Efficiency (overrides isentropic)
            'npt': Turbine Polytropic Efficiency (overrides isentropic)
            'pi_b': Total pressure ratio through the burner

        Returns
        -------
        None.

        '''
        # Stages
        # Atm moving
        # Inlet
        # Fan  (is a compressor)
        # Bypass
        # HP Compressor
        # Combustor
        # HP Turbine
        # Mixer
        # Afterburner
        # Nozzle
        self.inputs = kwargs.copy()
        self.UpdateInputs(ON_INITIATION=True, **kwargs)
        
        
        
    def UpdateInputs(self, ON_INITIATION=False, **new_inputs):
        for key in new_inputs.keys():
            self.inputs[key] = new_inputs[key]
        
        
        self.units = self.inputs.get('Units', 'SI')
        # Efficiencies
        eta_b = self.inputs.get('eta_b') # Cobustor
        eta_ab = self.inputs.get('eta_ab') # Afterburner
        eta_m = self.inputs.get('eta_m') # Mechanical
        npf = self.inputs.get('e_f') # Fan - Polytropic
        npc = self.inputs.get('e_c') # Compressor - Polytropic
        npt = self.inputs.get('e_t') # Turbine - Polytropic
        # Pressure Ratios/Relations
        pi_d = self.inputs.get('pi_d') # Diffuser total pressure ratio
        pi_b = self.inputs.get('pi_b') # Combustor total pressure ratio
        pi_f = self.inputs.get('pi_f') # Fan PR
        pi_overall = self.inputs.get('pi_c')   # Compressor PR
        pi_c  = pi_overall/pi_f
        pi_M  = self.inputs.get('pi_M') # Mixer total pressure ratio
        pi_n  = self.inputs.get('pi_n') # Nozzle total pressure ratio
        pi_AB = self.inputs.get('pi_ab') # Afterburner total ressure raito
        
        # Turbine Inlet / Combustor Outlet
        To_ti = self.inputs.get('Tt4') # K - Turbine inlet temp
        To_ab_e = self.inputs.get('Tt7')
        dTo_ab  = self.inputs.get('dTt_ab', None)
        M6 = self.inputs.get('M6') # Needs to be set to turbine exit mach
        # Air Mass flow
        mdot = self.inputs.get('mdot_a') # kg/s or lbm/s
        
        # Should have some overall engine restrictions
        BPR_Min = 0 
        Compressor_Max_Toe = 3000 
        Throttle_Setting = 1 # %
        
        self.gen_kwargs = {
            'Units': self.units, 
            'Ta': self.inputs.get('Ta'),
            'Pa': self.inputs.get('Pa'),
            'Vinf': self.inputs.get('Vinf'),
            'Minf': self.inputs.get('Minf'),
            'Q_fuel':  self.inputs.get('h_PR'),# kJ/kg
            'IS_IDEAL': self.inputs.get('IS_IDEAL')
            }
        cp_air = self.inputs.get('cp_c')
        cp_b = self.inputs.get('cp_t')
        cp_ab =  self.inputs.get('cp_ab')
        gam_air = self.inputs.get('Gamma_c')
        gam_b =  self.inputs.get('Gamma_t')
        gam_ab = self.inputs.get('Gamma_ab')
            
        
        # self.F = kwargs.get('F')
        
        
        if ON_INITIATION:
            # Define each stage and pass in parameters
            self.Inlet     = SS.Intake(**self.gen_kwargs, m_dot=mdot, pi=pi_d, cp_i=cp_air, Gamma_i=gam_air)
            self.Fan       = SS.Compressor(**self.gen_kwargs, pi=pi_f, np=npf)
            self.BP_duct   = SS.Duct(**self.gen_kwargs, pi=1)
            self.Compressor   = SS.Compressor(**self.gen_kwargs, pi=pi_c, np=npc, pi_overall=pi_overall)
            self.Combustor = SS.Combustor(**self.gen_kwargs, Toe=To_ti, pi=pi_b, ni=eta_b, cp_e=cp_b, Gamma_e=gam_b)
            self.Turbine   = SS.Turbine([self.Compressor, self.Fan], **self.gen_kwargs, nm=eta_m, np=npt)
            self.Mixer     = SS.Mixer(self.Turbine, self.BP_duct, **self.gen_kwargs, pi=pi_M, MIX_GAS_PROPERTIES=self.inputs.get('MIX_GAS_PROPERTIES'))
            self.Afterburner = SS.Combustor(**self.gen_kwargs, pi=pi_AB, ni=eta_ab, dTb=dTo_ab, cp_e=cp_ab, Gamma_e=gam_ab)
            self.Nozzle    = SS.Nozzle(air_type='hot',nozzle_type='CD',**self.gen_kwargs, pi=pi_n) 
            
            # Set names for easier readout checks
            self.Fan.StageName = 'Fan'
            self.Afterburner.StageName = 'Afterburner'
            
            # Set other conditions
            self.Combustor.Toe = To_ti # Set combustor outlet temperature
            self.Turbine.Me = M6
            self.Afterburner.Toe = To_ab_e 
            
            
            # Define all stages in engine to iterate through
            # Two dimensional since there is a bypass, ie one stage
            # passes params to two different stages
            self.AllStages = [[self.Inlet, None ],
                              [self.Fan, None], 
                              [self.Compressor,  self.BP_duct],
                              [self.Combustor,None],
                              [self.Turbine, None],
                              [self.Mixer, None],
                              [self.Afterburner, None],
                              [self.Nozzle, None]]
            
        else:
            # Only pass values into each stage rather than all
             self.Inlet      .UpdateInputs(**self.gen_kwargs, m_dot=mdot, pi=pi_d)
             self.Fan        .UpdateInputs(**self.gen_kwargs, pi=pi_f, np=npf)
             self.BP_duct    .UpdateInputs(pi=1)
             self.Compressor    .UpdateInputs(**self.gen_kwargs, pi=pi_c, np=npc,pi_overall=pi_overall)
             self.Combustor  .UpdateInputs(**self.gen_kwargs, Toe=To_ti, pi=pi_b, ni=eta_b)
             self.Turbine    .UpdateInputs([self.Compressor, self.Fan], **self.gen_kwargs, nm=eta_m, np=npt)
             self.Mixer      .UpdateInputs(self.Turbine, self.BP_duct, **self.gen_kwargs, pi=pi_M)
             self.Afterburner.UpdateInputs(**self.gen_kwargs, pi=pi_AB, ni=eta_ab)
             self.Nozzle     .UpdateInputs(**self.gen_kwargs, pi=pi_n) 
             
             # Set other conditions
             self.Combustor.Toe = To_ti # Set combustor outlet temperature
             self.Turbine.Me = M6
             self.Afterburner.Toe = To_ab_e 
        
             
    def calculate(self, printVals=True):
        '''
        Calculates the properties of the air throughout the engine.

        Parameters
        ----------
        printVals : Bool, optional
            When true, it will print out each components exit conditions. The default is True.

        Returns
        -------
        None.

        '''
        # Get gas props
        cp_a = self.Inlet.cp_i 
        cp_g = self.Combustor.cp_e
        gam_a = self.Inlet.gam_i 
        gam_g = self.Combustor.gam_e
        Ta =  self.inputs.get('Ta')
        To_ti = self.Combustor.Toe
        tau_lambda = cp_g*To_ti / (cp_a*Ta)
        tau_r = 1 + (self.inputs.get('Minf')**2)*(gam_a - 1)/2
        tau_c = self.Compressor.pi_overall**((gam_a-1)/(gam_a*self.Compressor.np))
        tau_f = self.Fan.pi**((gam_a-1)/(gam_a*self.Fan.np))
        
        # print(locals())
        # Get first-pass guess on alpha, will refine down the road
        f = cp_a*Ta *(tau_lambda - tau_r*tau_c) / self.inputs['h_PR']
        # alpha_f = (self.inputs['eta_m']*(1+f) * (tau_lambda/tau_c)*(1 - (self.Fan.pi/(self.Compressor.pi*self.Combustor.pi))**((gam_g-1)*self.Turbine.np/gam_g)) - (tau_c-1)) / (tau_f-1)
        # alpha = alpha_f*(1+f)   
        alpha = (tau_lambda*(tau_c-tau_f))/(tau_r*tau_c*(tau_f-1))  - (tau_c-1)/(tau_f-1)
        alpha_f = alpha / (1+f)
       
        
        # self.inputs['BPR'] = alpha # Need this here on case mdot=None
        # self.inputs['BPR_f'] = alpha_f 
        
        # # Pass into needed components
        # self.Fan.BPR = self.inputs['BPR']
        # self.Combustor.BPR = self.inputs['BPR']
        # self.Combustor.f = f 
        # self.Mixer.BPR_f = self.inputs['BPR_f']
        error = 0
        dPt = 1
        tol = 1e-4
        while abs(dPt) > tol:
            alpha += error 
            # check if alpha went < 0
            if alpha < 0:
                # print('Alpha set to 0')
                alpha = 0 
                alpha_f = 0 
            else:
                alpha_f = alpha / (1+f) if self.Combustor.f == None else alpha / (1+self.Combustor.f)
            # print(' f = {}\n alpha = {}\n alpha_f = {}'.format(f, alpha,alpha_f))
            
            self.inputs['BPR'] = alpha # Need this here on case mdot=None
            self.inputs['BPR_f'] = alpha_f 
            
            # Pass into needed components
            self.Fan.BPR = self.inputs['BPR']
            self.Combustor.BPR = self.inputs['BPR']
            # self.Combustor.f = f 
            self.Mixer.BPR_f = self.inputs['BPR_f']
            
            for i in range(0,len(self.AllStages)):
                # Calculate each row and print outputs
                self.AllStages[i][0].calculate()
                # Check if current stage has a parallel (ie, prev stage passes air to 2 stages)
                if self.AllStages[i][1] != None:
                    self.AllStages[i][1].calculate()
                    # if printVals: self.AllStages[i][1].printOutputs()
                    
                # Move forward/propogate
                if i != len(self.AllStages)-1: # It is not at the end, so forward
                    if self.AllStages[i+1][1] != None: 
                        # Means that this stage delivers to two stages: fan -> HPC & BP Noz
                        self.AllStages[i][0].forward(self.AllStages[i+1][0],self.AllStages[i+1][1])
                    else:
                        # Stage delivers to one stage
                        self.AllStages[i][0].forward(self.AllStages[i+1][0])
            
            # Check error
            if abs(alpha) < tol:
                dPt = 0 
                error = 0 
                # print('Warning: Bypass Ratio reached 0...')
            else:
                dPt = self.Turbine.Poe - self.BP_duct.Poe 
                error = dPt/(self.Turbine.Poe + self.BP_duct.Poe) #/self.Turbine.Poe 
            
            # print('dPt = {:.6f}   err = {:.6f}'.format(dPt,error))
            # if (abs(dPt) < tol) and alpha < 0:
            #     # Within error tolerance but with a negative BPR
            #     self.Fan.pi *= 0.99 # Since alpha negative, reduccing fan pi will increase alpha on recalc
            #     alpha += 0.1 # Fix BPR to positive
            #     dPt = 1 # Adjsut dPt to not exit loop
                
    
        
        if printVals:
            for i in range(0,len(self.AllStages)):   
                 self.AllStages[i][0].printOutputs()
                 if self.AllStages[i][1] != None:
                     self.AllStages[i][1].printOutputs()
            
    def RunParameterSweep(self, paramKey, paramList, perfFunctions=None, printVals=False):
        # Setup the output dictionary
        # Needs to contain an item for each perfFunction
        # Then the values within each perfFunction needs to be lists. 
        # Lengths of lists based on param list, this needs to be defined
        #   prior to running sweep
        
        arrayFormat = np.zeros((len(paramList),))
        # Make output dic
        finalOutputs = {paramKey: paramList}
        
        # If perf functions arent a list, make them one
        if type(perfFunctions) != list:
            perfFunctions = [perfFunctions]
        for func in perfFunctions:
            finalOutputs[func.__name__] = self.OutputsToEmptyOutputsArray(arrayFormat, func)
        
        
        
        for i, p in enumerate(paramList):
            # Update inputs with parameter
            # self.inputs[paramKey] = p
            self.UpdateInputs(**{paramKey: p})
            
            # Run Calculation Loop
            self.calculate(printVals=printVals)
            
            # Run output functions
            if perfFunctions != None:
                # Iterate through functions                  
                for func in perfFunctions:
                    # Get function outputs
                    func_outs = func()
                    # Run through each thing in outputs
                    for key in func_outs.keys():
                        if type(func_outs[key]) == dict:
                            # func outputs nested dict (assuming only 3 layers)
                            for key2 in func_outs[key].keys():
                                if type(func_outs[key][key2]) == dict:
                                    for key3 in func_outs[key][key2].keys():
                                        if type(func_outs[key][key2][key3]) != str:
                                            finalOutputs[func.__name__][key][key2][key3][i] = func_outs[key][key2][key3]
                                        else:
                                            finalOutputs[func.__name__][key][key2][key3] = func_outs[key][key2][key3]
                                
                                elif type(func_outs[key][key2]) != str:
                                    # Set an empty list as the value
                                    finalOutputs[func.__name__][key][key2][i] = func_outs[key][key2]
                                else:
                                    finalOutputs[func.__name__][key][key2] = func_outs[key][key2]
                        elif type(func_outs[key]) != str:
                            finalOutputs[func.__name__][key][i] = func_outs[key]
                        else:
                            # Is a string
                            finalOutputs[func.__name__][key] = func_outs[key]
                        
                                   
        return finalOutputs
        

    def OutputsToEmptyOutputsArray(self, arrayFormat, func):
        '''
        Returns a dictionary containing lists in the same format as
        the function 'func' outputs (which outputs a dict with single values)

        Parameters
        ----------
        arrayFormat : np Array
            An empty numpy array the same size as needed for param sweep.
        func : Function
            The performance/output function.

        Returns
        -------
        funcOuts : Dict
            Contains dict with values as lists.

        '''
        # Run the function
        self.calculate(False)
        funcOuts = func()
        for key in funcOuts.keys():
            if type(funcOuts[key]) == dict:
                # func outputs nested dict (assuming only 3 layers)
                for key2 in funcOuts[key].keys():
                    if type(funcOuts[key][key2]) == dict:
                        for key3 in funcOuts[key][key2].keys():
                            if type(funcOuts[key][key2][key3]) != str:
                                funcOuts[key][key2][key3] = arrayFormat.copy()
                            
                    elif type(funcOuts[key][key2]) != str:
                        # Set an empty list as the value
                        funcOuts[key][key2] = arrayFormat.copy()
                        
            elif type(funcOuts[key]) != str:
                funcOuts[key] = arrayFormat.copy()
                
        return funcOuts
                                  
               
    
    def printInputs(self):
        '''
        Prints out all of the kwargs entered on intialization

        Returns
        -------
        None.

        '''
        print('Turbofan Engine Inputs')
        for key,val in self.inputs.items():
            print('\t {}  =  {}'.format(key,val))
            
   
    def CalculatePerformanceParams(self):
        '''
        Need to get: 
        F/mdot_0
        SFC = f_tot / F/mdot0
        f
        f_AB
        f_0 = f_b / (1+ alpha) + f_ab 
        eta_T = 1 - 1 / (tau_r*tau_c)
        eta_P = 2*Minf * (V9/a0 - Minf) / ((V9/a0)^2 - M0^2) 
        eta_O = T*P
        alpha

        Returns
        -------
        Dict.
            'F_mdot': Specific Thrust
            'SFC': Thrust Specific Fuel Consumption
            'f': Fuel-Air Ratio of Combustor
            'f_AB': Fuel-Air Ratio of Afterburner
            'f_tot': FUel-Air Ratio of total engine
            'eta_T': Thermal Efficiency
            'eta_P': Propulsion Efficiency
            'eta_O': Overall Efficiency
            'alpha': Bypass Ratio (mdot_bypass / mdot_core)

        '''
        # Thrust, since Pe = Pa no need to incorperate 
        gam_a = self.inputs.get('Gamma_c')
        gc = self.Inlet.gc
        gf = self.Inlet.gf
        R = self.Inlet.R_i
        Ta = self.inputs.get('Ta')
        a0 = np.sqrt(gam_a*R*Ta*gc)
        Minf = self.inputs.get('Minf')
        V0 = Minf*a0
        h_PR = self.inputs.get('h_PR')
        
        F_mdot = (a0/gc)*(self.Nozzle.mdot_ratio*self.Nozzle.Ve/a0 - Minf)
        f_b = self.Combustor.f
        f_ab = self.Afterburner.f
        alpha = self.Fan.BPR 
        f_tot = f_b/(1 + alpha) + f_ab
        SFC = 3600*f_tot / F_mdot
        

        mdot_N = self.Nozzle.m_dot
        Thrust = None if self.Nozzle.m_dot == None else F_mdot * mdot_N
        Installed_Thrust =  None if self.Inlet.D_additive == None or Thrust == None else Thrust - self.Inlet.D_additive
        T_mdot = None if Installed_Thrust == None else Installed_Thrust / self.Inlet.m_dot
        mdot_f_tot =  None if self.Inlet.m_dot == None else f_tot * self.Inlet.m_dot
        TSFC =  None if Installed_Thrust == None else 3600*mdot_f_tot / Installed_Thrust
        
        eta_T = 0.5*(self.Nozzle.mdot_ratio*self.Nozzle.Ve**2 - V0**2) / (f_tot*h_PR*gf*gc)
        #1 - 1 / (self.Inlet.tau_r * self.Compressor.tau)
        eta_P = 2*Minf * (self.Nozzle.Ve/a0 - Minf) / ((self.Nozzle.Ve/a0)**2 - Minf**2)
        eta_O = eta_T * eta_P 
        
        TPR = self.Inlet.eta_r * self.Inlet.pi * self.Inlet.pi_r
        
        # STILL INCORRECT FOR SOME REASON
        # tau_lambda = (self.Turbine.cp_e/self.Compressor.cp_e)*(self.Combustor.Toe/self.inputs['Ta'])
        # alpha_2 = ((tau_lambda/(self.Inlet.tau_r*self.Fan.tau))*(1+f_b)*(self.Turbine.tau-1) + self.Fan.tau - self.Compressor.tau/self.Fan.tau)/(1-self.Fan.tau)
        outputs = {
            'F_mdot': F_mdot,
            'T_mdot': T_mdot,
            'F': Thrust,
            'T': Installed_Thrust,
            'S':SFC,
            'TSFC': TSFC, 
            'f_b':f_b,
            'f_ab':f_ab,
            'f_tot':f_tot,
            'eta_T':eta_T,
            'eta_P':eta_P,
            'eta_O':eta_O,
            'alpha':alpha,
            'TPR': TPR}
        # print('M0 = {:.2f}\tα1 = {:.3f}\tα2 = {:.3f}'.format(Minf, alpha, alpha_2))
        return outputs
    
    def getStageStagnationVals(self):
        StageStagnationProps = {}
        for i in range(0,len(self.AllStages)):   
             stageOuts = self.AllStages[i][0].StageValues()['outputs']
             StageStagnationProps[self.AllStages[i][0].StageName] = {'Poe': stageOuts['Poe'], 'Toe': stageOuts['Toe']}
             if self.AllStages[i][1] != None:
                 stageOuts = self.AllStages[i][1].StageValues()['outputs']
                 StageStagnationProps[self.AllStages[i][1].StageName] = {'Poe': stageOuts['Poe'], 'Toe': stageOuts['Toe']}
        return StageStagnationProps
    
    def getStageVals(self):
        StageProps = {} 
        
        for i in range(0,len(self.AllStages)):   
             stageOuts = self.AllStages[i][0].StageValues()
             StageProps[self.AllStages[i][0].StageName] = stageOuts
             
             if self.AllStages[i][1] != None:
                 stageOuts = self.AllStages[i][1].StageValues()
                 StageProps[self.AllStages[i][1].StageName] = stageOuts
        return StageProps