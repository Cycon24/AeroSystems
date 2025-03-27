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

import EnginePerformanceFunctions as EPF
import SimpleStages as SS


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
        self.UpdateInputs(**kwargs, ON_INITIATION=True)
        
        
        
    def UpdateInputs(self, new_inputs, ON_INITIATION=False):
        for key in new_inputs.keys():
            self.inputs[key] = new_inputs[key]
        
        
        self.units = self.inputs.get('units', 'SI')
        
        # Efficiencies
        eta_b = self.inputs.get('nb') # Cobustor
        eta_m = self.inputs.get('nm') # Mechanical
        npf = self.inputs.get('npf') # Fan - Polytropic
        npc = self.inputs.get('npc') # Compressor - Polytropic
        npt = self.inputs.get('npt') # Turbine - Polytropic
        # Pressure Ratios/Relations
        pi_d = self.inputs.get('pi_d') # Diffuser total pressure ratio
        pi_b = self.inputs.get('pi_b') # Combustor total pressure ratio
        pi_f = self.inputs.get('pi_f') # Fan PR
        pi_c   = self.inputs.get('pi_c')   # Compressor PR
        pi_M  = self.inputs.get('pi_M') # Mixer total pressure ratio
        pi_n  = self.inputs.get('pi_n') # Nozzle total pressure ratio
        pi_AB = self.inputs.get('pi_AB') # Afterburner total ressure raito
        
        # Turbine Inlet
        To_ti = self.inputs.get('T_turb_in') # K - Turbine inlet temp
        # Air Mass flow
        mdot = self.inputs.get('mdot_a') # kg/s or lbm/s
        
        self.gen_kwargs = {
            'Ta': self.inputs.get('Ta'),
            'Pa': self.inputs.get('Pa'),
            'Vinf': self.inputs.get('Vinf'),
            'Minf': self.inputs.get('Minf'),
            'Q_fuel':  self.inputs.get('Q_fuel')}# kJ/kg
            'cp_a'
            'cp_g'
            'cp_ab'
            'Gamma_a'
            'Gamma_g'
            'Gamma_ab'
            
            }
            
        
        # self.F = kwargs.get('F')
        
        
        if ON_INITIATION:
            # Define each stage and pass in parameters
            self.inlet     = SS.Intake(**gen_kwargs, m_dot=mdot, pi=pi_d)
            self.fan       = SS.Compressor(**gen_kwargs, pi=pi_f, np=npf)
            self.BP_duct   = SS.Duct(pi=1)
            self.HP_comp   = SS.Compressor(**gen_kwargs, pi=pi_c, np=npc)
            self.combustor = SS.Combustor(**gen_kwargs, Toe=To_ti, pi=pi_b, ni=eta_b)
            self.HP_turb   = SS.Turbine([self.HP_comp, self.fan], **gen_kwargs, nm=eta_m, np=npt)
            self.Mixer     = SS.Mixer(self.HP_turb, self.BP_duct, pi=pi_M)
            self.Afterburner = SS.Combustor(**gen_kwargs, )
            self.nozzle    = SS.Nozzle(**gen_kwargs) # Nozzle/Exhaust?
            
            # Set names for easier readout checks
            self.fan.StageName = 'Fan'
            self.nozzle.StageName = 'Hot Nozzle'
            self.HP_turb.StageName = 'Turbine'
            self.Afterburner.StageName = 'Afterburner'
            
            # Define all stages in engine to iterate through
            # Two dimensional since there is a bypass, ie one stage
            # passes params to two different stages
            self.AllStages = [[self.inlet, None ],
                              [self.fan, None], 
                              [self.HP_comp,  self.BP_duct],
                              [self.combustor,None],
                              [self.HP_turb, None],
                              [self.Mixer, None],
                              [self.Afterburner, None]
                              [self.nozzle, None]]
            
        else:
            # Only pass values into each stage rather than all
            
        
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
        cp_a = simp.cp_a 
        cp_g = simp.cp_g
        gam_a = simp.gam_a 
        gam_g = simp.gam_g
        Ta =  kwargs.get('Ta')
        tau_lambda = cp_g*To_ti / (cp_a*Ta)
        tau_r = 1 + (kwargs.get['Minf']**2)*(gam_a + 1)/2
        tau_c = pi_c**((gam_a-1)/(gam_a*npc))
        tau_f = pi_f**((gam_a-1)/(gam_a*npf))
        
        
        f = cp_a*Ta *(tau_lambda - tau_r*tau_c) / gen_kwargs['Q_fuel']
        alpha_f = (eta_m*(1+f) * (tau_lambda/tau_c)*(1 - (pi_f/(pi_c*pi_b))**((gam_g-1)*npt/gam_g)) - (tau_c-1)) / (tau_f-1)
        alpha = alpha_f*(1+f)    
        
        gen_kwargs['BPR'] = alpha # Need this here on case mdot=None
        
        
        
        
        
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
                # mdot = EPF.mdot_2(self.F, sel)
                return None 
            