# -*- coding: utf-8 -*-
"""
Created on Thu Aug 24 09:45:15 2023

@author: cycon

20240630: 
    Revised and seperated into different files, renamed to SimpleStages.py
    Moved engine definition portion into EngineCompiler.py
    Moved complex stage solving into ComplexStages.py
20250307:
    - Resolved 20250414: Maybe revsise stages to only house one set of gas properties as to unify equations 
    (but values of properties will be different in certain stages)
    - Resolved 20250321: Intake needs ram effect addition and efficiency
    - To Do: Maybe reconfigure isentropic efficiency to be None by default to allow for calc of 
    efficiencies from given temp/pressure values when needed (and force declaration of efficiency)
        20250321 - Resolved, but not complete changes in all stages (made ni=None by default)
20250320:
    - Done: Allow for changing of units to SI/Imperial
    - Done: Update Intake stage to utilize intake pressure ratio if given 
    - Done: Update compressor stage 
    - Resolved: Compressor stage uses an inputted BPR, need to determine a way to allow for mixer to determine BPR
        20250414 (resolved prior but forgot to comment)
        o Fixed on engine-level according to mixer inlet total pressure ratio is unity, aka
        the ratio of total pressures from each stream entering the mixer must be 1, driving the 
        calculation for the BPR.
        o This relation is used to iterate upon the BPR to reach this condition
    - Resolved 20250321: Add ram efficiency in Intake for M > 1 
    - Done: Updated compressor stage case structure to better handle varied inputs
20250407:
    - Done: Updated stages to utilize previous stage's gas properties and allowing for 
    the properties to be set externally (for the combustors/mixers). 
    
"""
import numpy as np
import engine_modules.EngineErrors as EngineErrors
onLaptop = False

import sys 
if onLaptop:
    pathstr = "C:\\Users\\cycon\\Documents\\Modules"
else:
    pathstr = 'C:\\Users\\Cycon Gaming Night\\Documents\\cModules'
sys.path.append(pathstr)
import _aerodynamics.GasDynamics as GD

'''
Naming Convention
X   - Static property X
Xa  - Static property X at atmospheric conditions
Xi  - Static property X at the inlet of a stage
Xe  - Static Property X at the exit of a stage
Xo  - Stagnation (total) property X 
Xoi - Stagnation proerty X at the inlet of a stage
Xoe - Stagnation property X at the exit of a stage


'''
# =============================================================================
# General Stage Description
# =============================================================================
# Use to define the general states/functions shared by each/most stages
class Stage():
    def __init__(self, **kwargs):
        '''
        A general stage class that serves as a baseline for every stage.
        Parameters
        ----------
        **kwargs : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        '''
        self.inputs = kwargs.copy()
        self.UpdateInputs_Gen(**kwargs)
        
    def UpdateInputs_Gen(self,**new_inputs):
        for key in new_inputs.keys():
            self.inputs[key] = new_inputs[key]
            
        # Units
        self.units = self.inputs.get('Units', 'SI')
        # Constants
        # self.gam_a = self.inputs.get('Gamma_a',1.4)                              # Gamma value of air
        # self.gam_g = self.inputs.get('Gamma_g',4/3)                              # Gamma value of gas (combustion products)
        # self.gam_ab = self.inputs.get('Gamma_ab', 1.3)
        # self.cp_a = self.inputs.get('cp_a',1.005 if self.units=='SI' else 0.240) # [kJ/kg*K] or [BTU/lbm*R]  cp of air
        # self.cp_g = self.inputs.get('cp_g',1.148 if self.units=='SI' else 0.295) # [kJ/kg*K] or [BTU/lbm*R]  cp of gas (combustion products)
        # self.cp_ab = self.inputs.get('cp_ab',1.148 if self.units=='SI' else 0.295) # [kJ/kg*K] or [BTU/lbm*R]  cp of gas after burner (combustion products)
        
        # self.R = self.inputs.get('R', self.cp_a*(self.gam_a - 1)/self.gam_a)        # [J/kg*K] or [ft*lbf/R*lbm]    R value of Air
        # self.R_g = self.inputs.get('R_g', self.cp_g*(self.gam_g - 1)/self.gam_g)
        # self.R_ab = self.inputs.get('R_ab', self.cp_ab*(self.gam_ab - 1)/self.gam_ab)
        
        # if self.units != 'SI':
        #     self.R *= 778.16 # Comnvert to [ft*lbf/R*lbm] 
        #     self.R_g *=  778.16
        #     self.R_ab *= 778.16
            
        # Maybe have gas properties passed along to better keep track of what gamma's are in each component
        self.gam_i = self.inputs.get('Gamma_i',None)      
        self.cp_i =  self.inputs.get('cp_i',None) 
        self.R_i =   self.inputs.get('R_i',None)  
        self.gam_e = self.inputs.get('Gamma_e',None) 
        self.cp_e =  self.inputs.get('cp_e',None) 
        self.R_e =   self.inputs.get('R_e',None)  
        
        self.g = self.inputs.get('g',9.81 if self.units=='SI' else 32.174)       # [m/s^2] or [ft/s^2]       Gravitational constant
        self.gc = self.inputs.get('gc',1 if self.units=='SI' else 32.174)        # [N/(kg*m/s^2)]  or [lbm*ft/lbf*s^2] Gravitational Conversion
        self.gf = self.inputs.get('gf',1 if self.units=='SI' else 778.16)         # [kg*m^2/s^2 / J] or [ft*lbf/BTU]
        
        # General properties that could be used by all stages 
        # so all components know the atm conditions
        self.Ta  = self.inputs.get('Ta') 
        self.Pa  = self.inputs.get('Pa')
        self.Vinf = self.inputs.get('Vinf')
        self.Minf = self.inputs.get('Minf')
        
        # All Stage Initial Conditions (not every one is needed)
        self.Toi = self.inputs.get('Toi')
        self.Poi = self.inputs.get('Poi')
        self.Ti  = self.inputs.get('Ti')
        self.Pi  = self.inputs.get('Pi')
        self.Mi  = self.inputs.get('Mi')
        self.Vi  = self.inputs.get('Vi')
        
        # All Stage Exit Conditions
        self.Toe = self.inputs.get('Toe')
        self.Poe = self.inputs.get('Poe')
        self.Te  = self.inputs.get('Te')
        self.Pe  = self.inputs.get('Pe')
        self.Me  = self.inputs.get('Me')
        self.Ve  = self.inputs.get('Ve')
        
        # Other properties
        self.m_dot = self.inputs.get('m_dot') # [kg/s] or [lbm/s] Stays constant through component
        self.mdot_ratio = 1 # used to track mass flow ratio through sections
        self.BPR = self.inputs.get('BPR',0) # mdot_bypass / mdot_core, 0 = no bypass
        
        # Isentropic Efficiency (basically stage efficiency)
        self.ni = self.inputs.get('ni', None) # Isentropic efficiency
        
        # Pressure Ratio and Temperature ratio
        self.pi = self.inputs.get('pi', None) # details total pressure losses in components
        self.tau = self.inputs.get('tau', None)
       
        # Other properties
        # self.StageName = ""
        self.Power = None       # W or BTU/s    - Power
        self.SpecPower = None  #  J/kg or BTU/lbm - Specific Power
        
        self.units_labels = {'V': 'm/s' if self.units=='SI' else 'ft/s',
                             'T': 'K'  if self.units=='SI' else 'R',
                             'P': 'Pa' if self.units=='SI' else 'lbf/ft^2',
                             'mdot': 'kg/s' if self.units=='SI' else 'lbm/s',
                             'Pow': 'W' if self.units=='SI' else 'BTU/s',
                             'Sp_P': 'J/kg' if self.units=='SI' else 'BTU/lbm',
                             'Cp': 'kJ/kg*K' if self.units=='SI' else 'BTU/lbm*R',
                             'R': 'J/kg*K' if self.units=='SI' else 'BTU/lbm*R'}
        
    def forward(self, next_Stage):
        '''
        Propogates the properties of this stage into the next stage.

        Parameters
        ----------
        next_Stage : Stage Object
            A custom child of the stage object that will receive and then proces the inlet properties.

        Returns
        -------
        None.

        '''
        next_Stage.Toi = self.Toe
        next_Stage.Poi = self.Poe
        next_Stage.Ti  = self.Te
        next_Stage.Pi  = self.Pe
        next_Stage.Mi  = self.Me
        next_Stage.Vi  = self.Ve
        next_Stage.m_dot = self.m_dot
        next_Stage.mdot_ratio = self.mdot_ratio
        
        next_Stage.gam_i = self.gam_e if next_Stage.gam_i == None else next_Stage.gam_i 
        next_Stage.cp_i = self.cp_e if next_Stage.cp_i == None else next_Stage.cp_i
        next_Stage.R_i = self.R_e if next_Stage.R_i == None else next_Stage.R_i 
        
    def printOutputs(self, form='{:9.3f}'):
        '''
        Print the outputs of a stage to the console for tracking values and debugging.

        Returns
        -------
        None.

        '''
        print('Stage: ', self.StageName)
        if self.Toe != None:
            print('\t Toe = {} {}'.format(form, self.units_labels['T']).format(self.Toe))
        if self.Poe != None:
            print('\t Poe = {} {}'.format(form, self.units_labels['P']).format(self.Poe))
        if self.Te != None:
            print('\t Te  = {} {}'.format(form, self.units_labels['T']).format(self.Te))
        if self.Pe != None:
            print('\t Pe  = {} {}'.format(form, self.units_labels['P']).format(self.Pe))
        if self.m_dot != None:
            print('\tmdot = {} {}'.format(form, self.units_labels['mdot']).format(self.m_dot))
        if self.mdot_ratio != None:
            print('\tmdot(x/0) = {} '.format(form).format(self.mdot_ratio))
        
        if self.Me != None:
            print('\t Me  = {}'.format(form).format(self.Me))
        if self.Ve != None:
            print('\t Ve  = {} {}'.format(form, self.units_labels['V']).format(self.Ve))
        if self.Power != None:
            print('\t Pow = {} {}'.format(form, self.units_labels['Pow']).format(self.Power))
        if hasattr(self,'specPower'):
            if self.specPower != None:
                print('\t Specific Pow = {} {}'.format(form, self.units_labels['Sp_P']).format(self.specPower))
        
        if self.gam_i == self.gam_e:
            print('\t γ  = {}'.format(form).format(self.gam_i))
        else:
            print('\t γ_i  = {}'.format(form).format(self.gam_i))
            print('\t γ_e  = {}'.format(form).format(self.gam_e))
        
        if self.cp_i == self.cp_e:
            print('\t cp  = {} {}'.format(form,self.units_labels['Cp']).format(self.cp_i))
        else:
            print('\t cp_i  = {} {}'.format(form,self.units_labels['Cp']).format(self.cp_i))
            print('\t cp_e  = {} {}'.format(form,self.units_labels['Cp']).format(self.cp_e))
        
        if self.R_i == self.R_e:
            print('\t R  = {} {}'.format(form,self.units_labels['R']).format(self.R_i))
        else:
            print('\t R_i  = {} {}'.format(form,self.units_labels['R']).format(self.R_i))
            print('\t R_e  = {} {}'.format(form,self.units_labels['R']).format(self.R_e))
        
        self.extraOutputs(form)
    
    def extraOutputs(self,form):
        # Overwrite this and put any extra outputs here within individual stages
        return None
    
    def StageValues(self):
        inputs = {
            'Vi':self.Vi,
            'Mi':self.Mi,
            'Ti':self.Ti,
            'Toi':self.Toi,
            'Pi': self.Pi,
            'Poi':self.Poi,
            'gam_i':self.gam_i, 
            'cp_i': self.cp_i,
            'R_i':self.R_i
            }
        
        outputs = {
            'Ve':self.Ve,
            'Me':self.Me,
            'Te':self.Te,
            'Toe':self.Toe,
            'Pe':self.Pe,
            'Poe':self.Poe,
            'gam_e':self.gam_e, 
            'cp_e': self.cp_e,
            'R_e':self.R_e}
            
        performance = {
            'mdot': self.m_dot,
            'mdot_ratio': self.mdot_ratio,
            'BPR': self.BPR,
            'Power':self.Power,
            'SpecPower': None if not hasattr(self,'specPower') else self.specPower,
            'ni':self.ni,
            'np': None if not hasattr(self,'np') else self.np,
            'pi': self.pi,
            'tau': self.tau,
            'f':None if not hasattr(self,'f') else self.f,
            'BPR': self.BPR
            }
        
        stageVals = {
            'Stage': self.StageName,
            'Units': self.units_labels,
            'inputs': inputs,
            'outputs': outputs,
            'performance': performance}
        return stageVals

# =============================================================================
# Compressor Simple Stage Description
# =============================================================================

class Intake(Stage):
    def __init__(self, **kwargs):
        Stage.__init__(self, **kwargs)
        self.StageName = "Intake"
        # NOTE: Ram efficiency ~= Isentropic Efficiency
        
    def UpdateInputs(self, **kwargs):
        Stage.UpdateInputs_Gen(self, **kwargs)
        
    def calculate(self):
        # Always assume Pi/Pa and Ti/Ta are given (atmos conditions)
        # Cant always assume this, there is a ram affect and streamtubes
        #   occuring outside the intake face
        
        # Check if gas properties were given to intake
        if self.gam_i == None and self.gam_e == None:
            if self.cp_i == None and self.cp_e == None:
                raise Warning('Warning: No gas constants inputted to Intake. Assuming Air')
            else:
                self.gam_i = 1.4  # Gamma value of air
                self.cp_i = 1.005 if self.units=='SI' else 0.240 # [kJ/kg*K] or [BTU/lbm*R]  cp of air
                self.R_i = self.gf*self.cp_i*(self.gam_i - 1)/self.gam_i        # [J/kg*K] or [ft*lbf/R*lbm]    R value of Air
                
        
        self.R_i = self.gf*self.cp_i*(self.gam_i - 1)/self.gam_i        # [J/kg*K] or [ft*lbf/R*lbm]    R value of Air
    
        # Run external conditions
        external = self.external_conditions(self.Minf)
        self.D_additive = external.get('D_add', None)
        self.m_dot = external.get('mdot', self.m_dot)
        self.Mi = external.get('M1', self.Mi)
        self.A0 = external.get('A0', None)
        
        # Ram Affects
        self.eta_r = external['eta_r'] # Pt1/Pt0 
        self.tau_r = 1 + (self.Minf**2)*(self.gam_i - 1)/2       # Totoal/static (Tt0/T0)
        self.pi_r  = self.tau_r**(self.gam_i / (self.gam_i - 1)) # Total/Static  (Pt0/P0)
        # print('RAM: \n\t eta {:.4f}\n\t tau {:.4f}\n\t pi {:.4f}'.format(self.eta_r, self.tau_r,self.pi_r))
        # NOTE: tau_r and pi_r are the only ratios that are
        # static and stagnation ratios
        
        # Correct the diffusor pressure ratio
        self.pi = 1 if self.pi == None else self.pi # Pt2/Pt1
        self.pi_d = self.pi * self.eta_r # total pressure loss from freestream to diffuser exit: Pt2/Pt0
        
        
        # self.Vi = self.Vinf
        self.Mi = self.Minf if self.Mi == None else self.Mi
        
        # If no vel or mach num inputted, assume stationary
        if self.Mi == None:
            if self.Vi == None:
                self.Mi = 0
                self.Vi = 0 
            else:
                self.Mi = self.Vi/np.sqrt(self.gam_i*self.R_i*self.Ti*self.gc)
        # else:
        #     if self.Vi == None:
        #         self.Vi = self.Mi*np.sqrt(self.gam_i*self.R_i*self.Ti*self.gc)
        
        # Now we should have mach num no matter what
        # and the static props (atm props)
        # Find inlet stag props        
        self.Toi = self.Ta * self.tau_r
        self.Poi = self.Pa * self.pi_r * self.eta_r # At station 1 (inlet plane)
        
        self.Pi = self.Poi / GD.Po_P_ratio(self.Mi, self.gam_i) # Corresponds to P1 (static pressure at inlet face, not freestream)
        self.Ti = self.Toi / GD.To_T_ratio(self.Mi, self.gam_i) # Corresponds to T1 (static pressure at inlet face, not freestream
        
        # Find outlet stag props
        self.Toe = self.Toi # Adiabatic
        self.Poe = self.pi * self.Poi # Using only pi_dmax because it is the pressure loss within diffuser. 
                                      # Poi already contains pressure losses from external phenomena
        
        self.tau = self.Toe/self.Toi
        
        self.gam_e = self.gam_i 
        self.cp_e = self.cp_i 
        self.R_e = self.R_i
        
    def external_conditions(self, Minf):
        '''
        Calculates the Ram Efficiency eta_r from the Mach number since
        the equation used will vary depending on freestream Mach. Can be 
        replaced by an outside function that uses the same inputs/outputs
        
        Parameters
        ----------
        Minf : Float
            Freastream Mach Number.
        
        Returns
        -------
        Outputs : Dict
            Necessary:
                'eta_r' - Ram Efficiency
            Recomended:
                'mdot'  - Mass Flow rate into inlet
                'D_add' - Streamtube Additive Drag
                'M1'    - Mach Number at Inlet plane
        
        '''
        if Minf < 0:
            EngineErrors.MissingValue("Minf is Negative", self.StageName)
        elif Minf < 1:
            eta_r = 1 
        elif Minf < 5:
            eta_r = 1 - 0.075*(Minf - 1)**1.35
        else:
            eta_r = 800 / (Minf**4 + 935)
            
        return {'eta_r': eta_r}
    
    def StageValues(self):
        outs = Stage.StageValues(self)
        outs['performance']['tau_r'] = self.tau_r
        outs['performance']['pi_r'] = self.pi_r
        outs['performance']['eta_r'] = self.eta_r
        outs['performance']['pi_d'] = self.pi_d
        outs['performance']['pi_max'] = self.pi
        outs['performance']['A0'] = self.A0
        outs['performance']['D_add'] = self.D_additive
        
        return outs
        
# =============================================================================
# Compressor Simple Stage Description
# =============================================================================
class Compressor(Stage):
    def __init__(self, **kwargs):
        Stage.__init__(self, **kwargs)
        self.UpdateInputs(**kwargs)
        self.StageName = "Compressor"
        # Adding PR and BPR
        # self.r = kwargs.get('rc') # Pressure Ratio of stage (pi)
        
   
    def UpdateInputs(self, **kwargs):
        Stage.UpdateInputs_Gen(self, **kwargs)
        self.BPR = self.inputs.get('BPR', 1) # Bypass Ratio (α): bypass mass flow (air)/mass flow through core (air)
        self.np = self.inputs.get('np') # Polytropic efficiency
        self.mdot_ratio = 1 # Starts as 1 for fan, will be updated by prior comp if
                            # different from 1 from forward section
        self.pi_overall = self.inputs.get('pi_overall')
        
    def calculate(self):
        # Should always have input To and Po, need to calculate power
        # and output To and Po. r will always be given, BPR will affect output 
        # to next stage
        if self.pi == None:
            if self.Poi != None and self.Poe != None:
                self.pi = self.Poe / self.Poi 
            else:
                raise EngineErrors.MissingValue('pi(π) - Pressure Ratio','Compressor')
        
        if self.np == None:
            # If only isentropic efficiency is assumed to be 1 if not entered, send warning
            if self.ni == None:
                # No poly or isen efficiency, check if exit conditions are given
                if self.tau == None and self.Toe == None:
                    # No temp ratio or exit temp given. No way to proceed
                    EngineErrors.MissingValue('ni or np (isentropic or polytropic efficiencies), or tau_c/pi_c.', self.StageName) 
                else: 
                    # have tau or Toe
                    if self.tau == None:
                        self.tau = self.Toe/self.Toi 
                    
                    self.np = (self.gam_i - 1) / (self.gam_i*np.log(self.tau) / np.log(self.pi))
            else:
                # Have isentropic efficiency, no exit conditions
                # Canc calculate np from ni and pi.
                self.np = ((self.gam_i-1)/self.gam_i)*np.log(self.pi) / \
                            np.log( (self.pi**((self.gam_i-1)/self.gam_i) - 1)/self.ni + 1)
                                    
        # Now have np 
        
        n_frac =  (self.gam_i-1)/(self.gam_i*self.np)
        self.Toe = self.Toi + self.Toi*(self.pi**n_frac - 1)
        self.Poe = self.pi*self.Poi
        
        if self.m_dot == None:
            self.specPower = self.mdot_ratio*self.cp_i*(self.Toe-self.Toi)
        else:
            self.Power = self.m_dot*self.cp_i*(self.Toe-self.Toi)
            self.specPower = self.Power/self.m_dot
            
        self.tau = self.Toe/self.Toi if self.tau == None else self.tau
        self.ni = self.calculate_ni_c(self.np, self.gam_i) if self.ni == None else self.ni 
        
        self.gam_e = self.gam_i 
        self.cp_e = self.cp_i 
        self.R_e = self.R_i
        # Done
        
        
    def forward(self, next_Stage_hot, next_Stage_cold=None):
        next_Stage_hot.Toi = self.Toe
        next_Stage_hot.Poi = self.Poe
        next_Stage_hot.Ti  = self.Te
        next_Stage_hot.Pi  = self.Pe
        next_Stage_hot.Mi  = self.Me
        next_Stage_hot.Vi  = self.Ve
        
        next_Stage_hot.gam_i = self.gam_e if next_Stage_hot.gam_i == None else next_Stage_hot.gam_i 
        next_Stage_hot.cp_i = self.cp_e if next_Stage_hot.cp_i == None else next_Stage_hot.cp_i
        next_Stage_hot.R_i = self.R_e if next_Stage_hot.R_i == None else next_Stage_hot.R_i 
        
        # Split airflow if there is bypass after this component
        if next_Stage_cold == None:
            # No Bypass
            next_Stage_hot.m_dot = self.m_dot
            next_Stage_hot.mdot_ratio = self.mdot_ratio
        else:
            # Yes Bypass
            if self.BPR == None:
                raise EngineErrors.MissingValue('BPR','Compressor')
            else:
                if self.m_dot != None:
                    # We have actual mass flow 
                    # Note alpha = BPR = mdot_bypass / mdot_core
                    # mdot = mdot_core + mdot_bypass
                    # mdot_core = mdot_bypass / BPR
                    # mdot_core = (mdot - mdot_core) / BPR
                    # mdot_core*BPR = mdot - mdot_core
                    # mdot_core(BPR + 1) = mdot
                    # mdot_core = mdot/(BPR + 1)
                    m_dot_h = self.m_dot/(self.BPR+1) # Core 
                    m_dot_c = self.m_dot - m_dot_h
                    
                    next_Stage_hot.m_dot = m_dot_h
                    next_Stage_cold.m_dot = m_dot_c
                    
                    # mdot Ratio Calcs
                    mdot_ratio_h = 1/(self.BPR + 1) # Denotes mdot_core / mdot_total
                    mdot_ratio_c = 1 - mdot_ratio_h
                    
                    next_Stage_hot.mdot_ratio = mdot_ratio_h
                    next_Stage_cold.mdot_ratio = mdot_ratio_c

                else:
                    # No inputted mdot
                    mdot_ratio_h = 1/(self.BPR + 1) # Denotes mdot_core / mdot_total
                    mdot_ratio_c = 1 - mdot_ratio_h
                    
                    next_Stage_hot.mdot_ratio = mdot_ratio_h
                    next_Stage_cold.mdot_ratio = mdot_ratio_c
   
                    # Dont need to send mdot ratio to cold section
                    
                next_Stage_cold.Toi = self.Toe
                next_Stage_cold.Poi = self.Poe
                next_Stage_cold.Ti  = self.Te
                next_Stage_cold.Pi  = self.Pe
                next_Stage_cold.Mi  = self.Me
                next_Stage_cold.Vi  = self.Ve
                
                next_Stage_cold.gam_i = self.gam_e if next_Stage_cold.gam_i == None else next_Stage_cold.gam_i 
                next_Stage_cold.cp_i = self.cp_e if next_Stage_cold.cp_i == None else next_Stage_cold.cp_i
                next_Stage_cold.R_i = self.R_e if next_Stage_cold.R_i == None else next_Stage_cold.R_i 
    
        
    def calculate_ni_c(self, np, gamma=1.4):
        '''
        Calculates isentropic efficiency of the compressor from the polytropic efficiency
        Parameters
        ----------
        np : Float, 0-1
            Polytropic efficiency.
        gamma : Float, optional
            Gamma. The default is 1.4.

        Returns
        -------
        Isentropic Efficiency.

        '''
        ni_c = ( self.pi**((self.gam_i-1)/self.gam_i) - 1 ) / ( self.pi**((self.gam_i-1)/(self.gam_i*np)) - 1 )
        return ni_c
        
    def StageValues(self):
        outs = Stage.StageValues(self)
        outs['performance']['pi_overall'] = self.pi_overall
        return outs

    
# =============================================================================
# Combustor Simple Stage Description  
# =============================================================================
class Combustor(Stage):
    def __init__(self, **kwargs):
        Stage.__init__(self, **kwargs)
        self.UpdateInputs(**kwargs)
        self.StageName = "Combustor"
    
    def UpdateInputs(self, **kwargs):
        self.UpdateInputs_Gen(**kwargs) 
        self.dTo = self.inputs.get('dTb')
        # Commenting out to use the pressure ratio pi instead
        # self.dPo = kwargs.get('dPb_dec', 0) # the pressure loss within the compressor as a decimal (0.05 = 5% loss)
        self.f  = self.inputs.get('f') # actual/real fuel-air-ratio 
        self.Q  = self.inputs.get('Q_fuel')
        self.ni = self.inputs.get('nb', 1) # Combustor efficiency
        self.IS_IDEAL = self.inputs.get('IS_IDEAL')
                
    def calculate(self):
        # Assuming we have the Toi and Poi from compressor/prev stage
        # We need to have the exit 
        # Assuming that we have initial gas properties and exit gas properties were set from engine
        if self.R_e == None:
            self.R_e = self.gf*self.cp_e*(self.gam_e - 1)/self.gam_e        # [J/kg*K] or [ft*lbf/R*lbm]    R value of Air
            
        
        
        if self.Toe == None: 
            # No Turbine inlet temp given
            if self.dTo == None: 
                # No combustor increase temp given
                if self.f == None and self.Q == None:
                    # No air-fuel ratio given, cant calculate temps
                    raise EngineErrors.MissingValue('Toe, dTo, or f&Q',self.StageName)
                else: 
                    # We have f and Q to calculate exit temp
                    f_ideal = self.f*self.ni # inputted f would be actual
                    self.Toe = (f_ideal*self.Q + self.cp_i*self.Toi)/(self.cp_e*(1+f_ideal))
            else:
                # We dont have exit temp, but do have temp increase
                self.Toe = self.Toi + self.dTo
         
         # Now we have exit temperature   
        self.pi = 1 if self.pi == None else self.pi 
        self.tau = self.Toe / self.Toi
        self.Poe = self.Poi*self.pi 
        self.dTo = self.Toe - self.Toi # will use later for f calcs
        
        # if self.f == None:
        if self.Q != None: 
            # Assuming non-ideal, will calculate f and use in mass fuel flow
            self.f = (self.cp_e*self.Toe - self.cp_i*self.Toi) / (self.ni*(self.Q - self.cp_e*self.Toe))
                
        # Note: f = mdot_fuel / mdot_core_air
        if not self.IS_IDEAL:
            # Only change mass flow rate if not ideal case
            if self.m_dot != None and self.f != None:
                self.m_dot += self.f*self.m_dot
                self.mdot_ratio += self.f*self.mdot_ratio
            elif self.m_dot == None and self.f != None:
                self.mdot_ratio += self.f*self.mdot_ratio
                
        # final mdot_ratio is mdot_combustor exit / mdot_total_air
        # BPR = alpha = mdot_bypass / mdot_core 
        # its current mdot_ratio = mdot_core / mdot_total_air 
        # 
        # f = mdot_fuel / mdot_core_air
        # mdot_be = mdot_c + f*mdot_c 
        # mdot_be/mdot_air_total = (1+f)mdot_core/mdot_air = (1+f)/(1 + BPR)
            
    def extraOutputs(self,form):
        if self.f != None:
            print('\t f = {}'.format(form).format(self.f))
        
            
# =============================================================================
# Turbine Simple Stage Description
# =============================================================================
class Turbine(Stage):
    def __init__(self, Comp_to_power, **kwargs):
        Stage.__init__(self, **kwargs)
        self.UpdateInputs(Comp_to_power, **kwargs)
        self.StageName = "Turbine"

    def UpdateInputs(self, Comp_to_power, **kwargs):
        self.UpdateInputs_Gen(**kwargs) 
        self.np = self.inputs.get('np') # Polytropic efficiency
        self.nm = self.inputs.get('nm',1) # Mechanical efficiency
        self.Compressor = Comp_to_power # Could be list
        # Will have inlet temp, compressor power
        #self.r  = kwargs.get('rt') # Add for later, not used now
        # this will be for generators or when turbine pressure ratio is specified
        

    def calculate(self):
        if self.m_dot != None:
            if type(self.Compressor) == list:
                com_power = 0
                for i in range(0,len(self.Compressor)):
                    com_power += self.Compressor[i].Power
                self.Power = com_power/self.nm 
            else:
                self.Power = self.Compressor.Power/self.nm 
            # Calculate exit temp
            self.Toe = self.Toi - self.Power/(self.m_dot*self.cp_i)
            self.specPower = self.Power/self.m_dot
        else:
            # No m_dot is given, need to power balance based on 
            # BPR ratios instead
            if type(self.Compressor) == list:
                com_power = 0
                for i in range(0,len(self.Compressor)):
                    com_power += self.Compressor[i].specPower
                self.specPower = com_power/self.nm 
            else:
                self.specPower = self.Compressor.specPower/self.nm 
            # Calculate exit temp
            self.Toe = self.Toi - self.specPower/(self.mdot_ratio*self.cp_i)
        
            
        if self.np == None:
            if self.pi != None:
                # Calculate np
                self.np = np.log(1- self.ni*(1 - self.pi**((self.gam_i-1)/self.gam_i)))
                self.np /= np.log(self.pi)*(self.gam_i-1)/self.gam_i
            else:
                print('Warning: insufficient parameters given to turbine')
                print('Continuing assuming polytropic efficiency = 1')
                self.np = 1
                
        m_frac = self.np*(self.gam_i-1)/self.gam_i
        self.Poe = self.Poi*(1- (self.Toi-self.Toe)/self.Toi )**(1/m_frac)
        
        self.pi = self.Poe/self.Poi 
        self.tau = self.Toe/self.Toi 
        
        self.gam_e = self.gam_i 
        self.cp_e = self.cp_i 
        self.R_e = self.R_i
    # Done
        
    
    
# =============================================================================
# Mizer Simple Stage Description    
# =============================================================================
class Mixer(Stage):
    def __init__(self, CoreMixStage, BypassMixStage, **kwargs):
        Stage.__init__(self, **kwargs)
        self.UpdateInputs(CoreMixStage, BypassMixStage, **kwargs)
        self.StageName = "Mixer"

    def UpdateInputs(self, CoreMixStage, BypassMixStage, **kwargs):
        self.CoreMix = CoreMixStage
        self.BypassMix = BypassMixStage
        self.UpdateInputs_Gen(**kwargs) 
        self.BPR_f = self.inputs.get('BPR_f')
        self.MIX_GAS_PROPERTIES = self.inputs.get('MIX_GAS_PROPERTIES', True)
        
    def calculate(self):
        if self.Mi == None and self.CoreMix.Me == None:
            raise EngineErrors.MissingValue('Missing Coreflow Mach Number',self.StageName)
            
        # FOR PROJECT 1:
        # We have M6 (core exit Mach) and pi_M (pressure ratio for mixer)
        # Assumes that gamma for bypass is the same for air 
        
        self.mdot_ratio = self.CoreMix.mdot_ratio + self.BypassMix.mdot_ratio
        if self.m_dot != None:
            self.m_dot = self.CoreMix.m_dot + self.BypassMix.m_dot
        
        # M_16 calculation
        # Gonna just rename values to make double checking equation easier
        M_6 = self.CoreMix.Me 
        gam_6 = self.CoreMix.gam_e
        gam_16 = self.BypassMix.gam_e 
        cp_6 = self.CoreMix.cp_e
        cp_16 = self.BypassMix.cp_e 
        R_6   = self.CoreMix.R_e 
        R_16  = self.CoreMix.R_e
        P_t6 = self.CoreMix.Poe 
        P_t16 = self.BypassMix.Poe 
        T_t6 = self.CoreMix.Toe 
        T_t16 = self.BypassMix.Toe 
        
        M_16 = np.sqrt( (2/(gam_16 - 1)) * (( (P_t16/P_t6)*(1 + (M_6**2)*(gam_6 - 1)/2)**(gam_6/(gam_6-1))  )**((gam_16 - 1)/gam_16)  - 1)  )
        if abs(self.BPR) < 1e-4:
            M_16 = 0.01
        self.BypassMix.Me = M_16 
        
        self.cp_6a = (cp_6 + self.BPR_f*cp_16) / (1 + self.BPR_f)
        self.R_6a = (R_6 + self.BPR_f*R_16) / (1 + self.BPR_f)
        self.gam_6a = self.cp_6a*self.gf / (self.cp_6a*self.gf - self.R_6a)
        
        self.tau = (cp_6 / self.cp_6a) * (1 + self.BPR_f * (cp_16/cp_6)*(self.BypassMix.Toe/self.CoreMix.Toe)) / (1 + self.BPR_f)
        
        # if self.pi==None:
            # Use ideal pi 
        # A16_A6 = self.BPR_f*(P_t6/P_t16)*np.sqrt(T_t16/T_t6)*(GD.MassFlowParam_norm(M_6 ,gam_6)/(GD.MassFlowParam_norm(M_16,gam_16)))
        # Calculate A16/A6
        num = gam_16*R_6 * (1 + (M_16**2)*(gam_16-1)/2)
        den = gam_6*R_16 * (1 + (M_6**2)*(gam_6-1)/2)
        A16_A6 = self.BPR_f*np.sqrt(T_t16/T_t6)*(M_6/M_16) / np.sqrt(num/den)
        self.A16_A6 = A16_A6
        # Calculate phi(M6A, gam_6a):
        num = 1 + self.BPR_f 
        den = 1/np.sqrt(GD.Rayleigh_phi_MS(M_6,gam_6)) + self.BPR_f*np.sqrt(R_16*gam_6*(T_t16/T_t6)/(R_6*gam_16*GD.Rayleigh_phi_MS(M_16,gam_16)))
        phi = ((num/den)**2)* (self.R_6a*gam_6*self.tau)/(R_6*self.gam_6a)
        # Calculate M_6A
        M_6A = np.sqrt(2*phi / (1 - 2*self.gam_6a*phi + np.sqrt(1 - 2*(self.gam_6a+1)*phi)))
       
        self.pi_ideal = ((1+self.BPR_f)*np.sqrt(self.tau)/(1+ A16_A6))*(GD.MassFlowParam_norm(M_6 ,gam_6)*np.sqrt(self.gc/R_6) / (GD.MassFlowParam_norm(M_6A,self.gam_6a)*np.sqrt(self.gc/self.R_6a)))

        self.pi_M = self.pi_ideal *self.pi
        self.Me = M_6A
            
        self.Toe = self.tau * self.CoreMix.Toe 
        self.Poe = self.pi * self.CoreMix.Poe
       
        if self.MIX_GAS_PROPERTIES:
            self.gam_e = self.gam_6a
            self.cp_e = self.cp_6a
            self.R_e = self.R_6a
        else:
            self.gam_e = gam_6 
            self.cp_e = cp_6 
            self.R_e = R_6 
        ## NEED TO GET ALL VALUES FOR THIS
        # alpha_m = (eta_m * (1 + f) * (tau_t / tau_c) * (1 - (pi_f/(pi_c*pi_b))**((gamma_t - 1)*e_t/gamma_t)) - (tau_c - 1)) / (tau_f - 1)
        
        return None
    
    def StageValues(self):
        outs = Stage.StageValues(self)
        outs['performance']['pi_ideal'] = self.pi_ideal
        outs['performance']['pi_max'] = self.pi
        outs['performance']['pi_M'] = self.pi_M
        outs['performance']['A16/A6'] = self.A16_A6
        return outs
    





# =============================================================================
# Nozzle Simple Stage Description
# =============================================================================
class Nozzle(Stage):
    def __init__(self, nozzle_type='CD', **kwargs):
        Stage.__init__(self, **kwargs)
        # if air_type == 'hot':
        #     self.gam = self.gam_g
        #     self.R = self.R_g
        # else:
        #     self.gam = self.gam_a
        #     self.R = self.R
        self.nozzle_type = nozzle_type # 'C' for Converging, 'CD' for Conv-Div
        self.UpdateInputs(**kwargs)
        self.StageName = "Nozzle"

    def UpdateInputs(self, **kwargs):
        self.UpdateInputs_Gen(**kwargs) 
        
    def calculate(self):
        # Check if choked
        Tc = self.Toi*(2/(self.gam_i+1))
        ni = 1 # TEMPORARY BC CURRENTLY NOT SUPPORTED TO NOT HAVE IT AN INPUT
        Pc = self.Poi*(1 - (1/ni)*(1-Tc/self.Toi))**(self.gam_i/(self.gam_i-1))
        
        P_rat = self.Poi/self.Pa
        P_crit = self.Poi/Pc
        if P_rat > P_crit:
            # Nozzle is choked
            if self.nozzle_type == 'C':
                self.Pe = Pc # For Converging nozzle
            else:
                # CD Nozzle
                self.Pe = self.Pa 
                # is because Pe = Pa for a fully expanded CD Nozzle
                # This part may not be correct 
                # self.Te = self.Toi*(1-self.ni*(1-(self.Pe/self.Poi)**((self.gam-1)/self.gam)))
        else:
            # Nozzle is not choked
            self.Pe = self.Pa
        
        self.NPR = self.Poi / self.Pa
        if self.pi != None:
            # Can calculate isentropic efficiency
            self.ni = ((self.NPR*(self.Pe/self.Pa))**((self.gam_i-1)/self.gam_i) - self.pi**((1-self.gam_i)/self.gam_i)) \
                / ((self.NPR*(self.Pe/self.Pa))**((self.gam_i-1)/self.gam_i) - 1)
        elif self.pi == None and self.ni != None:
            self.pi =  ((self.NPR*(self.Pe/self.Pa))**((self.gam_i-1)/self.gam_i) - (self.ni*((self.NPR(self.Pe/self.Pa))**((self.gam_i-1)/self.gam_i) - 1))) ** (self.gam_i/(1-self.gam_i))
            
        else:
            EngineErrors.IncompleteInputs("Missing eta_n or pi_n", self.StageName)
        
        self.Poe = self.pi * self.Poi
        self.Me = GD.Mach_at_PR(self.Poe/self.Pe,Gamma = self.gam_i)
        
        self.Te = self.Toi*(1-self.ni*(1-(self.Pe/self.Poi)**((self.gam_i-1)/self.gam_i)))
        self.Ve = self.Me*np.sqrt(self.gam_i*self.R_i*self.Te*self.gc)
        
        # Stag props at exit
        self.Toe = self.Toi
        self.Poe = self.Pe * (1 + (self.gam_i -1)*(self.Me**2)/2)**(self.gam_i/(self.gam_i -1))
        self.tau = self.Toe/self.Toi
        
        self.gam_e = self.gam_i 
        self.cp_e = self.cp_i 
        self.R_e = self.R_i 
        
   

class Duct(Stage):
    def __init__(self, **kwargs): 
        # Adiabatic duct without friction currently
        Stage.__init__(self, **kwargs)
        self.UpdateInputs(**kwargs)
        self.StageName = 'Duct'
        
    def UpdateInputs(self,**kwargs):
        self.UpdateInputs_Gen(**kwargs) 
    
    def calculate(self):
        self.Poe = self.pi * self.Poi 
        self.Toe = self.Toi
        self.tau = 1 
        
        self.gam_e = self.gam_i 
        self.cp_e = self.cp_i 
        self.R_e = self.R_i
        
    
