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
20251122:
    - Noticed that there are errors in Nozzle setup, specifically for Pc calculation
    that doesnt take into account ni which is a function of pi if provided.
    - Done: Added ability to directly define "NextStage" and "BypassStages" for more complex engine geometries
    
    
"""
import numpy as np
import engine_modules.EngineErrors as EngineErrors
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

20251009 - Need to add more structure and separate things better. Beter documentation
            and explanation on each input for each stage
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
        self._update_inputs_base(**kwargs)
        
    def _update_inputs_base(self,**new_inputs):
        '''
        Re-gathers inputs and saves them into self.inputs dict to allow for 
        property variation and propagation. _update_inputs_base not called
        outside of specified stages in which it is called. To update a specific
        stage use UpdateInputs(). This is to ensure stage-specific variables
        that are inputed are correctly handled

        Parameters
        ----------
        **new_inputs : dict or kwargs
            Any new or updated kargs for the inputs of the stage.

        Returns
        -------
        None.

        '''
        # Save/update all inputs
        for key in new_inputs.keys():
            self.inputs[key] = new_inputs[key]
            
        # =============================================================================
        #  General Properties (same for all stages, unless directly altered)      
        # =============================================================================
        # Units
        self.units = self.inputs.get('Units', 'SI')
        self.IS_IDEAL = self.inputs.get('IS_IDEAL', None)
        
        # General properties that could be used by all stages 
        # so all components know the atm conditions
        self.Ta  = self.inputs.get('Ta') 
        self.Pa  = self.inputs.get('Pa')
        self.Vinf = self.inputs.get('Vinf')
        self.Minf = self.inputs.get('Minf')
                
        # =============================================================================
        #   Gas constants      
        # =============================================================================
        # Gas constants initialized as none and passed through each stage
        # IS_IDEAL will ensure gas properties do not change.
        self.gam_i = self.inputs.get('Gamma_i',None)      
        self.cp_i =  self.inputs.get('cp_i',None) # [kJ/kg*K]
        self.R_i =   self.inputs.get('R_i',None)  # [J/kg*K]
        self.gam_e = self.inputs.get('Gamma_e',None) 
        self.cp_e =  self.inputs.get('cp_e',None) 
        self.R_e =   self.inputs.get('R_e',None)  
        
        self.g = self.inputs.get('g',9.81 if self.units=='SI' else 32.174)       # [m/s^2] or [ft/s^2]       Gravitational constant
        self.gc = self.inputs.get('gc',1 if self.units=='SI' else 32.174)        # [N/(kg*m/s^2)]  or [lbm*ft/lbf*s^2] Gravitational Conversion
        self.gf = self.inputs.get('gf',1000 if self.units=='SI' else 778.16)         # [kg*m^2/s^2 / J] or [ft*lbf/BTU]
        # DEBUG WARNING gf above default value may need to be restored to 1, unsure
        # what other calculations are included but this will convert R calculated later into
        # J rather than kJ (seems to be okay - 20251009)
        
       
        # =============================================================================
        #  Stage Properties    
        # =============================================================================
        # All Stage Initial Conditions, static properties not typically calculated
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
        self.mdot = self.inputs.get('mdot') # [kg/s] or [lbm/s] Stays constant through component
        self.mdot_ratio = 1 # used to track mass flow ratio through sections
        self.BPR = self.inputs.get('BPR',0) # mdot_bypass / mdot_core, 0 = no bypass
        
        
        # Isentropic Efficiency (basically stage efficiency)
        self.ni = self.inputs.get('ni', None) # Isentropic efficiency
        
        # Pressure Ratio and Temperature ratio
        self.pi = self.inputs.get('pi', None) # details total pressure losses in components
        self.tau = self.inputs.get('tau', None)
       
        # Other properties
        self.StageName = self.inputs.get('StageName', "Stage")
        self.StageID = self.inputs.get('StageID', "s")
        # self.StageType = "Stage"
        self.Power = None       # W or BTU/s    - Power
        self.SpecPower = None  #  J/kg or BTU/lbm - Specific Power
        
        # Physical Dimensions
        self.Ai = self.inputs.get('Ai', None)
        self.Ae = self.inputs.get('Ae', None)
        self.Ai_mdota = None 
        self.Ae_mdota = None
        
        # Stage passes, to overwrite in forward
        self.NextStage = self.inputs.get("NextStage", None)
        self.BypassStages= self.inputs.get("BypassStages", [None])
        
        
        # =============================================================================
        #   Unit Dict to be referenced if needed       
        # =============================================================================
        self.units_labels = {'V': 'm/s' if self.units=='SI' else 'ft/s',
                             'T': 'K'  if self.units=='SI' else 'R',
                             'P': 'Pa' if self.units=='SI' else 'lbf/ft^2',
                             'mdot': 'kg/s' if self.units=='SI' else 'lbm/s',
                             'Pow': 'W' if self.units=='SI' else 'BTU/s',
                             'Sp_P': 'J/kg' if self.units=='SI' else 'BTU/lbm',
                             'Cp': 'kJ/kg*K' if self.units=='SI' else 'BTU/lbm*R',
                             'R': 'J/kg*K' if self.units=='SI' else 'BTU/lbm*R',
                             'A': 'm^2' if self.units=='SI' else 'ft^2'}
        
    def forward(self, core_Stage: object, bypass_Stages: list = [None]):
        '''
        Propogates the properties of this stage into the next stage/stages.

        Parameters
        ----------
        core_Stage : Stage Object
            A custom child of the stage object that will receive and then proces the inlet properties.
        bypass_stages : list of Stage Objects
            A list of stages that air is bypassed to, if len > 1 then multiple BPRs must be defined
        Returns
        -------
        None.

        '''
        # How each property type is passed along:
        #  - Stagnation: Same to all next stages
        #  - Static: Currently not passed (20251009 - I think its best if its calculated internally)
        #  - Mach Number: Currently not passed, same as above 
        #  - Mdot/ratio: Passed to each based on BPR for that stage
        #  - gas properties: Now passed unconditionally (20251009 - previously dependend on if next
        #                   stage defined them, but input gammas for a stage will not be)
        
        # Check if next and bypass stages were defined, overwrite if so:
        if self.NextStage != None:
            core_Stage = self.NextStage 
        if self.BypassStages != [None]:
            bypass_Stages = self.BypassStages 
        
        # pass to Core flow
        core_Stage.Toi = self.Toe
        core_Stage.Poi = self.Poe
        core_Stage.gam_i = self.gam_e 
        core_Stage.cp_i = self.cp_e 
        core_Stage.R_i = self.R_e 
        
        if not isinstance(bypass_Stages,list):
            bypass_Stages = [bypass_Stages]
        
        # Check if theres bypass stages
        if bypass_Stages[0] == None:
            # No bypass, pass as regular
            core_Stage.mdot = self.mdot
            core_Stage.mdot_ratio = self.mdot_ratio
        else:
            # Theres bypass, ensure BPR is a list
            if type(self.BPR) != list:
                self.BPR = [self.BPR] # Turn into list
            # Check if lengths match:
            if len(self.BPR) != len(bypass_Stages):
                raise ValueError(f"[Error]\t {self.StageName}: List mismatch with BPRs provided for each bypass stage")
            # Check if any BPRs are None
            if None in self.BPR:
                raise ValueError(f"[Error]\t {self.StageName}: Bypass Ratio for one or more bypass stages is None")
            
            # Now everything has been checked, loop through to calculate mass flow rate to pass into each stage
            # calculate the sum of BPRs
            sum_BPRs = 0
            for alpha_n in self.BPR: sum_BPRs += alpha_n 
            
            
            # Pass information into core stream
            core_Stage.mdot = self.mdot / (1 + sum_BPRs) if self.mdot != None else None
            core_Stage.mdot_ratio = self.mdot_ratio / (1 + sum_BPRs)
            
            
            # Iterate through each bypass stage and pass needed variables
            for n, bypass in enumerate(bypass_Stages):
                bypass.mdot = self.BPR[n] * core_Stage.mdot if self.mdot != None else None
                bypass.mdot_ratio = self.BPR[n] * core_Stage.mdot_ratio
                bypass.Toi = self.Toe
                bypass.Poi = self.Poe
                bypass.gam_i = self.gam_e 
                bypass.cp_i = self.cp_e 
                bypass.R_i = self.R_e 
                    
                    
                    
        # next_Stage.Ti  = self.Te
        # next_Stage.Pi  = self.Pe
        # next_Stage.Mi  = self.Me
        # next_Stage.Vi  = self.Ve
        # next_stage.gam_i = self.gam_e if next_Stage.gam_i == None else next_Stage.gam_i 
        
        return None
        
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
        if self.mdot != None:
            print('\tmdot = {} {}'.format(form, self.units_labels['mdot']).format(self.mdot))
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
    
    def calcDetailProps_i(self, AssumeSubsonic=True):
        '''
        Calculates and updates all static, area, or Mach number values at the
        inlet of the stage.
        
        Currently inefficient but should work
        
        Returns
        -------
        None.

        '''

        
        # Always assume we have stagnation props, impossible to run outside code if not
        
        # Check if we have gas properties, if not return none
        haveGasProps = self.gam_i != None and self.cp_i != None and self.R_i != None 
        if not haveGasProps: return None 
        
        # Define this so we can avoid running extra calculation 
        HAVE_FLOW_PROPS = False
        
        # Handle handle if we have static properties
        if (self.Pi != None or self.Ti != None) and self.Mi == None:
            # Calculate Mach number from pressure ratio
            T_To = None if self.Ti == None else self.Ti/self.Toi  
            P_Po = None if self.Pi == None else self.Pi/self.Poi 
            max_Mach = 1 if AssumeSubsonic else None
            flow_props = GD.Isentropic_Flow(T_To=T_To, P_Po=P_Po, Gamma=self.gam_i, root_max_search=max_Mach)
            HAVE_FLOW_PROPS = True 
            # Only get Mach for now
            self.Mi = flow_props['Mach']
            
        
    
        # Check if we have area and mass flow rate without Mach number
        if self.Ai != None and self.mdot!=None and self.Mi == None: 
            # Also need mass flow rate
            if abs(self.mdot) < 1e-6:
                self.Mi = 0.0 
            else:
                try:
                    self.Mi = GD.Mach_at_mdot(self.mdot, self.Poi, self.Toi, self.Ai, Gamma=self.gam_i, R=self.R_i, gc=self.gc)
                except ValueError as e:
                    print(f'[Error]\t{self.StageName}_i at M0={self.Minf:.3f}: {e}')
                    # print(f'mdot={self.mdot:.6f}, Pt={self.Poi:.3f}, Tt={self.Toi:.3f}, A={self.Ai:.4f}, Gamma={self.gam_i:.3f}, R={self.R_i:.2f}')
                    self.Mi = None
           
            
        # Check if we have Mach number now
        if self.Mi != None:
            # Calculate if we havent already
            if not HAVE_FLOW_PROPS:
                flow_props = GD.Isentropic_Flow(Mach=self.Mi, Gamma=self.gam_i)
            # Calculate static properties
            self.Ti = flow_props['T_To'] * self.Toi
            self.Pi = flow_props['P_Po'] * self.Poi 
            
            # Calculate velocity
            self.Vi = self.Mi * np.sqrt(self.Ti*self.gam_i*self.R_i*self.gc)
        
            # Calculate area
            mdot_A = GD.mdot(self.Poi, self.Toi, 1, Mach=self.Mi, Gamma=self.gam_i, R=self.R_i, gc=self.gc)
            
            if abs(mdot_A) > 1e-8:
                self.Ai = None if self.mdot == None else self.mdot / mdot_A 
                self.Ai_mdota = self.mdot_ratio / mdot_A 
        
        return None
    
    def calcDetailProps_e(self, AssumeSubsonic=True):
        '''
        Calculates and updates all static, area, or Mach number values at the
        outlet of the stage.
        
        Currently inefficient but should work

        Returns
        -------
        None.

        '''

        
        # Always assume we have stagnation props, impossible to run outside code if not
        
        # Check if we have gas properties, if not return none
        haveGasProps = self.gam_e != None and self.cp_e != None and self.R_e != None 
        if not haveGasProps: return None 
        
        # Define this so we can avoid running extra calculation 
        HAVE_FLOW_PROPS = False
        
        # Handle handle if we have static properties
        if (self.Pe != None or self.Te != None) and self.Me == None:
            # Calculate Mach number from pressure ratio
            T_To = None if self.Te == None else self.Te/self.Toe  
            P_Po = None if self.Pe == None else self.Pe/self.Poe 
            max_Mach = 1 if AssumeSubsonic else None
            flow_props = GD.Isentropic_Flow(T_To=T_To, P_Po=P_Po, Gamma=self.gam_e, root_max_search=max_Mach)
            HAVE_FLOW_PROPS = True
            
            # Only set Mach for now, can save on extra compute time
            self.Me = flow_props['Mach']
            
           
        
        # Calc Mach if we have area and mass flow rate but not Mach
        if self.Ae != None and self.mdot!=None and self.Me == None: 
            # Also need mass flow rate
            if abs(self.mdot) < 1e-6:
                self.Me = 0.0 
            else:
                try:
                    self.Me = GD.Mach_at_mdot(self.mdot, self.Poe, self.Toe, self.Ae, Gamma=self.gam_e, R=self.R_e, gc=self.gc)
                except ValueError as e:
                    print(f'[Error]\t{self.StageName}_e at M0={self.Minf:.3f}: {e}')
                    self.Me = None
           
        # Check if we have Mach number now
        if self.Me != None:
            if not HAVE_FLOW_PROPS:
                flow_props = GD.Isentropic_Flow(Mach=self.Me, Gamma=self.gam_e)
            
            # Calculate static properties
            self.Te = flow_props['T_To'] * self.Toe
            self.Pe = flow_props['P_Po'] * self.Poe 
            
            # Calculate velocity
            self.Ve = self.Me * np.sqrt(self.Te*self.gam_e*self.R_e*self.gc)
        
            # Calculate area
            mdot_A = GD.mdot(self.Poe, self.Toe, 1, Mach=self.Me, Gamma=self.gam_e, R=self.R_e, gc=self.gc)
            self.Ae = None if (self.mdot == None) or abs(mdot_A) < 1e-6 else self.mdot / mdot_A 
            self.Ae_mdota = self.mdot_ratio / mdot_A 
        
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
            'R_i':self.R_i,
            'Ai': self.Ai,
            "Ai_mdota": self.Ai_mdota
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
            'R_e':self.R_e,
            'Ae': self.Ae,
            "Ae_mdota": self.Ae_mdota
            }
        BPRs = {}
        if isinstance(self.BPR, list):
            if len(self.BPR) == 1:
                BPRs["BPR"] = self.BPR[0]
            else:
                for i, alpha in enumerate(self.BPR):
                    BPRs[f"BPR_{i+1}"] = alpha
        else:
            BPRs = {'BPR': self.BPR}
            
        performance = {
            'mdot': self.mdot,
            'mdot_ratio': self.mdot_ratio,
            **BPRs,
            'Power':self.Power,
            'SpecPower': None if not hasattr(self,'specPower') else self.specPower,
            'ni':self.ni,
            'np': None if not hasattr(self,'np') else self.np,
            'pi': self.pi,
            'tau': self.tau,
            'f':None if not hasattr(self,'f') else self.f,
            }
        
        stageVals = {
            'Stage': self.StageName,
            'Units': self.units_labels,
            'inputs': inputs,
            'outputs': outputs,
            'performance': performance}
        return stageVals
    

    def __str__(self):
        return self.StageType
    
    def _getID(self):
        return self.StageIdentifier 

# =============================================================================
# Compressor Simple Stage Description
# =============================================================================

class Intake(Stage):
    def __init__(self, **kwargs):
        Stage.__init__(self, **kwargs)
        self.StageType = "Intake"
        self.UpdateInputs(StageName = "Intake", StageID="d", AssumeSubsonicExit=True, **kwargs)
        # NOTE: Ram efficiency ~= Isentropic Efficiency
        
    def UpdateInputs(self, **kwargs):
        Stage._update_inputs_base(self, **kwargs, )
        self.AssumeSubsonicExit = self.inputs['AssumeSubsonicExit']
        
        
    def calculate(self):
        # Always assume Pi/Pa and Ti/Ta are given (atmos conditions)
        # Cant always assume this, there is a ram affect and streamtubes
        #   occuring outside the intake face
        
        
        # Check if gas properties were given to intake
        if self.gam_i == None and self.gam_e == None:
            if self.cp_i == None and self.cp_e == None:
                print('[Warning]\t No gas constants inputted to Intake. Assuming Air')
        
                self.gam_i = 1.4  # Gamma value of air
                self.cp_i = 1.005 if self.units=='SI' else 0.240 # [kJ/kg*K] or [BTU/lbm*R]  cp of air
                self.R_i = self.gf*self.cp_i*(self.gam_i - 1)/self.gam_i        # [J/kg*K] or [ft*lbf/R*lbm]    R value of Air
                
        # print(f"{self.gf}, {self.cp_i}, {self.gam_i}")
        self.R_i = self.gf*self.cp_i*(self.gam_i - 1)/self.gam_i        # [J/kg*K] or [ft*lbf/R*lbm]    R value of Air
    
        # Ram Affects
        self.tau_r = 1 + (self.Minf**2)*(self.gam_i - 1)/2       # Totoal/static (Tt0/T0)
        self.pi_r  = self.tau_r**(self.gam_i / (self.gam_i - 1)) # Total/Static  (Pt0/P0)
        Pt0 = self.Pa*self.pi_r 
        Tt0 = self.Ta*self.tau_r
        self.pi = 1 if self.pi == None else self.pi # Pt2/Pt1
        
        # Run external conditions
        external = self.external_conditions(self.Minf, P0=self.Pa, T0=self.Ta, mdot=self.mdot, pi_dmax=self.pi, BPRs=self.BPR)
        self.D_additive = external.get('D_add', None)
        self.D_nacelle = external.get("D_nac", None)
        self.mdot = external.get('mdot', self.mdot)
        self.Mi = external.get('M1', self.Mi)
        self.A0 = external.get('A0', None)
        self.eta_r = external['eta_r'] # Pt1/Pt0 
        
        # Correct the diffusor pressure ratio
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
                
        # Now we should have mach num no matter what
        # and the static props (atm props)
        # Find inlet stag props        
        self.Toi = self.Ta * self.tau_r
        self.Poi = self.Pa * self.pi_r * self.eta_r # At station 1 (inlet plane)
        
        self.Pi = self.Poi / GD.Po_P_ratio(self.Mi, self.gam_i) # Corresponds to P1 (static pressure at inlet face, not freestream)
        self.Ti = self.Toi / GD.To_T_ratio(self.Mi, self.gam_i) # Corresponds to T1 (static pressure at inlet face, not freestream
        
        # Calculate detailed values
        self.calcDetailProps_i()
        
        # Find outlet stag props
        self.tau = 1 if self.tau == None else self.tau
        self.Toe = self.Toi * self.tau # Typically Adiabatic (tau =1)
        self.Poe = self.pi * self.Poi # Using only pi_dmax because it is the pressure loss within diffuser. 
                                      # Poi already contains pressure losses from external phenomena
        
        
        self.gam_e = self.gam_i 
        self.cp_e = self.cp_i 
        self.R_e = self.R_i
        
        # Calculate detailed values
        self.calcDetailProps_e(self.AssumeSubsonicExit) # Assume subsonic at diffuser exit?
        # Done
        
        
    def external_conditions(self, Minf, **kwargs):
        '''
        Calculates the Ram Efficiency eta_r from the Mach number since
        the equation used will vary depending on freestream Mach. Using
        the milspec experimentally based case structure.
        
        Can be replaced by an outside function that uses the same inputs/outputs
        
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
                'A0'    - Stream tube capture area
        
        '''
        if Minf < 0:
            EngineErrors.MissingValue("Minf is Negative", self.StageName)
        elif Minf < 1: # or self.IS_IDEAL:
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
        outs['performance']['D_nac'] = self.D_nacelle
        return outs
        
# =============================================================================
# Compressor Simple Stage Description
# =============================================================================
class Compressor(Stage):
    def __init__(self, **kwargs):
        Stage.__init__(self, **kwargs)
        self.StageType = "Compressor"
        self.UpdateInputs(StageName = "Compressor", StageID="c", **kwargs)
        # Adding PR and BPR
        # self.r = kwargs.get('rc') # Pressure Ratio of stage (pi)
        
   
    def UpdateInputs(self, **kwargs):
        Stage._update_inputs_base(self, **kwargs)
        self.BPR = self.inputs.get('BPR', None) # Bypass Ratio (α): bypass mass flow (air)/mass flow through core (air)
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
        
        if self.mdot == None:
            self.specPower = self.mdot_ratio*self.cp_i*(self.Toe-self.Toi)
        else:
            self.Power = self.mdot*self.cp_i*(self.Toe-self.Toi)
            self.specPower = self.Power/self.mdot
            
        self.tau = self.Toe/self.Toi if self.tau == None else self.tau
        self.ni = self.calculate_ni_c(self.np, self.gam_i) if self.ni == None else self.ni 
        
        self.gam_e = self.gam_i 
        self.cp_e = self.cp_i 
        self.R_e = self.R_i
        
        # Calculate detailed props
        self.calcDetailProps_i()
        self.calcDetailProps_e()
        
        # Done
   
        
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
        self.StageType = "Combustor"
        self.UpdateInputs(StageName = "Combustor", StageID="b", **kwargs)
        
    
    def UpdateInputs(self, **kwargs):
        self._update_inputs_base(**kwargs) 
        self.dTo = self.inputs.get('dTb')
        # Commenting out to use the pressure ratio pi instead
        # self.dPo = kwargs.get('dPb_dec', 0) # the pressure loss within the compressor as a decimal (0.05 = 5% loss)
        self.f  = self.inputs.get('f') # actual/real fuel-air-ratio 
        self.Q  = self.inputs.get('Q_fuel')
        self.pi_off = self.inputs.get('pi_off', None)
        # self.ni = self.inputs.get('nb', 1) # Combustor efficiency
        # self.IS_IDEAL = self.inputs.get('IS_IDEAL') # Shouldnt need? Defined in stage
                
    def calculate(self):
        # Assuming we have the Toi and Poi from compressor/prev stage
        # We need to have the exit 
        # Assuming that we have initial gas properties and exit gas properties were set from engine
        if self.R_e == None:
            try:
                self.R_e = self.gf*self.cp_e*(self.gam_e - 1)/self.gam_e        # [J/kg*K] or [ft*lbf/R*lbm]    R value of Air
            except:
                self.R_e = self.R_i 
                self.cp_e = self.cp_i 
                self.gam_e = self.gam_i 
                print("[Warning]\t Combustor gas exit properties set to inputs.")
        
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
        self.tau = self.Toe / self.Toi
        
        # Check if combustor is on/off (tau!=1)
        if abs(1-self.tau) < 1e-6:
            # burner is off 
            self.pi = 1 if self.pi == None else self.pi 
            self.pi = self.pi_off if self.pi_off != None else self.pi
        
        self.Poe = self.Poi*self.pi 
        self.dTo = self.Toe - self.Toi # will use later for f calcs
        
        # if self.f == None:
        if self.Q != None: 
            # Assuming non-ideal, will calculate f and use in mass fuel flow
            self.f = (self.cp_e*self.Toe - self.cp_i*self.Toi) / (self.ni*self.Q - self.cp_e*self.Toe)
            if self.IS_IDEAL:
                self.f = self.cp_i*(self.Toe-self.Toi) / self.Q
            
        # Note: f = mdot_fuel / mdot_core_air
        if not self.IS_IDEAL:
            # Only change mass flow rate if not ideal case
            if self.mdot != None and self.f != None:
                self.mdot += self.f*self.mdot
                self.mdot_ratio += self.f*self.mdot_ratio
            elif self.mdot == None and self.f != None:
                self.mdot_ratio += self.f*self.mdot_ratio
                
        # Calculate details if available
        self.calcDetailProps_i(True)
        self.calcDetailProps_e(True)
            
    def extraOutputs(self,form):
        if self.f != None:
            print('\t f = {}'.format(form).format(self.f))
        
    def StageValues(self):
        outs = Stage.StageValues(self)
        outs['performance']['dTo'] = self.dTo
        return outs
            
# =============================================================================
# Turbine Simple Stage Description
# =============================================================================
class Turbine(Stage):
    def __init__(self, Comp_to_power, **kwargs):
        Stage.__init__(self, **kwargs)
        self.StageType = "Turbine"
        self.Compressor = Comp_to_power
        self.UpdateInputs(StageName = "Turbine", StageID='t', **kwargs)


    def UpdateInputs(self, **kwargs):
        self._update_inputs_base(**kwargs) 
        self.np = self.inputs.get('np') # Polytropic efficiency
        self.nm = self.inputs.get('nm',1) # Mechanical efficiency
         # Could be list
        # Will have inlet temp, compressor power
        #self.r  = kwargs.get('rt') # Add for later, not used now
        # this will be for generators or when turbine pressure ratio is specified
        

    def calculate(self):
        if self.mdot != None:
            if type(self.Compressor) == list:
                com_power = 0
                for i in range(0,len(self.Compressor)):
                    com_power += self.Compressor[i].Power
                self.Power = com_power/self.nm 
            else:
                self.Power = self.Compressor.Power/self.nm 
            # Calculate exit temp
            self.Toe = self.Toi - self.Power/(self.mdot*self.cp_i)
            self.specPower = self.Power/self.mdot
        else:
            # No mdot is given, need to power balance based on mass-flow ratios instead
            if type(self.Compressor) == list:
                com_power = 0
                for i in range(0,len(self.Compressor)):
                    com_power += self.Compressor[i].specPower
                self.specPower = com_power/self.nm 
               
            else:
                self.specPower = self.Compressor.specPower/self.nm 
            # Calculate exit temp
            '''IS THIS ASSUMING CONSTANT CP? 20251010'''
            self.Toe = self.Toi - self.specPower/(self.mdot_ratio*self.cp_i)
        
            
        if self.np == None:
            if self.pi != None:
                # Calculate np
                self.np = np.log(1- self.ni*(1 - self.pi**((self.gam_i-1)/self.gam_i)))
                self.np /= np.log(self.pi)*(self.gam_i-1)/self.gam_i
            else:
                print('[Warning]\t Insufficient parameters given to turbine')
                print('[Warning]\t Continuing assuming polytropic efficiency = 1')
                self.np = 1
                
        m_frac = self.np*(self.gam_i-1)/self.gam_i
        self.Poe = self.Poi*(1- (self.Toi-self.Toe)/self.Toi )**(1/m_frac)
        
        self.pi = self.Poe/self.Poi 
        self.tau = self.Toe/self.Toi 
        
        self.gam_e = self.gam_i 
        self.cp_e = self.cp_i 
        self.R_e = self.R_i
        
        self.ni = (1 - self.pi**((self.gam_e-1)*self.np/self.gam_e)) / (1 - self.pi**((self.gam_e-1)/self.gam_e))
    
        # Calculate details if available
        self.calcDetailProps_i(True)
        self.calcDetailProps_e(True)
    # Done
        
    
    
# =============================================================================
# Mizer Simple Stage Description    
# =============================================================================
class Mixer(Stage):
    def __init__(self, CoreMixStage, BypassMixStage, **kwargs):
        Stage.__init__(self, **kwargs)
        self.StageType = "Mixer"
        self.CoreMix = CoreMixStage
        self.BypassMix = BypassMixStage
        self.UpdateInputs(StageName = "Mixer",StageID='M', **kwargs)

    def UpdateInputs(self, **kwargs):
        self._update_inputs_base(**kwargs) 
        self.BPR_f = self.inputs.get('BPR_f')
        self.MIX_GAS_PROPERTIES = self.inputs.get('MIX_GAS_PROPERTIES', True)
        
    def calculate(self):
          
        # FOR PROJECT 1:
        # We have M6 (core exit Mach) and pi_M (pressure ratio for mixer)
        # Assumes that gamma for bypass is the same for air 
        
        self.mdot_ratio = self.CoreMix.mdot_ratio + self.BypassMix.mdot_ratio
        if self.mdot != None:
            self.mdot = self.CoreMix.mdot + self.BypassMix.mdot
            
            
        if self.IS_IDEAL:
            # Calculate bypass ratio from core flow ratio
            raise Warning("[Warn]\t Ideal Mixer currently only support a single bypass.")
            # NOTE: Only applicable with one bypass
            alpha = 1/self.CoreMix.mdot_ratio - 1
            self.tau = (1 + alpha * (self.BypassMix.Toe/self.CoreMix.Toe)) / (1 + alpha)
            self.pi = 1 if self.pi == None else self.pi
            
            self.Toe = self.tau * self.CoreMix.Toe 
            self.Poe = self.pi * self.CoreMix.Poe
            
            self.gam_e =  self.CoreMix.gam_e
            self.cp_e = self.CoreMix.cp_e
            self.R_e = self.CoreMix.R_e 
            # print(f"Mix Re = {self.R_e}")
            
            # Calculate details if available
            self.calcDetailProps_i(True)
            self.calcDetailProps_e(True)
        
            return None
        
        
            
        
        # Check to ensure we have required mach numbers
        if self.CoreMix.Me == None:
            if abs(self.mdot_ratio) < 1e-8:
                self.CoreMix.Me = 0.0 # Try this 
            else:
                raise EngineErrors.MissingValue('Missing Coreflow Mach Number',self.StageName)
          
        # M_16 calculation
        # Gonna just rename values to make double checking equation easier
        M_6 = self.CoreMix.Me 
        M_16 = self.BypassMix.Me 
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
        
        self.gam_i = self.CoreMix.gam_e
        self.R_i = self.CoreMix.R_e
        self.cp_i = self.CoreMix.cp_e
        self.Poi = self.CoreMix.Poe 
        self.Toi = self.CoreMix.Toe
        
        
        # Check if Mach number at bypass was calculated:
        if M_16 == None:
            M_16 = np.sqrt( (2/(gam_16 - 1)) * (( (P_t16/P_t6)*(1 + (M_6**2)*(gam_6 - 1)/2)**(gam_6/(gam_6-1))  )**((gam_16 - 1)/gam_16)  - 1)  )
            self.BypassMix.Me = M_16 
        
        # Check if mdot from bypass or M_16 == 0
        NO_BYPASS = False
        if abs(self.BypassMix.mdot_ratio) < 1e-5:
            NO_BYPASS = True
        if abs(M_16) < 1e-5:
            NO_BYPASS = True 
        
        # BPR of traditional mixer = mdot_bp / (mdot_core + mdot_fuel), 
        # so BPR_f = mdot_rat_bp / mdot_rat_core 
        if abs(self.BypassMix.mdot_ratio) < 1e-8 or abs(self.CoreMix.mdot_ratio) < 1e-8:
            self.BPR_f = 0 
        else:
            self.BPR_f = self.BypassMix.mdot_ratio / self.CoreMix.mdot_ratio
        
        # gas props will equal station 6 (core) if BPR_F = 0
        self.cp_6a = (cp_6 + self.BPR_f*cp_16) / (1 + self.BPR_f)
        self.R_6a = (R_6 + self.BPR_f*R_16) / (1 + self.BPR_f)
        self.gam_6a = self.cp_6a*self.gf / (self.cp_6a*self.gf - self.R_6a)
        
        # Will eqal 1 if BPR_f = 0
        self.tau = (cp_6 / self.cp_6a) * (1 + self.BPR_f * (cp_16/cp_6)*(self.BypassMix.Toe/self.CoreMix.Toe)) / (1 + self.BPR_f)
        
        # if self.pi==None:
            # Use ideal pi 
        # A16_A6 = self.BPR_f*(P_t6/P_t16)*np.sqrt(T_t16/T_t6)*(GD.MassFlowParam_norm(M_6 ,gam_6)/(GD.MassFlowParam_norm(M_16,gam_16)))
        
        # Calculate A16/A6
        if self.CoreMix.Ae == None or self.BypassMix.Ae == None:
            num = gam_16*R_6 * (1 + (M_16**2)*(gam_16-1)/2)
            den = gam_6*R_16 * (1 + (M_6**2)*(gam_6-1)/2)
            if NO_BYPASS:
                A16_A6 = 0.0
            else:
                A16_A6 = self.BPR_f*np.sqrt(T_t16/T_t6)*(M_6/M_16) / np.sqrt(num/den)
        else:
            A16_A6 = self.BypassMix.Ae / self.CoreMix.Ae
        self.A16_A6 = A16_A6
        # Calculate phi(M6A, gam_6a):
        num = 1 + self.BPR_f 
        if NO_BYPASS:
            den = 1/np.sqrt(GD.Rayleigh_phi_MS(M_6,gam_6)) + 0 # This avoids divide by 0 error
        else:
            den = 1/np.sqrt(GD.Rayleigh_phi_MS(M_6,gam_6)) + self.BPR_f*np.sqrt(R_16*gam_6*(T_t16/T_t6) / (R_6*gam_16*GD.Rayleigh_phi_MS(M_16,gam_16)))
        
        phi = ((num/den)**2)* (self.R_6a*gam_6*self.tau)/(R_6*self.gam_6a)
        # Calculate M_6A, should be = M6 if BPR_f=0
        M_6A = np.sqrt(2*phi / (1 - 2*self.gam_6a*phi + np.sqrt(1 - 2*(self.gam_6a+1)*phi)))
       
        self.pi_ideal = ((1+self.BPR_f)*np.sqrt(self.tau)/(1+ A16_A6))*(GD.MassFlowParam_norm(M_6 ,gam_6)*np.sqrt(self.gc/R_6) / (GD.MassFlowParam_norm(M_6A,self.gam_6a)*np.sqrt(self.gc/self.R_6a)))
        
        self.pi = 1 if self.pi == None else self.pi
        self.pi_M = self.pi_ideal * self.pi
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
            
        # Calculate details if available
        self.calcDetailProps_i(True)
        self.calcDetailProps_e(True)
        ## NEED TO GET ALL VALUES FOR THIS
        # alpha_m = (eta_m * (1 + f) * (tau_t / tau_c) * (1 - (pi_f/(pi_c*pi_b))**((gamma_t - 1)*e_t/gamma_t)) - (tau_c - 1)) / (tau_f - 1)
        
        return None
    
    def StageValues(self):
        outs = Stage.StageValues(self)
        if not self.IS_IDEAL:
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
        
        self.StageType = "Nozzle"
        self.nozzle_type = nozzle_type # 'C' for Converging, 'CD' for Conv-Div
        self.UpdateInputs(StageName = "Nozzle",StageID='n', **kwargs)


    def UpdateInputs(self, **kwargs):
        self._update_inputs_base(**kwargs) 
        self.CV = self.inputs.get("CA", None) 
        self.CA = self.inputs.get("CA", None) 
        self.CD = self.inputs.get("CD", None)
        self.At = self.inputs.get("At", None)
        self.At_max = self.inputs.get("At_max", None)
        self.At_min = self.inputs.get("At_min", None)
        self.Ae_max = self.inputs.get("Ae_max", None)
        self.Ae_min = self.inputs.get("Ae_min", None)
        
        
    def calculate(self):
        # Note: Assume all pressure loss occurs after throat for a CD Nozzle
        self.IsChoked, Pc = self._check_nozzle_choke()
        if self.IsChoked:
            # Nozzle is choked
            # Check nozzle type
            if self.nozzle_type == 'C':
                # For Converging nozzle
                self.Pe = Pc 
            else:
                # CD Nozzle
                if self.Pe == None: 
                    # No Value provided
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
            raise ValueError(f"[Error] Missing eta_n or pi_n in {self.StageName}")
        
        # Stag props at exit
        self.Poe = self.pi * self.Poi
        self.Toe = self.Toi
        
        try:
            Me_mdot = GD.Mach_at_mdot(self.mdot, self.Poe, self.Toe, self.Ae_max, self.gam_i, self.R_i, self.gc, forceSupersonic=self.IsChoked)
        except:
            Me_mdot = 10000 
    
        Me_max = GD.Mach_at_PR(self.Poe/self.Pe, Gamma = self.gam_i)
        self.Me = min(Me_max, Me_mdot) 
        self.Pe = self.Poe/GD.Po_P_ratio(self.Me, self.gam_i)
        self.Te = self.Toi*(1-self.ni*(1-(self.Pe/self.Poi)**((self.gam_i-1)/self.gam_i)))
        self.Ve = self.Me*np.sqrt(self.gam_i*self.R_i*self.Te*self.gc)
        
        
        
        # self.Poe = self.Pe * (1 + (self.gam_i -1)*(self.Me**2)/2)**(self.gam_i/(self.gam_i -1))
        self.tau = self.Toe/self.Toi
        
        self.gam_e = self.gam_i 
        self.cp_e = self.cp_i 
        self.R_e = self.R_i 
        
        
        # Calculate details if available
        self.calcDetailProps_i()
        self.calcDetailProps_e()
        
        # Calculate Nozzle Drag
        self._calculate_drag()
        
        self._zeroPropsIfNoFlow()
        
    def _check_nozzle_choke(self):
        '''
        Checks if the nozzle is choked. 
        Process (mainly for CD nozzle):
            - Check if mass-flow rate is provided
                - If so, check if throat area is provided and if mdot there = Mach 1
                - If throat area not provided, calculate it and check if within bounds for throat areas
                
    

        Returns
        -------
        IsChoked
            DESCRIPTION.
        Pc : TYPE
            DESCRIPTION.

        '''
        isChoked = False
        # First check if mass flow is provided:
        # if self.mdot != None:
        #     # Calculate necessary throat area
        #     Astar = GD.A_throat(self.mdot, self.Poi, self.Toi, self.gam_i, self.R_i, self.gc)
        #     # Check if throat area provided:
        #     if self.At != None:
        #         # Check if At ~ Astar
        #         if abs(self.Astar - self.At) < 1e-5:
        #             # Is choked and areas match
        #             IsChoked = True
        #         elif self.Astar > self.At:
        #             # Check if necessary area is larger (would reduce mass flow and cause surge)
        #             raise Warning(f"[Warning]\t {self.StageName} choke reducing mass flow.")
        #         else: # self.Astar < self.At:
        #             # Not choked
        #             IsChoked = False
            
        #     else: # throat area not provided
        #         # Check if min and max are provided
        #         if self.At_max != None and self.At_min !=None:
        #             # Check if between values
        #             if self.At_max > Astar and self.At_min < Astar:
        #                 # It is, so set throat area to Astar
        #                 self.At = Astar 
        #             else:
        #                 raise Warning(f"[Warning]\t {self.StageName} - Required throat area outside bounds of minmumt and maximum At\n\t" + 
        #                               f"A*={Astar:.4f}, At_min={self.At_min:.4f}, At_max={self.At_max:.4f}")
                    
                
        
        
        # Check if choked
        Tc = self.Toi*(2/(self.gam_i+1))
        ni = 1 if self.ni == None else self.ni # TEMPORARY BC CURRENTLY NOT SUPPORTED TO NOT HAVE IT AN INPUT
        Pc = self.Poi*(1 - (1/ni)*(1-Tc/self.Toi))**(self.gam_i/(self.gam_i-1))
        
        P_rat = self.Poi/self.Pa
        P_crit = self.Poi/Pc
        if P_rat > P_crit:
            IsChoked = True
        else:
            IsChoked = False
        
        return IsChoked, Pc
    
    def _calculate_drag(self):
        if not self.IS_IDEAL:
            if self.ni != None:
                self.CV = np.sqrt(self.ni)
            else:
                self.CV = 1
            self.CA = 1 if self.CA == None else self.CA
            self.CD = 1 if self.CD == None else self.CD
            self.Cfg = self.CD*self.CV*np.sqrt((1 - ((self.Pe/self.Pa)/(self.NPR*self.CD**2))**((self.gam_e-1)/self.gam_e))/(1-self.NPR**(-(self.gam_e-1)/self.gam_e))) * \
                (self.CA + (((self.gam_e - 1)/(2*self.gam_e))*(1 - self.Pa/self.Pe))/((self.pi*self.NPR*self.Pa/self.Pe)**((self.gam_e-1)/self.gam_e) -1))
            if self.mdot != None:
                self.Fg_ideal = self.mdot*self.Ve + self.Ae*(self.Pe - self.Pa)
                self.Drag = self.Fg_ideal*(1-self.Cfg)
            else: 
                self.Fg_ideal = None
                self.Drag = None
        
            self.Fg_ideal_mdota = self.mdot_ratio *self.Ve + self.Ae_mdota *(self.Pe - self.Pa)
            self.Drag_mdota = self.Fg_ideal_mdota*(1-self.Cfg)
        else:
            self.Cfg = None
            self.Fg_ideal_mdota = None
            self.Drag_mdota = None
            self.Fg_ideal = None
            self.Drag = None
        return None
        
    def _zeroPropsIfNoFlow(self):
        if abs(self.mdot_ratio) < 1e-8:
            self.Me = 0 
            self.Ve = 0 
            self.Drag = 0 
            self.Drag_mdota = 0 
            self.Pe = self.Pa 
            self.Te = self.Ta 
            self.Ae = 0
            self.Ae_mdota = 0
            
    
    def StageValues(self):
        outs = Stage.StageValues(self)
        outs['performance']['C_V'] = self.CV
        outs['performance']['C_A'] = self.CA
        outs['performance']['C_D'] = self.CD
        outs['performance']['C_fg'] = self.Cfg
        outs['performance']['Fg_ideal'] = self.Fg_ideal
        outs['performance']['D_noz'] = self.Drag
        outs['performance']['Fg_ideal/mdot0'] = self.Fg_ideal_mdota
        outs['performance']['D_noz/mdot0'] = self.Drag_mdota
        
        return outs
        
        
   

class Duct(Stage):
    def __init__(self, **kwargs): 
        # Adiabatic duct without friction currently
        Stage.__init__(self, **kwargs)
        self.StageType = "Duct"
        self.UpdateInputs(StageName = 'Duct',StageID='bp', **kwargs)
        
        
    def UpdateInputs(self,**kwargs):
        self._update_inputs_base(**kwargs) 
    
    def calculate(self):
        self.tau = 1 if self.tau == None else self.tau
        self.pi = 1 if self.pi == None else self.pi
        
        self.Poe = self.pi * self.Poi 
        self.Toe = self.tau * self.Toi
         
        
        self.gam_e = self.gam_i 
        self.cp_e = self.cp_i 
        self.R_e = self.R_i
        
        # Calculate details if available
        self.calcDetailProps_i()
        self.calcDetailProps_e()
        
        
        
# =============================================================================
# Testing Main
# =============================================================================

if __name__ == "__main__":
    turb = Compressor() 
    listTest = [[1, [1,2]], [2, None], [3,None], [4, [1,2,3]]]
    print(len(listTest[0]))
    print(type(turb))
    stageType = str(turb)
    print(stageType)
