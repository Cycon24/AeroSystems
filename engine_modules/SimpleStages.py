# -*- coding: utf-8 -*-
"""
Created on Thu Aug 24 09:45:15 2023

@author: cycon

20240630: 
    Revised and seperated into different files, renamed to SimpleStages.py
    Moved engine definition portion into EngineCompiler.py
    Moved complex stage solving into ComplexStages.py
"""
import numpy as np
import EngineErrors as EngineErrors

'''
Coding stuff
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
        # Constants
        self.R = 287 # J/kg*K    R value of Air
        self.gam_a = 1.4        # Gamma value of air
        self.gam_g = 4/3        # Gamma value of gas (combustion products)
        self.cp_a = 1.005 # kJ/kg*K     cp of air
        self.cp_g = 1.148 # kJ/kg*K     cp of gas (combustion products)
        self.g = 9.81 # m/s^2   # gravitational constant
        
        
        # General properties that could be used by all stages 
        # so all components know the atm conditions
        self.Ta  = kwargs.get('Ta') 
        self.Pa  = kwargs.get('Pa')
        self.Vinf = kwargs.get('Vinf')
        self.Minf = kwargs.get('Minf')
        
        # All Stage Initial Conditions (not every one is needed)
        self.Toi = kwargs.get('Toi')
        self.Poi = kwargs.get('Poi')
        self.Ti  = kwargs.get('Ti')
        self.Pi  = kwargs.get('Pi')
        self.Mi  = kwargs.get('Mi')
        self.Vi  = kwargs.get('Vi')
        
        # All Stage Exit Conditions
        self.Toe = kwargs.get('Toe')
        self.Poe = kwargs.get('Poe')
        self.Te  = kwargs.get('Te')
        self.Pe  = kwargs.get('Pe')
        self.Me  = kwargs.get('Me')
        self.Ve  = kwargs.get('Ve')
        
        # Other properties
        self.m_dot = kwargs.get('m_dot') # Stays constant through component
        self.mdot_ratio = 1 # used to track mass flow ratio through sections
        self.BPR = kwargs.get('BPR',1)
        
        # Isentropic Efficiency (basically stage efficiency)
        self.ni   = kwargs.get('ni', 1) # Isentropic efficiency
       
        # Other properties
        self.StageName = ""
        self.Power = None       # W    - Power
        self.SpecPower = None  #  J/kg - Specific Power
        
    def forward(self, next_Stage):
        '''
        Propogates the properties of this stage into the next stage.

        Parameters
        ----------
        next_Stage : Object
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
    
    def printOutputs(self):
        '''
        Print the outputs of a stage to the console for tracking values and debugging.

        Returns
        -------
        None.

        '''
        form ='{:9.3f}'
        print('Stage: ', self.StageName)
        if self.Toe != None:
            print('\t Toe = {} K'.format(form).format(self.Toe))
        if self.Poe != None:
            print('\t Poe = {} Pa'.format(form).format(self.Poe))
        if self.Te != None:
            print('\t Te  = {} K'.format(form).format(self.Te))
        if self.Pe != None:
            print('\t Pe  = {} Pa'.format(form).format(self.Pe))
        if self.m_dot != None:
            print('\tmdot = {} kg/s'.format(form).format(self.m_dot))
        if self.Me != None:
            print('\t Me  = {}'.format(form).format(self.Me))
        if self.Ve != None:
            print('\t Ve  = {} m/s'.format(form).format(self.Ve))
        if self.Power != None:
            print('\t Pow = {} W'.format(form).format(self.Power))
        if self.SpecPower != None:
            print('\t Specific Pow = {} J/kg'.format(form).format(self.specPower))
        self.extraOutputs()
    
    def extraOutputs(self):
        # Overwrite this and put any extra outputs here within individual stages
        return None

# =============================================================================
# Compressor Simple Stage Description
# =============================================================================

class Intake(Stage):
    def __init__(self, **kwargs):
        Stage.__init__(self, **kwargs)
        self.StageName = "Intake"
        # NOTE: Ram efficiency ~= Isentropic Efficiency
        
    def calculate(self):
        # Always assume Pi/Pa and Ti/Ta are given (atmos conditions)
        self.Pi = self.Pa
        self.Ti = self.Ta
        self.Vi = self.Vinf
        self.Mi = self.Minf
        # If no vel or mach num inputted, assume stationary
        if self.Mi == None:
            if self.Vi == None:
                self.Mi = 0
            else:
                self.Mi = self.Vi/np.sqrt(self.gam_a*self.R*self.Ti)
        else:
            if self.Vi == None:
                self.Vi = self.Mi*np.sqrt(self.gam_a*self.R*self.Ti)
        
        # Now we should have mach num no matter what
        # and the static props (atm props)
        self.Toe = self.Ti * (1 + (self.gam_a-1)*(self.Mi**2)/2)
        self.Poe = self.Pi * (1 + self.ni*(self.Mi**2)*(self.gam_a-1)/2)**(self.gam_a/(self.gam_a-1))
        
        

# =============================================================================
# Compressor Simple Stage Description
# =============================================================================
class Compressor(Stage):
    def __init__(self, **kwargs):
        Stage.__init__(self, **kwargs)
        self.StageName = "Compressor"
        # Adding PR and BPR
        self.r = kwargs.get('rc') # Pressure Ratio of stage
        self.BPR = kwargs.get('BPR', 1) # Bypass Ratio: total mass flow (air)/mass flow through core
        self.np = kwargs.get('np') # Polytropic efficiency
        self.mdot_ratio = 1 # Starts as 1 for fan, will be updated by prior comp if
                            # different from 1 from forward section
    def calculate(self):
        # Should always have input To and Po, need to calculate power
        # and output To and Po. r will always be given, BPR will affect output 
        # to next stage
        if self.r == None:
            raise EngineErrors.MissingValue('R-Press. Ratio','Compressor')
        elif self.np == None:
            self.np = ((self.gam_a-1)/self.gam_a)*np.log(self.r) / \
                        np.log( (self.r**((self.gam_a-1)/self.gam_a) - 1)/self.ni + 1)
        
        n_frac =  (self.gam_a-1)/(self.gam_a*self.np)
        self.Toe = self.Toi + self.Toi*(self.r**n_frac - 1)
        self.Poe = self.r*self.Poi
        
        if self.m_dot == None:
            self.specPower = self.mdot_ratio*self.cp_a*(self.Toe-self.Toi)
        else:
            self.Power = self.m_dot*self.cp_a*(self.Toe-self.Toi)
        # Done
        
        
    def forward(self, next_Stage_hot, next_Stage_cold=None):
        next_Stage_hot.Toi = self.Toe
        next_Stage_hot.Poi = self.Poe
        next_Stage_hot.Ti  = self.Te
        next_Stage_hot.Pi  = self.Pe
        next_Stage_hot.Mi  = self.Me
        next_Stage_hot.Vi  = self.Ve
        next_Stage_hot.mdot_ratio = self.mdot_ratio
        
        if next_Stage_cold == None:
            next_Stage_hot.m_dot = self.m_dot
        else:
            if self.BPR == None:
                raise EngineErrors.MissingValue('BPR','Compressor')
            else:
                if self.m_dot != None:
                    m_dot_h = self.m_dot/(1 + self.BPR)
                    m_dot_c = self.m_dot - m_dot_h
                    
                    next_Stage_hot.m_dot = m_dot_h
                    next_Stage_cold.m_dot = m_dot_c
                else:
                    # No inputted mdot
                    mdot_ratio_h = 1/(1+self.BPR)
                    
                    next_Stage_hot.m_dot = None
                    next_Stage_hot.mdot_ratio = mdot_ratio_h
                    next_Stage_cold.m_dot = None
                    # Dont need to send mdot ratio to cold section
                    
                next_Stage_cold.Toi = self.Toe
                next_Stage_cold.Poi = self.Poe
                next_Stage_cold.Ti  = self.Te
                next_Stage_cold.Pi  = self.Pe
                next_Stage_cold.Mi  = self.Me
                next_Stage_cold.Vi  = self.Ve
    
        
    def calculate_nc(self, np, gamma=1.4):
        '''

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
        nc = ( self.r**((self.gam_a-1)/self.gam_a) - 1 ) / ( self.r**((self.gam_a-1)/(self.gam_a*np)) - 1 )
        return nc
        
  

    
# =============================================================================
# Combustor Simple Stage Description  
# =============================================================================
class Combustor(Stage):
    def __init__(self, **kwargs):
        Stage.__init__(self, **kwargs)
        self.StageName = "Combustor"
        self.dTo = kwargs.get('dTb')
        self.dPo = kwargs.get('dPb_dec', 0) # the pressure loss within the compressor as a decimal (0.05 = 5% loss)
        self.f  = kwargs.get('f')
        self.Q  = kwargs.get('Q_fuel')
        self.nb = kwargs.get('nb', 1) # Combustor efficiency
        
    def calculate(self):
        # Assuming we have the Ti and Pi from compressor/prev stage
        # We need to have the exit 
        if self.Toe == None: 
            # No Turbine inlet temp given
            if self.dTo == None: 
                # No combustor increase temp given
                if self.f == None and self.Q == None:
                    # No air-fuel ratio given, cant calculate temps
                    raise EngineErrors.MissingValue('Toe, dTo, or f&Q','Combustor')
                else: 
                    # We have f and Q to calculate exit temp
                    f_ideal = self.f*self.nb # inputted f would be actual
                    self.Toe = (f_ideal*self.Q + self.cp_a*self.Toi)/(self.cpg(1+f_ideal))
            else:
                # We dont have exit temp, but do have temp increase
                self.Toe = self.Toi + self.dTo
         # else: Dont need to use since we have what we need
             # We have turbine inlet temp (Te)
             
        self.Poe = self.Poi*(1-self.dPo)
        self.dTo = self.Toe - self.Toi # will use later for f calcs
        
        if self.f == None:
            if self.Q != None: 
                # Assuming non-ideal, will calculate f and use in mass fuel flow
                self.f = (self.cp_g*self.Toe - self.cp_a*self.Toi) / (self.ni*(self.Q - self.cp_g*self.Toe))
                
        if self.m_dot != None and self.f != None:
            self.m_dot += self.f*self.m_dot
        elif self.m_dot == None and self.f != None:
            self.mdot_ratio = (1+self.f)/(self.BPR + 1)
        
            
# =============================================================================
# Turbine Simple Stage Description
# =============================================================================
class Turbine(Stage):
    def __init__(self, Comp_to_power, **kwargs):
        Stage.__init__(self, **kwargs)
        self.StageName = "Turbine"
        self.np = kwargs.get('np') # Polytropic efficiency
        self.nm = kwargs.get('nm',1)
        self.Compressor = Comp_to_power # Could be list
        # Will have inlet temp, compressor power
        self.r  = kwargs.get('rt') # Add for later, not used now
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
            self.Toe = self.Toi - self.Power/(self.m_dot*self.cp_g)
        else:
            # No m_dot is given, need to power balance based on 
            # BPR ratios instead
            if type(self.Compressor) == list:
                com_power = 0
                for i in range(0,len(self.Compressor)):
                    com_power += self.Compressor[i].Power_ratio
                self.specPower = com_power/self.nm 
            else:
                self.specPower = self.Compressor.Power_ratio/self.nm 
            # Calculate exit temp
            self.Toe = self.Toi - self.specPower/(self.mdot_ratio*self.cp_g)
        
            
        if self.np == None:
            if self.r != None:
                # Calculate np
                self.np = np.log(1- self.ni*(1 - self.r**((self.gam_g-1)/self.gam_g)))
                self.np /= np.log(self.r)*(self.gam_g-1)/self.gam_g
            else:
                print('Warning: insufficient parameters given to turbine')
                print('Continuing assuming polytropic efficiency = 1')
                self.np = 1
                
        m_frac = self.np*(self.gam_g-1)/self.gam_g
        self.Poe = self.Poi*(1- (self.Toi-self.Toe)/self.Toi )**(1/m_frac)
    
        
# =============================================================================
# Nozzle Simple Stage Description
# =============================================================================
class Nozzle(Stage):
    def __init__(self, air_type='hot', **kwargs):
        Stage.__init__(self, **kwargs)
        self.StageName = "Nozzle"
        if air_type == 'hot':
            self.gam = self.gam_g
        else:
            self.gam = self.gam_a
        
    def calculate(self):
        # Check if choked
        Tc = self.Toi*(2/(self.gam_g+1))
        Pc = self.Poi*(1 - (1/self.ni)*(1-Tc/self.Toi))**(self.gam/(self.gam-1))
        
        P_rat = self.Poi/self.Pa
        P_crit = self.Poi/Pc
        if P_rat > P_crit:
            # Nozzle is choked
            if self.Pe == None:
                self.Pe = Pc
            else:
                # We are given exit pressure. For now assuming this
                # is because Pe = Pa for a fully expanded CD Nozzle
                self.Te = self.Toi*(1-self.ni*(1-(self.Pe/self.Poi)**((self.gam-1)/self.gam)))
        else:
            # Nozzle is not choked
            self.Pe = self.Pa
        
        self.Te = self.Toi*(1-self.ni*(1-(self.Pe/self.Poi)**((self.gam-1)/self.gam)))
        self.Me = np.sqrt((2/(self.gam-1))*(self.Toi/self.Te - 1))
        self.Ve = self.Me*np.sqrt(self.gam*self.R*self.Te)
        
        # Stag props at exit
        self.Toe = self.Toi
        self.Poe = self.Pe * (1 + (self.gam -1)*(self.Me**2)/2)**(self.gam/(self.gam -1))
       

        
