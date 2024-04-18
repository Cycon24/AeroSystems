# -*- coding: utf-8 -*-
"""
Created on Thu Aug 24 09:45:15 2023

@author: cycon
"""
import sys
import numpy as np
import EngineErrors as EngineErrors
import EnginePerformanceFunctions as EPF
from WorkDoneFactor import interp_wdf
import pandas as pd

from freeVortexCompressor import freeVortexCompressorMeanLine

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
        self.R = 287 # J/kg*K
        self.gam_a = 1.4
        self.gam_g = 4/3
        self.cp_a = 1.005 # kJ/kg*K
        self.cp_g = 1.148 # kJ/kg*K
        # General properties that could be used by all stages 
        # so all components know the atm conditions
        self.Ta  = kwargs.get('Ta') 
        self.Pa  = kwargs.get('Pa')
        self.Vinf = kwargs.get('Vinf')
        self.Minf = kwargs.get('Minf')
        
        self.Toi = kwargs.get('Toi')
        self.Poi = kwargs.get('Poi')
        self.Ti  = kwargs.get('Ti')
        self.Pi  = kwargs.get('Pi')
        self.Mi  = kwargs.get('Mi')
        self.Vi  = kwargs.get('Vi')
        
        self.m_dot = kwargs.get('m_dot') # Stays constant through component
        self.mdot_ratio = 1 # used to track mass flow ratio through sections
        self.ni   = kwargs.get('ni', 1) # Isentropic efficiency
        self.BPR = kwargs.get('BPR',1)
        
        self.Toe = kwargs.get('Toe')
        self.Poe = kwargs.get('Poe')
        self.Te  = kwargs.get('Te')
        self.Pe  = kwargs.get('Pe')
        self.Me  = kwargs.get('Me')
        self.Ve  = kwargs.get('Ve')
        
        self.StageName = ""
        self.Power = None
        self.SpecPower = None
        
    def forward(self, next_Stage):
        next_Stage.Toi = self.Toe
        next_Stage.Poi = self.Poe
        next_Stage.Ti  = self.Te
        next_Stage.Pi  = self.Pe
        next_Stage.Mi  = self.Me
        next_Stage.Vi  = self.Ve
        next_Stage.m_dot = self.m_dot
        next_Stage.mdot_ratio = self.mdot_ratio
    
    def printOutputs(self):
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
       
            
class compressor_stage():
    # Rotor and Stator, use 1 to denote inflow to rotor, 2 for between rotor
    # and stator, 3 for exit of stator (and therefore stage)
    def __init__(self, stage_number, Lamda, dTo_s, **inputs):
        self.r_m = inputs.get('r_m')
        # self.w = inputs.get('omega')
        self.stage_num = stage_number
        self.dTo = dTo_s
        self.Ca = inputs.get('Ca')
        self.nc = inputs.get('nc_s')
        self.U_m = inputs.get('U_m')
        self.mdot = inputs.get('mdot_c')
        self.N = inputs.get('N') # Rev/s
        
        # Limit Props
        self.de_Haller_min = 0.65
        self.tip_Mach_max = 1.2
        self.deHaller_Failed_Flag = False
        self.MachExceeded_Flag = False
        
        # Gas Props
        self.cp = 1005 # J/kg*K
        self.R = 287   # J/kg*K
        self.gam = 1.4
        
        # Inlet conditions
        self.Toi = inputs.get('Toi')
        self.Poi = inputs.get('Poi')
        self.Toe = None
        self.Poe = None
        
        self.Lamda = Lamda
        self.WDF = inputs.get('WorkDoneFactor')
        
        # self.Ca1 = inputs.get('Ca1')
        # self.Ca2 = inputs.get('Ca2')
        # self.Ca3 = inputs.get('Ca3')
        
        self.C1 = inputs.get('C1')
        self.C2 = inputs.get('C1')
        self.C3 = inputs.get('C3')
        self.Cw1 = inputs.get('Cw1')
        self.Cw2 = inputs.get('Cw2')
        self.Cw3 = inputs.get('Cw3')
        
        self.V1 = inputs.get('V1')
        self.V2 = inputs.get('V2')
        # self.V3 = inputs.get('V3')
        self.Vw1 = inputs.get('Vw1')
        self.Vw2 = inputs.get('Vw2')
        # self.Vw3 = inputs.get('Vw3')
        
        self.alpha_1 =inputs.get('alpha_1')
        self.alpha_2 =inputs.get('alpha_2')
        self.alpha_3 =inputs.get('alpha_3')
        self.beta_1 =inputs.get('beta_1')
        self.beta_2 =inputs.get('beta_2')
        # self.beta_3 =inputs.get('beta_3')
        
        
    # Velocity Triangle
    #                    /|\  
    #                   /β|α\
    #               V1 /  |  \ C1 
    #                 /   |   \   
    #                /    |Ca  \
    #               /     |     \
    #              <------|------>
    #                Vw1     Cw1
    #              -------------->
    #                     U 
    #           /P
    #          //       Rotor 
    #          ||       ----> U
    #          ||
    #           \
    def calculate(self):
        # Currently assuming Constant Axial Velocity
        wdf = interp_wdf(self.stage_num) if self.WDF == None else self.WDF
        self.WDF = wdf
        To1 = self.Toi
        Po1 = self.Poi
        if self.Lamda == None:
            # First stage, need to find lamda, assuming Cw1 = 0 (no inlet swirl)
            self.Cw1 = 0 if self.Cw1 == None else self.Cw1
            # Solve for whirl difference and whirl vel
            dCw = self.cp*self.dTo / (wdf*self.U_m) # Cw2 - Cw1 = dCw

            self.Cw2 = dCw + self.Cw1
            # Solve for Gas angles
            self.beta_1 = np.arctan((self.U_m - self.Cw1) / self.Ca)
            self.beta_2 = np.arctan((self.U_m - self.Cw2) / self.Ca)
            self.alpha_1 = np.arctan(self.Cw1/self.Ca) # Should be 0 for first stage
            self.alpha_2 = np.arctan(self.Cw2/self.Ca)
            # Calculate Lamda
            self.Lamda = 1 - (self.Cw1 + self.Cw2)/(2*self.U_m)
        else: 
            # We already know lamda
            # Equations to help solve for betas
            eq1= (self.cp*self.dTo) / (wdf*self.U_m*self.Ca) # = tan(B1) - tan(B2)
            eq2= 2*self.Lamda*self.U_m/self.Ca # = tan(B1) + tan(B2)
            
            # Solve for gas angles
            self.beta_1 = np.arctan((eq1+eq2)/2)
            self.beta_2 = np.arctan(eq2 - np.tan(self.beta_1))
            self.alpha_1 = np.arctan(self.U_m/self.Ca  - np.tan(self.beta_1))
            self.alpha_2 = np.arctan(self.U_m/self.Ca  - np.tan(self.beta_2))
            # Solve for whirl velocities
            self.Cw1 = self.Ca*np.tan(self.alpha_1)
            self.Cw2 = self.Ca*np.tan(self.alpha_2)
            
        # The rest should be the same wether first, last, or intermediate stage    
        # Solve for absolute gas velocities
        self.C1 = np.sqrt(self.Cw1**2 + self.Ca**2)
        self.C2 = np.sqrt(self.Cw2**2 + self.Ca**2)
        # Solve for other gas-triangle properties
        # self.Vw1 = self.U_m - self.Cw1
        # self.Vw2 = self.U_m - self.Cw2
        # self.V1 = np.sqrt(self.U_m**2 - self.C1**2)
        # self.V2 = np.sqrt(self.U_m**2 - self.C2**2)
        # Solve for static params
        self.T1 = To1 - (self.C1**2)/(2*self.cp)
        self.P1 = Po1*(self.T1/To1)**(self.gam/(self.gam-1))
        self.rho1 = self.P1/(self.R*self.T1)
        # Solve for blade dimensions
        self.h = self.mdot/(2*np.pi * self.rho1 * self.Ca * self.r_m)
        self.r_t = self.r_m + self.h/2
        self.r_r = self.r_m - self.h/2


        obj = freeVortexCompressorMeanLine(self.beta_1, self.beta_2, self.alpha_1, self.alpha_2, self.Ca, self.Lamda, self.r_m, To1, self.N, self.dTo, self.WDF)

        try:
            print("Stage", self.stage_num)
            self.mean = obj.calculate(self.r_m, "mean")
            mean_data = self.mean.toDataFrame()
            # self.mean.print()
            self.tip  = obj.calculate(self.r_t, "tip")
            tip_data = self.tip.toDataFrame()
            # self.tip.print()
            self.data = pd.concat([mean_data, tip_data], keys=2*[f"Stage {self.stage_num}"])
            self.root = obj.calculate(self.r_r, "root")
            root_data = self.root.toDataFrame()
            # root_data.
            self.data = pd.concat([root_data, mean_data, tip_data], keys=3*[f"Stage {self.stage_num}"])
            #print(stage_data)
            #stage_data = pd.DataFrame(stage_data, index=[f"Stage {self.stage_num}"])
            # self.root.print()
        except ValueError as e:
            print(f"While running {self.stage_num}")
            print("  "+str(e))
            exit()

        # # Checking tip mach
        # self.Cw1_t = self.Cw1 * self.r_m / self.r_t
        # self.C1_t = np.sqrt(self.Cw1_t**2 + self.Ca**2)
        # T1_t = To1 - (self.C1_t**2)/(2*self.cp)
        # self.V1_t = np.sqrt((2*np.pi*self.N*self.r_t - self.Cw1_t)**2 + self.Ca**2)
        # self.M1_t = self.V1_t / np.sqrt(self.gam*self.R*T1_t)
        
        # if self.M1_t > self.tip_Mach_max:
        #     print('Max tip Mach exceeded in stage {}, M_t = {:.4f}'.format(self.stage_num, self.M1_t))
        #     self.MachExceeded_Flag = True
        
        # R1_tm = self.r_t/self.r_m 
        # R1_rm = self.r_r/self.r_m
        
        # # Get Degree of Reactions at root and tip
        # self.Lamda_t = 1 - (1 - self.Lamda)/R1_tm**2
        # self.Lamda_r = 1 - (1 - self.Lamda)/R1_rm**2

        # if self.Lamda_r < 0:
        #     print('root degree of reaction check failed in stage {}, Lambda_r = {:.4f}'.format(self.stage_num, self.Lamda_r))
        # if self.Lamda_t < 0:
        #     print('tip degree of reaction check failed in stage {}, Lambda_t = {:.4f}'.format(self.stage_num, self.Lamda_t))
        
        # # Check de Haller Criteria
        # deHaller = np.cos(self.beta_1)/np.cos(self.beta_2)
        # if deHaller < self.de_Haller_min:
        #     print('de Haller check failed in stage {}, V2/V1 = {:.4f}'.format(self.stage_num, deHaller))
        #     self.deHaller_Failed_Flag = True
        
        # # Ca = const for free vortex so Ca_t = Ca
        # self.alpha_1t = np.arctan(self.Cw1_t/self.Ca)
        # self.U_t = 2*np.pi*self.r_t
        # self.Vw1_t = -1*(self.U_t - self.Cw1_t)
        # self.beta_1t = np.arctan(self.Vw1_t/self.Ca)

        # self.Cw2_t = self.Cw2 * self.r_m / self.r_t
        # self.C2_t = np.sqrt(self.Cw2_t**2 + self.Ca**2)
        # T2_t = To1 - (self.C2_t**2)/(2*self.cp)
        # self.V2_t = np.sqrt((2*np.pi*self.N*self.r_t - self.Cw2_t)**2 + self.Ca**2)
        # self.M2_t = self.V2_t / np.sqrt(self.gam*self.R*T2_t)

        # self.alpha_2t = np.arctan(self.Cw2_t/self.Ca)
        # self.Vw2_t = -1*(self.U_t - self.Cw2_t)
        # self.beta_2t = np.arctan(self.Vw2_t/self.Ca)

        # # Check de Haller Criteria
        # deHaller = np.cos(self.beta_1t)/np.cos(self.beta_2t)
        # if deHaller < self.de_Haller_min:
        #     print('de Haller check at the tip failed in stage {}, V2/V1 = {:.4f}'.format(self.stage_num, deHaller))
        #     self.deHaller_Failed_Flag = True

        # if self.M2_t > self.tip_Mach_max:
        #     print('Max tip Mach exceeded in exit of stage {}, M_2t = {:.4f}'.format(self.stage_num, self.M2_t))
        #     self.MachExceeded_Flag = True

        # # Checking tip mach
        # self.Cw1_r = self.Cw1 * self.r_m / self.r_r
        # self.C1_r = np.sqrt(self.Cw1_r**2 + self.Ca**2)
        # T1_r = To1 - (self.C1_r**2)/(2*self.cp)
        # self.V1_r = np.sqrt((2*np.pi*self.N*self.r_r - self.Cw1_r)**2 + self.Ca**2)
        # self.M1_r = self.V1_r / np.sqrt(self.gam*self.R*T1_r)

        # # Ca = const for free vortex so Ca_t = Ca
        # self.alpha_1r = np.arctan(self.Cw1_r/self.Ca)
        # self.U_r = 2*np.pi*self.r_r
        # self.Vw1_r = -1*(self.U_r - self.Cw1_r)
        # self.beta_1r = np.arctan(self.Vw1_r/self.Ca)

        # self.Cw2_r = self.Cw2 * self.r_m / self.r_r
        # self.C2_r = np.sqrt(self.Cw2_r**2 + self.Ca**2)
        # T2_r = To1 - (self.C2_r**2)/(2*self.cp)
        # self.V2_r = np.sqrt((2*np.pi*self.N*self.r_r - self.Cw2_r)**2 + self.Ca**2)
        # self.M2_r = self.V2_r / np.sqrt(self.gam*self.R*T2_r)

        # self.alpha_2r = np.arctan(self.Cw2_r/self.Ca)
        # self.Vw2_r = -1*(self.U_t - self.Cw2_r)
        # self.beta_2r = np.arctan(self.Vw2_r/self.Ca)

        # # Calcualte exit stag properties
        self.Poe = Po1*(1 + self.nc*self.dTo/To1)**(self.gam/(self.gam - 1))
        self.Toe = To1 + self.dTo
        
        return self.MachExceeded_Flag, self.deHaller_Failed_Flag
    
    
    def forward(self, next_stage):
        next_stage.Toi = self.Toe
        next_stage.Poi = self.Poe

    # def calcVelocityTriangleParams(self, radius)
        
    def printVelocityTrianges(self):
        num1s = 1 + (self.stage_num - 1)*3
        num2s = num1s + 1 
        num3s = num2s + 1
        print('----------------------')
        print('Stage ',self.stage_num, 'λ = {:4.3f}'.format(self.WDF))
        print('  Mean Line Properties')
        print('  Tip Properties')
        print('  Root Properties')
        print('----------------------')
        # print('  Gas Angles')
        # print('\t α_{:.0f} {:6.3f}°'.format(num1s, np.degrees(self.alpha_1)))
        # print('\t Β_{:.0f} {:6.3f}°'.format(num1s, np.degrees(self.beta_1)))
        # print('\t α_{:.0f} {:6.3f}°'.format(num2s, np.degrees(self.alpha_2)))
        # print('\t Β_{:.0f} {:6.3f}°'.format(num2s, np.degrees(self.beta_2)))
        # print('\t Λ_{:.0f} {:6.3%}'.format(self.stage_num, self.Lamda))
        # print('  Absolute Vels:')
        # print('\t C_{:.0f}  {:6.3f} m/s'.format(num1s, self.C1))
        # print('\t Cw_{:.0f} {:6.3f} m/s'.format(num1s, self.Cw1))
        # print('\t C_{:.0f}  {:6.3f} m/s'.format(num2s, self.C2))
        # print('\t Cw_{:.0f} {:6.3f} m/s'.format(num2s, self.Cw2))
        # print('  Relative Vels:')
        # print('\t V_{:.0f}  {:6.3f} m/s'.format(num1s, self.V1))
        # print('\t Vw_{:.0f} {:6.3f} m/s'.format(num1s, self.Vw1))
        # print('\t V_{:.0f}  {:6.3f} m/s'.format(num2s, self.V2))
        # print('\t Vw_{:.0f} {:6.3f} m/s'.format(num2s, self.Vw2))
        # print(' Tip Properties')
        # print('  Gas Angles')
        # print('\t α_{:.0f} {:6.3f}°'.format(num1s, np.degrees(self.alpha_1t)))
        # print('\t Β_{:.0f} {:6.3f}°'.format(num1s, np.degrees(self.beta_1t)))
        # print('\t α_{:.0f} {:6.3f}°'.format(num2s, np.degrees(self.alpha_2t)))
        # print('\t Β_{:.0f} {:6.3f}°'.format(num2s, np.degrees(self.beta_2t)))
        # print('\t Λ_{:.0f} {:6.3%}'.format(self.stage_num, self.Lamda_t))
        # print('  Absolute Vels:')
        # print('\t C_{:.0f}  {:6.3f} m/s'.format(num1s, self.C1_t))
        # print('\t Cw_{:.0f} {:6.3f} m/s'.format(num1s, self.Cw1_t))
        # print('\t C_{:.0f}  {:6.3f} m/s'.format(num2s, self.C2_t))
        # print('\t Cw_{:.0f} {:6.3f} m/s'.format(num2s, self.Cw2_t))
        # print('  Relative Vels:')
        # print('\t V_{:.0f}  {:6.3f} m/s'.format(num1s, self.V1_t))
        # print('\t Vw_{:.0f} {:6.3f} m/s'.format(num1s, self.Vw1_t))
        # print('\t V_{:.0f}  {:6.3f} m/s'.format(num2s, self.V2_t))
        # print('\t Vw_{:.0f} {:6.3f} m/s'.format(num2s, self.Vw2_t))
        # print(' Root Properties')
        # print('  Gas Angles')
        # print('\t α_{:.0f} {:6.3f}°'.format(num1s, np.degrees(self.alpha_1r)))
        # print('\t Β_{:.0f} {:6.3f}°'.format(num1s, np.degrees(self.beta_1r)))
        # print('\t α_{:.0f} {:6.3f}°'.format(num2s, np.degrees(self.alpha_2r)))
        # print('\t Β_{:.0f} {:6.3f}°'.format(num2s, np.degrees(self.beta_2r)))
        # print('\t Λ_{:.0f} {:6.3%}'.format(self.stage_num, self.Lamda_r))
        # print('  Absolute Vels:')
        # print('\t C_{:.0f}  {:6.3f} m/s'.format(num1s, self.C1_r))
        # print('\t Cw_{:.0f} {:6.3f} m/s'.format(num1s, self.Cw1_r))
        # print('\t C_{:.0f}  {:6.3f} m/s'.format(num2s, self.C2_r))
        # print('\t Cw_{:.0f} {:6.3f} m/s'.format(num2s, self.Cw2_r))
        # print('  Relative Vels:')
        # print('\t V_{:.0f}  {:6.3f} m/s'.format(num1s, self.V1_r))
        # print('\t Vw_{:.0f} {:6.3f} m/s'.format(num1s, self.Vw1_r))
        # print('\t V_{:.0f}  {:6.3f} m/s'.format(num2s, self.V2_r))
        # print('\t Vw_{:.0f} {:6.3f} m/s'.format(num2s, self.Vw2_r))

        
class turbine_free_vortex():
    def __init__(self, CycleTurbineObject, numStages, psi_m1, phi_m1, Lamda, npt, N, PHI_MIN=0.78, PSI_MAX=3.3):
        # Retrieve Tubine inlet and outlet values
        self.Toi = CycleTurbineObject.Toi
        self.Toe = CycleTurbineObject.Toe
        self.Poi = CycleTurbineObject.Poi
        self.Poe = CycleTurbineObject.Poe
        self.dTo_cycle = self.Toi - self.Toe
        self.mdot = CycleTurbineObject.m_dot
        self.PHI_MIN = PHI_MIN
        self.PSI_MAX = PSI_MAX
        
        self.psi_m = psi_m1
        self.phi_m = phi_m1
        self.Lamda = Lamda
        self.nt_p = npt
        self.N = N
        
        self.gam = CycleTurbineObject.gam_g
        self.cp = CycleTurbineObject.cp_g * 1000 # convert to J/kg*K
        self.R = CycleTurbineObject.R 
        
        # Assumes even work distribustion between stages
        self.dTo_stage = self.dTo_cycle/numStages # K
        
        self.stages = []
        for i in range(0, numStages):
            if i == 0:
                self.stages.append(turbine_stage(self.Poi, self.Toi, self, i+1))
            else:
                self.stages.append(turbine_stage(self.stages[i-1].Po3, self.stages[i-1].To3, self, i+1))
    def calculate(self):
        for stage in self.stages:
            stage.calculate()
                
class turbine_stage():
    def __init__(self, Poi, Toi, TurbineUnit_FV, stageNum):
        self.stageNum = stageNum
        self.gam = TurbineUnit_FV.gam
        self.cp = TurbineUnit_FV.cp # J/kg*K
        self.R = TurbineUnit_FV.R 
        nt_p = TurbineUnit_FV.nt_p
        self.Lamda = TurbineUnit_FV.Lamda
        self.N = TurbineUnit_FV.N
        self.mdot = TurbineUnit_FV.mdot
        self.PHI_MIN= TurbineUnit_FV.PHI_MIN
        self.PSI_MAX= TurbineUnit_FV.PSI_MAX
        
        self.dTo_s = TurbineUnit_FV.dTo_stage # Still assuming even work distribution between stages
        self.phi_m = TurbineUnit_FV.phi_m 
        self.psi_m = TurbineUnit_FV.psi_m
        self.Po1 = Poi
        self.To1 = Toi
        self.Po3 = Poi*(1 - self.dTo_s/Toi)**(self.gam/(nt_p*(self.gam - 1))) 
        self.To3 = Toi - self.dTo_s
        
    def calculate(self):
        name = f"Stage {self.stageNum}"
        # Velocities at Mean
        U_m = np.sqrt(2*self.cp*self.dTo_s/self.psi_m)
        Ca = self.phi_m * U_m
        
        
        # MEAN LINE
        # Gas angles
        phi = self.phi_m
        psi = self.psi_m
        beta_2m  = np.arctan(0.5/phi * (psi/2 - 2*self.Lamda)) # rad
        beta_3m  = np.arctan(0.5/phi * (psi/2 + 2*self.Lamda)) # rad
        alpha_2m = np.arctan(np.tan(beta_2m) + 1/phi) # rad
        alpha_3m = np.arctan(np.tan(beta_3m) - 1/phi) # rad
        
        # Mean Radius
        r_m = U_m / (2*np.pi*self.N) # m
        
        # Obtaining Bade Height and radii
        C_3m          = Ca/np.cos(alpha_3m) # m/s
        T_3m          = self.To3 - C_3m**2 / (2*self.cp) # K
        P_3m          = self.Po3 * (T_3m/self.To3)**(self.gam/(self.gam - 1))
        rho_3m        = P_3m / (self.R*T_3m) # kg/m^3 - Assumed pressure in Pa
        h           = self.mdot / (2*np.pi*rho_3m*Ca*r_m) # m
        r_t         = r_m + h/2 # m
        r_r         = r_m - h/2 # m
        # note: r2's = r3's
        
        U_r = 2*np.pi*self.N*r_r
        U_t = 2*np.pi*self.N*r_t
        
        psi_r = 2*self.cp*self.dTo_s / U_r**2 
        psi_t = 2*self.cp*self.dTo_s / U_t**2 
        phi_r = Ca/U_r 
        phi_t = Ca/U_t
        
        # Tip calculations
        C_w2m        = self.cp * self.dTo_s / U_m # m/s
        C_w2t        = C_w2m * r_m/r_t # m/s
        C_2t         = np.sqrt(C_w2t**2 + Ca**2) # m/s
        
        V_2t         = np.sqrt((C_w2t - U_t)**2 + Ca**2) # m/s
        beta_2t      = np.arccos(Ca / V_2t) # rad
        alpha_2t     = np.arctan(C_w2t / Ca) # rad
        T_2t         = self.To1 - C_2t**2 / (2*self.cp) # K
        M_2t         = V_2t / np.sqrt(self.gam * self.R * T_2t)
        # ADD THIS SOMEWHERE ELSE LATER
        self.Mach_max = 1.2
        
        
        # if M_2t > self.Mach_max:
        #     raise ValueError('Max Mach exceeded at {} : M_1 = {:.4f}'.format(name, M_2t))
        # if psi_r > self.PSI_MAX:
        #     raise ValueError('Min Psi exceeded at {} : psi_r = {:.4f}'.format(name, psi_r))
        # if phi_r < self.PHI_MIN:
        #     raise ValueError('Max Phi exceeded at {} : phi_r = {:.4f}'.format(name, phi_r))
        # if psi_t > self.PSI_MAX:
        #     raise ValueError('Min Psi exceeded at {} : psi_t = {:.4f}'.format(name, psi_t))
        # if phi_t < self.PHI_MIN:
        #     raise ValueError('Max Phi exceeded at {} : phi_t = {:.4f}'.format(name, phi_t))
        
        
        
        
class combustor_component():
    def __init__(self):
        # Possible inputs/Outputs:
        # mdot
        # mdot_ratio ?
        # mdot_fuel
        # A_inlet
        # A_mean
        # A_exit
        # phi - Overall equivalence ratio
        # phi_local - Local equivalence ratio near flame
        # mdot_local - Local mass flow rate near flame
        # dPo_full - Full load pressure loss
        # dPo - Pressure loss at current conditions
        # V_inlet - Average Inlet velocity
        # V_exit - Average Exit velocity
        # To_i - Initial stagnation temperature
        # To_e - Exit stagnation temperature
        # Po_i - Inlet stagnation pressure
        # Po_e - Exit stagnation pressure
        # f_stoich - Stoicheometric fuel/air ratio
        # f_ideal  - Ideal fuel/air ratio
        # f_actual - Actual fuel/air ratio
        # f_local  - Local fuel/air ratio
        # nb - combustor efficiency
        
        # Intermediates
        # K1
        # K2 
        # X - 1/phi
        
        # General Values
        self.R = 287 # J/kg*K
        self.cpa = 1.005 # kJ/kg
        self.cpg = 1.148 # kJ/kg
        self.gam_a = 1.4 
        self.gam_g = float(4/3)
        self.MW_air = 137.3653 # kg/kmol
        self.MW_C = 12.011 # kg/kmol
        self.MW_H = 1.0078 # kg/kmol
        
        # Desired Calculations:
        # Pressure Loss
        # Pressure Loss at part-load
        # Inlet Area
        # Exit Area
        # Mean Area
        # Inlet Velocity
        # Exit Velocity
        # Fuel to air ratio
        # Overall equivalence ratio
        # Local mass flow rate
        # Local equivalence ratio
        
    def rho(self, To, Po, mdot, A, Gamma=None):
        Gamma = self.gam_a if Gamma==None else Gamma
        P = self.SOMEITERATIVESOLVER(Po, To, mdot, A)
        a = To 
        b = -P/self.R 
        c = - (Gamma-1)*mdot / (2*Gamma*self.R*A**2)
        rho = (-b + np.sqrt(b**2 - 4*a*c)) / (2*a) # Dont need the - for the quad
        return rho
    
    def rho(self, To, Po, Mach, Gamma=None):
        Gamma = self.gam_a if Gamma==None else Gamma 
        rho_o = Po / (self.R*To)
        rho = rho_o / ((1 + (Gamma-1)*(Mach**2)/2)**(1/(Gamma-1)) )
        return rho
     
    def SOMEITERATIVESOLVER(self, Po, To, mdot, A):
        print("Working on derivation (its uggy)")
        
class Turbofan_SingleSpool():
    def __init__(self, **kwargs):
        '''
        A signgle spool turbofan which has one turbine to power the
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
            'Minf': 0.85, # Or Vinf, if none its assumed stationary
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
        self.inlet = Intake(**gen_kwargs,ni=ni,m_dot=mdot)
        self.fan = Compressor(**gen_kwargs, rc=rfan, np=npf, ni=nf)
        self.BP_nozzle = Nozzle('cold',**gen_kwargs, ni=nj)
        self.HP_comp = Compressor(**gen_kwargs, rc=rc, np=npc, ni=nc)
        self.combustor = Combustor(**gen_kwargs, Toe=To_ti, dPb_dec=dP_b, ni=nb)
        self.HP_turb = Turbine([self.fan, self.HP_comp], **gen_kwargs, nm=nm, ni=nt, np=npt)
        self.nozzle = Nozzle(**gen_kwargs, ni=nj) # Nozzle/Exhaust?
        
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
        self.inlet = Intake(**gen_kwargs,ni=ni,m_dot=mdot)
        self.fan = Compressor(**gen_kwargs, rc=rfan, np=npf, ni=nf)
        self.BP_nozzle = Nozzle('cold',**gen_kwargs, ni=nj)
        self.HP_comp = Compressor(**gen_kwargs, rc=rc, np=npc, ni=nc)
        self.combustor = Combustor(**gen_kwargs, Toe=To_ti, dPb_dec=dP_b, ni=nb)
        self.HP_turb = Turbine(self.HP_comp, **gen_kwargs, nm=nm, ni=nt, np=npt)
        self.LP_turb = Turbine(self.fan, **gen_kwargs, nm=nm, ni=nt_lp, np=npt_lp)
        self.nozzle = Nozzle(**gen_kwargs, ni=nj) # Nozzle/Exhaust?
        
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