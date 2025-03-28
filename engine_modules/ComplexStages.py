# -*- coding: utf-8 -*-
"""
Created on Sun Jun 30 07:38:38 2024

@author: cycon

20240630:
    Created from previously named file EngineModule.py, which is renamed to SimpleStages.py
    
"""
import numpy as np
# from engine_modules.WorkDoneFactor import interp_wdf
import pandas as pd

# from freeVortexCompressor import freeVortexCompressorMeanLine
            
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