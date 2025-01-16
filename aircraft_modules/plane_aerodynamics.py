# -*- coding: utf-8 -*-
"""
Created on Sat Feb 25 15:25:12 2023

@author: cycon
Update-Fixes Log:
    08/24/23- I dont like how I have the planform setup within the plane class
                Id rather only have one set of variables, may just need to go through
                and replace all times the plane's vars are used with the planform.
                Also means I need to replace in the creation of a plane
            - Would like to add a function called checks, that runs any time a function
            is requested by the user. It will check which variables are still at None
            and will calculate anything for which it has the suitable values to use\
            - Make a function that calculates all flight characteristics at a particular 
            altitude and load factor an return a dictionary (aka kwargs). This will allow 
            the use of a single function to hanle the calculations, then pull the needed
            values for graphing
    01/08/25 - Finsihed HW8 calculations (Through range factor)
            - added "CHECK" comments to refer to later
            - Checked and fixed end factor and range factor calcs and verified against HW plots
            - Finsihed HW8 calcs in Calculate_FlightCharacteristics but still need to add functions
            to utilize the values
            
    01/11/25 
            - Addded HW10 plotting function, will need revised to better place subplots (include multiple cols) for 
            multiple altitude inputs (currently only adds rows)
            - Resolved some calculation issues
            
To Do:
    - Make functions for HW8 use cases (fuel burned/left)
    - Make functions for HW9 use cases (max alt or speed attainable at energy (h and v))
    - Make KEAS to V function?  KTAS = KEAS/sqrt(DR),  V= KTAS*(6076.12/3600)
    - Make function 
    - Revise Plot_Thrust_vs_Mach to better place subplots 
            
            
"""
import numpy as np
import math
import matplotlib.pyplot as plt

if __name__=='__main__':
    import sys 
    from pathlib import Path
    current_file_path = Path(__file__).resolve()
    parent_directory = current_file_path.parent.parent
    
    sys.path.append(str(parent_directory))
    
from atmosphere_modules.AtmosphereModule import Atmosphere
# Plane Class
# '''
# Homeworks of interest
# HW  5 - Low - For the most part done. Could include more weight options if wanted to?
# HW  7 - High - details (idk what I meant by this) - Done
# HW  8 - High - Done (but needsd range and endurance time calculators for X lb fuel burned)
# HW  9 - Med  - Done (but needs specific V and h combos for specific outputs, aka with certain flying conditions what are max altitudew and max speed)
#                - also need to add function to calculate speed difference between 2 altitudes (or vice versa) with a constant energy
# HW 10 - High - Done 
# HW 11 - High
# HW 12 - High
# HW 13 - High 
# '''

class planform():
    
    def __init__(self):
        # Parameters that matter to a planform
        # ___ Wing Design Parameters ___ 
        # Inputs
        self.S_ref = None # ft^2 - Wing Reference area
        self.b = None     # ft   - Wing Span
        self.LE_Sweep = None # deg - Sweep of leading edge - Λ
        self.c_t = None   # ft - Tip Chord
        self.c_r = None   # ft - Root Chord
        
        # Calculated
        self.AR_ref = None  # Aspect Ration using S_ref
        self.TaperRatio = None   # Taper Ratio - λ
        self.MAC = None     # ft - Mean Aerodynamic Chord
        self.y_mac = None   # ft - istance from centerline of MAC

    
    def Calculate(self):
        # Calculations
        S_trap = (self.c_r + self.c_t)*self.b / 2
        self.AR_ref = self.b**(2/self.S_ref)
        AR_trap = self.b**(2/S_trap)
        self.TaperRatio = self.c_t / self.c_r
        self.MAC = (2*self.c_r/3)*(1 + self.TaperRatio \
                                   + self.TaperRatio**2)/(1+self.TaperRatio)
        self.y_mac =self.b/6 * (1 + 2*self.TaperRatio) / (1 + self.TaperRatio)
        
    def Plot(self):
        # Sweep Angles: x/c
        xc = np.arange(0,1,0.1)
        sweep_xc = np.zeros(np.size(xc))
        for i in range(0,np.size(xc)):
            sweep_xc[i] = 180*math.atan(np.tan(self.LE_Sweep*np.pi/180) \
                            -xc[i]*(2*self.c_r*(1-self.TaperRatio)/self.b))/np.pi
        
        # Spanwise View
        x = np.zeros(5)
        y_pos = x.copy()
        y_neg = x.copy()
        
        x[1] = self.c_r
        x[2] = self.b/2*np.tan(np.radians(sweep_xc[0])) + self.c_t
        x[3] = self.b/2*np.tan(np.radians(sweep_xc[0]))
        
        print(x)
        
        y_pos[2:4] = self.b/2
        y_neg = -y_pos
        
        fig, ax = plt.subplots(1, 1)
        ax.plot(x,y_pos, x,y_neg)
        ax.set_xlim([0,self.c_r+0.3*self.c_r])
        ax.set_ylim([-self.b/2-0.1*self.b,self.b/2+0.1*self.b])
        ax.set_xlabel('Along Centerline')
        ax.set_ylabel('Along Span')
        ax.set_aspect('equal','box')
        plt.tight_layout()
        plt.grid()
    
    def ChordFunct(self):
        # y is location along span from -b/2 to b/2
        # will return a function that inputs y and returns chord length
        
        # Sweep Angles: x/c
        xc = np.arange(0,1,0.1)
        sweep_xc = np.zeros(np.size(xc))
        for i in range(0,np.size(xc)):
            sweep_xc[i] = 180*math.atan(np.tan(self.LE_Sweep*np.pi/180) \
                            -xc[i]*(2*self.c_r*(1-self.TaperRatio)/self.b))/np.pi
        
        # Spanwise View
        x_LE = np.zeros(2)
        x_TE = x_LE.copy()
        y = x_LE.copy()
        
        y = np.array([0, self.b/2])
        x_LE[0] = 0
        x_LE[1] = self.b/2*np.tan(np.radians(sweep_xc[0]))
        x_TE[0] = self.c_r
        x_TE[1] = self.b/2*np.tan(np.radians(sweep_xc[0])) + self.c_t
        
        # Make line equations (For Check)
        # def LE_x(yi):
        #     return (x_LE[1]/y[1])*yi
        
        # def TE_x(yi):
        #     return x_TE[0] + yi*(x_TE[0] - x_TE[1])/(y[0] - y[1])
        
        def chord_f(yi):
            yi = abs(yi)
            LE_x = (x_LE[1]/y[1])*yi
            TE_x = x_TE[0] + yi*(x_TE[0] - x_TE[1])/(y[0] - y[1])
            return TE_x - LE_x
        
        return chord_f
    
class Plane():
    # Possible changes/Improvements:
        # i would like to change the planform to have
        # a wing and a tail rather than having two
        # sets of variables for root cords, tip cords,
        # taper ratios, etc.
        # S_ref should be the only one that is included 
        # in plane variables since it is used in many
        # plnae characteristics
        
    
    def __init__(self):
        # ___ Aircraft Parameters ______
        # Limits
        self.max_n = None  # g's
        self.min_n = None # g's
        self.max_V_fts = None # ft/s - Max Vel
        self.max_V_KEAS = None # Max KEAS
        self.max_M = None     # Max Mach Number
        self.max_h = None  # ft - Maximum altitude
        
        # Weights
        self.W = None  # lbfs - Performance Calculation default Weight
        self.FuelCap = None # lbfs - Fuel Capacity
        self.TOGW = None # lbfs - Gross Take Off Weight
        self.OpWeight = None # Operating Weight Empty
        
        # Performance
        self.Ta_SL = None # lbf - Max Thrust at Sea Level
        self.sfc_SL = None   # lbm/(lbf-hr) - Specific Fuel Consumption (Fuel Flow / Thrust)
        self.CL_max = None # Max Lift Coefficient steady flight
        self.CL_max_TO = None # Max Lift Coefficient with Take Off flaps
        self.CL_max_L = None  # Max Lift Coefficient with Landing flaps
        self.CD0 = None # Zero-Lift Drag Coefficient
        self.K = None   # Drag Polar coef
        self.LD_max = None # Max L/D 
        self.TrW_min = None # Min thrust required/weight (is at V_LD_max)
        self.D_min = None # Min drag aka Thrust required (is at V_LD_max)
        self.V_LD_max = None # Velocity at max L/D
        
        # Cd = Cd0 + K*CL^2
        # CL = n*W / (0.5*rho*S*V^2) = n*W / (q*S)
        
        # Propulsion Type
        self.PropulsionTypes = {
            "Piston": None, 
            "Turboprop": None,
            "HighBPTurboFan":None,
            "LowBPTurboFan": None,
            "MedBPTurboFan": None,
            "Turbojet": None,
            "Afterburner": None}
        
        self.numEngines = None
        self.T_SL_perEngine = None
        
        
        '''
        Piston
            Ta = SHP_SL*(η_p / V)(ρ / ρ_SL) # Thrust Available
            FFR = SHP*sfc                   # Fuel Flow Rate
        Turboprop
            Ta = ESHP_SL(η_p / V)(ρ / ρ_SL)
            FFR = ESHP*sfc
        
        sfc  in (lb/HP*hr)
        sfct in (lb/lb_t*hr)
        
        HighBPTurboFan
            Ta = Ta_SL * (0.1/M) * (ρ / ρ_SL)
            FFR = Ta *sfct_SL * (a / a_SL)
        LowBPTurboFan, MedBPTurboFan, & Turbojet
            Ta = Ta_SL * (ρ / ρ_SL)
            FFR = Ta * sfct_SL * (a / a_SL)
        Afterburner
            Ta = Ta_SL * (ρ / ρ_SL) * (1 + 0.7*M)
            FFR = Ta * sfct_SL * (a / a_SL)
        '''
        
        
        # ___ Wing Design Parameters ___ 
        # Plane planform
        self.planform = planform()
        
        # Inputs
        self.S_ref = None  # ft^2 - Wing Reference area
        self.b = None    # ft   - Wing Span
        self.LE_Sweep = None # deg - Sweep of leading edge - Λ
        self.c_t = None   # ft - Tip Chord
        self.c_r = None   # ft - Root Chord
        
        # Calculated
        self.AR_ref = None  # Aspect Ration using S_ref
        self.TaperRatio = None   # Taper Ratio - λ
        self.MAC = None     # ft - Mean Aerodynamic Chord
        self.y_mac = None   # ft - istance from centerline of MAC
        
        
        # ____ Other Calculation assistors
        self.V_stall_fraction = 0.75
        
        # Constants
        self.g = 32.174 # ft/s^2
        #______________________________________________________________________
        
    def Planform_Pull(self):
        '''
        Pulls wing planform values from planiform and sacves to aircraft.

        Returns
        -------
        None.

        '''
        # Pulls values from planform to plane
        # Inputs
        self.S_ref = self.planform.S_ref() # ft^2 - Wing Reference area
        self.b = self.planform.b     # ft   - Wing Span
        self.LE_Sweep = self.planform.LE_Sweep # deg - Sweep of leading edge - Λ
        self.c_t = self.planform.c_t   # ft - Tip Chord
        self.c_r = self.planform.c_r   # ft - Root Chord
        
        # Calculated
        self.AR_ref = self.planform.AR_ref  # Aspect Ration using S_ref
        self.TaperRatio = self.planform.TaperRatio   # Taper Ratio - λ
        self.MAC = self.planform.MAC     # ft - Mean Aerodynamic Chord
        self.y_mac = self.planform.y_mac   # ft - istance from centerline of MAC
    
    def Planform_Push(self):
        '''
        Pushes wing planform values from aircraft and saves them to the planform.
        
        Returns
        -------
        None.

        '''
        # Pulls values from plane and sets planform
        # Inputs
        self.planform.S_ref = self.S_ref() # ft^2 - Wing Reference area
        self.planform.b = self.b     # ft   - Wing Span
        self.planform.LE_Sweep = self.LE_Sweep # deg - Sweep of leading edge - Λ
        self.planform.c_t = self.c_t   # ft - Tip Chord
        self.planform.c_r = self.c_r   # ft - Root Chord
        
        # Calculated
        self.planform.AR_ref = self.AR_ref  # Aspect Ration using S_ref
        self.planform.TaperRatio = self.TaperRatio   # Taper Ratio - λ
        self.planform.MAC = self.MAC     # ft - Mean Aerodynamic Chord
        self.planform.y_mac = self.y_mac   # ft - istance from centerline of MAC
    
    def Check(self):
        if self.CD0 != None and self.K != None:
            if self.LD_max == None:
                self.LD_max = np.sqrt(1/(4*self.CD0*self.K))
            if self.TrW_min == None:
                self.TrW_min = np.sqrt(4*self.CD0*self.K)
        
        
        return None
    
    def Calculate_FlightCharacteristics(self, h, n, W, atm):
        """
        Calculates various flight characteristics for a given altitude, load factor, and aircraft weight.
    
        This method computes a range of aerodynamic and performance metrics, including stall velocities, 
        flight envelopes, lift/drag characteristics, and endurance factors. It is designed for single-altitude 
        analysis and is not optimized to handle a range of altitudes.
    
        Parameters
        ----------
        h : float
            Altitude at which to calculate the characteristics (in feet or meters, depending on the atmospheric model).
        n : float
            Load factor in g's, typically 1 for level flight, >1 for turns, and <1 for negative G maneuvers.
        W : float
            Weight of the aircraft (in consistent units with the aerodynamic model).
        atm : Atmosphere object
            Atmospheric model object containing methods and properties, such as:
            - `atm.linterp_h(h, property)`: Interpolates a given atmospheric property (e.g., speed of sound, air density).
            - `atm.a`, `atm.QMS`, `atm.rho`: Arrays of atmospheric properties across altitudes.
        
        Returns
        -------
        dict
            A dictionary containing various flight characteristics:
            - ThrustAvailable: float
                The trust available at the current altitude (max T_a = T_max*rho/rho_SL)
            - MachRange: numpy.ndarray
                  Array of Mach numbers for the flight envelope.
            - VelocityRange: numpy.ndarray
                  Array of velocities (in ft/s or m/s) corresponding to Mach numbers.
            - n_ClimbLimit: numpy.ndarray
                  Load factor limits for climbing flight at given velocities.
            - n_DescendLimit: numpy.ndarray
                  Load factor limits for descending flight at given velocities.
                  
            - StallVelocity: float
                  Stall velocity (ft/s or m/s) at the given altitude.
            - StallMach: float
                  Stall Mach number at the given altitude.
            - MinVelocity-T=D: float
                  Minimum velocity where thrust equals drag.
            - MaxVelocity-T=D: float
                  Maximum velocity where thrust equals drag.
            - MaxVelocityAtAltitude: float
                  Maximum velocity limited by KEAS at the given altitude.
            - MaxMach-T=D: float
                  Maximum Mach where thrust equals drag.
            - MaxMachAtAltitude: float
                  Maximum Mach limited by KEAS at the given altitude.
                  
            - LiftCoefficientRange: numpy.ndarray
                  Lift coefficients over the velocity range.
            - DragCoefficientRange: numpy.ndarray
                  Drag coefficients over the velocity range.
            - DragRange: numpy.ndarray
                  Drag forces corresponding to the velocity range (same units as weight).
            - LiftRange: numpy.ndarray
                  Lift forces corresponding to the velocity range (same units as weight).
            - LDRatioRange: numpy.ndarray
                  Lift-to-drag ratios over the velocity range.
                  
            - EnduranceFactorRange: numpy.ndarray
                  Endurance factors in hours, based on lift-to-drag ratio and fuel efficiency.
            - RangeFactorRange: numpy.ndarray
                  Range factors in nautical miles, based on lift-to-drag ratio and velocity.
            - SpecificEnergyRange: numpy.ndarray
                  Specific energy levels (altitude + kinetic energy) for the velocity range (in feet).
            - SpecificExcessPowerRange: numpy.ndarray
                  Specific excess power at various velocities (in ft/s).
            - MaxVelocityAtSpecificEnergy: numpy.ndarray
                  Maximum velocities for specific energy levels.
            - MaxMachAtSpecificEnergy: numpy.ndarray
                  Maximum Mach numbers for specific energy levels.
    
        Notes
        -----
        - The function assumes the aircraft operates within a single altitude.
        - The calculations depend heavily on the atmospheric model's accuracy and interpolation methods.
        - The flight envelope is constrained by maximum Mach number (`self.max_M`) and aerodynamic limits.
        - The method includes default aerodynamic parameters (`self.CL_max`, `self.CD0`, `self.K`) and geometric values (`self.S_ref`) that must be defined in the class.

        """

        # Will only be used for singlular altitudes, not built to handle a range of h's
        flightChars = {}
        
        # Save the n in case it gets overwritten within a section
        n_def = n
        
        # Atmosphereic values
        a = atm.linterp_h(h, atm.a) # Speed of sound
        a_SL = atm.linterp_h(0, atm.a) # Speed of sound at sea level
        QMS = atm.linterp_h(h, atm.QMS)
        rho = atm.linterp_h(h, atm.rho)
        rho_SL = atm.linterp_h(0,atm.rho)
        
        max_q = atm.QMS[0] * (self.max_V_KEAS/atm.VELA[0])**2
        max_M_q = np.sqrt(max_q / QMS)
        
        T_W_at_h = (rho/rho_SL)*self.Ta_SL/self.W   # Thrust/Weight ratio at altitude
        Ta_at_h = T_W_at_h*self.W                   # Thrust available at altitude (assumes things CHECK)
        W_S      = self.W / self.S_ref              # Weight to Reference Area ratio
        
        # LHS Bound of Envelope, CL_max aka V stall
        V_stall = np.sqrt((2/rho) * (W_S) * (n/self.CL_max))
        V_stall_M = np.divide(V_stall, a)
        
        # Flight envelope (Need to adjust to handle array altitudes)
        
        # ______________ V_n Diagram __________________________________________
        # At Altitude values
        Lim_M = max_M_q if max_M_q < self.max_M else self.max_M
        
        Ms = np.arange(V_stall_M*self.V_stall_fraction, Lim_M, 0.0001) # Array of Mach numbers
        Vs = Ms * a # Array of airspeeds
        
        LoadF = self.CL_max * 0.5 * rho * np.power(Vs,2)*self.S_ref / W
        
        n_ClimbLim = np.zeros(np.size(Ms))
        n_DescendLim = np.zeros(np.size(Ms))
        
        for i, n in enumerate(LoadF):
            n_ClimbLim[i] = n if n < self.max_n else self.max_n 
            n_DescendLim[i] = -n if -n > self.min_n else self.min_n
        
        n = n_def # Reset load factor for later calculations
        
        # ____________ Thrust Required/Drag ___________________________________
        # Velocity depends on air density aka altitude, points at Tr/W min
        self.V_LD_max = np.sqrt((2/rho)*(W/self.S_ref)*np.sqrt(self.K/self.CD0)) # ft/s
        M_LD_max = np.divide(self.V_LD_max, a)
                
        
        # ____________ Mach and Lift/Drag rabges ______________________________
        M_range = Ms
        CL_range = 1*W/(self.S_ref*0.5*rho*np.power(M_range*a,2))
        CD_range = self.CD0 + self.K*np.power(CL_range,2)
        
        D_range = W*np.divide(CD_range,CL_range)
        L_range = np.multiply(D_range, np.divide(CL_range,CD_range))
        LD_range = np.divide(CL_range,CD_range)
        
        # ____________ Endurance and Range Factors ____________________________
        end_factor_range = LD_range / (self.sfc_SL * a / a_SL) # [hr]
        # CHECK: End_factor calculation
        ran_factor_range = (3600/6076.12)*np.multiply(Vs,LD_range) / (self.sfc_SL * a / a_SL) # [NM]
        
        
        # ____________ Energies _______________________________________________
        spec_energy_range = h + np.power(Vs,2) / 2*self.g # [ft], also max altitude at this energy level
        spec_excess_power_range = np.multiply(self.Ta_SL * (rho/rho_SL) - D_range, Vs) / self.W  # [ft/s]
        max_vel_at_energy = np.sqrt(spec_energy_range*2*self.g) # [ft/s] Speed attained for h-> 0 ft (all PE -> KE)
        max_mach_at_energy = max_vel_at_energy / a 
        
        # ____________ Min and Max Velocities_____ ____________________________
        min_vel_T_eq_D = np.sqrt( (T_W_at_h * W_S - W_S*np.sqrt(T_W_at_h**2 - 4*self.CD0*self.K)) / (rho*self.CD0) )
        # Vmin from stall calculated as V_stall above
        max_vel_T_eq_D = np.sqrt( (T_W_at_h * W_S + W_S*np.sqrt(T_W_at_h**2 - 4*self.CD0*self.K)) / (rho*self.CD0) )
        max_vel_at_h   = (self.max_V_KEAS/np.sqrt(rho/rho_SL))*6076.12/3600 # Max KEAS -> Max KTAS (at h) -> Max ft/s (at h)
        
        # Insert values into dictionary
        flightChars['ThrustAvailable'] = Ta_at_h
        flightChars['MachRange'] = Ms
        flightChars['VelocityRange'] = Vs
        
        flightChars['n_ClimbLimit'] = n_ClimbLim
        flightChars['n_DescendLimit'] = n_DescendLim
        flightChars['StallVelocity'] = V_stall
        flightChars['StallMach'] = V_stall_M
        flightChars['MinVelocity-T=D'] = min_vel_T_eq_D
        flightChars['MaxVelocity-T=D'] = max_vel_T_eq_D
        flightChars['MaxVelocityAtAltitude'] = max_vel_at_h
        flightChars['MaxMach-T=D'] = max_vel_T_eq_D/a 
        flightChars['MaxMachAtAltitude'] = max_vel_at_h/a 
        
        flightChars['LiftCoefficientRange'] = CL_range
        flightChars['DragCoefficientRange'] = CD_range
        flightChars['DragRange'] = D_range
        flightChars['LiftRange'] = L_range
        flightChars['LDRatioRange'] = LD_range
        
        flightChars['EnduranceFactorRange'] = end_factor_range
        flightChars['RangeFactorRange'] = ran_factor_range
        
        flightChars['SpecificEnergyRange'] = spec_energy_range
        flightChars['SpecificExcessPowerRange'] = spec_excess_power_range
        flightChars['MaxVelocityAtSpecificEnergy'] = max_vel_at_energy
        flightChars['MaxMachAtSpecificEnergy'] = max_mach_at_energy
        
        
 
        return flightChars
    
    # From HW7:
    def PlotFlightEnvelope(self, n = 1, W = None):
        '''
        Parameters
        ----------
        n : Float,
            Load Factor in g's. The default is 1.
        W : TYPE, optional
            The weight at which the envelope is calculated. The default is the Aircrafts's performance weight.

        Returns
        -------
        Plots the fight envelope. (May eventually have option to return figure)

        '''
        self.Check()
        W = self.W if W == None else W
        
        alt = np.arange(0,self.max_h+1000, 1000)
        atm = Atmosphere()
        
        a = atm.linterp_h(alt, atm.a)
        QMS = atm.linterp_h(alt, atm.QMS)
        rho = atm.linterp_h(alt, atm.rho)
        
        # LHS Bound of Envelope, CL_max aka V stall
        V_stall = np.sqrt((2/rho) * (W/self.S_ref) * (n/self.CL_max))
        V_stall_M = np.divide(V_stall, a)
        
        TopLine = np.array([[V_stall_M[-1], self.max_h], [self.max_M,self.max_h]], dtype=float)
        
        V_max_M = np.zeros(np.size(alt))
        for i, a_h in enumerate(a):
            # Based on q_max
            V_max_1 = np.sqrt((2/rho[i]) * atm.QMS[0] * (self.max_V_KEAS/atm.VELA[0])**2 )
            # Based on M max
            V_max_2 = self.max_M*a_h
            
            V_max_M[i] = V_max_1/a_h if V_max_1 < V_max_2 else V_max_2/a_h
    
        # V-n diagram - given weight and Altitude
        # V_stall_pos[i] = np.sqrt((2/rho) * (W/self.S_ref) * (n_max[i]/self.CL_max))
        plt.figure()
        plt.plot(V_stall_M,alt, TopLine[:,0],TopLine[:,1], V_max_M, alt)
        plt.title("{} Flight Envelope\nn = {} g's     W = {} lbs".format(self.Name,n,W))
        plt.grid()
        plt.xlabel('Mach Number')
        plt.ylabel('Altitude')
        
        
    def Plot_V_n_Diagram(self, h=0, W=None):
        '''
        Parameters
        ----------
        h : Float/Integer/Array,
            Altitude in feet. The default is sea level (0). If array given, all altitudes will be plotted on same plot
        W : TYPE, optional
            The weight at which the V-n diagram is calculated. The default is the Aircrafts's performance weight.

        Returns
        -------
        Plots the V-n Diagram. (May eventually have option to return figure)

        '''
        self.Check()
        W = self.W if W == None else W
        
        # Setup to handle a vector of altitudes:
        if type(h) == np.ndarray or type(h) == list:
            altitudes = h
        else:
            altitudes = [h]
        
        atm = Atmosphere()
        
        plt.figure("V_n Diagram")
        for h in altitudes:
            h = float(h)
            
            flightChars = self.Calculate_FlightCharacteristics(h,1,W,atm)
        
            Ms = flightChars['MachRange'] 
            n_ClimbLim = flightChars['n_ClimbLimit']
            n_DescendLim = flightChars['n_DescendLimit'] 
            
            # Create the Line to plot
            x_points = []
            y_points = []
            
            for i in range(0,len(Ms)):
                x_points.append(Ms[i])
                y_points.append(n_ClimbLim[i])
            
            x_points.append(Ms[-1])
            y_points.append(self.min_n)
            
            for i in range(len(Ms)-1, -1, -1):
                x_points.append(Ms[i])
                y_points.append(n_DescendLim[i])
            
            plt.plot(x_points,y_points, label='h={:.0f} ft'.format(h))
        plt.title("{} V-n Diagram\n W = {} lbs".format(self.Name,h, W))
        plt.legend()
        plt.xlabel("Mach Number")
        plt.ylabel("Load Factor (g's)")
        plt.show()
        plt.grid()
   
    
    # From HW8:
    # Min drag occurs at:
        # CD0 = K*CL^2
        # max L/D    
    def PlotThrustRequired(self, h=0, W=None, n=1):
        '''
        Parameters
        ----------
        h : Float/Integer/Array,
            Altitude in feet. The default is sea level (0). If array given, all altitudes will be plotted on same plot
        W : Float, optional
            The weight at which the V-n diagram is calculated. The default is the Aircrafts's performance weight.
        n : Float, optional
            Load Factor. The default is 1g.

        Returns
        -------
        None. Generates the Thrust Required(Drag) plot at one or a range of altitudes

        '''
        self.Check()
        #Need to change to set up handle multiple altitudes
        W = self.W if W == None else W
        
        # Calcualtions (assume all values are inputted)
        atm = Atmosphere()
        
        # Setup to handle a vector of altitudes:
        if type(h) == np.ndarray or type(h) == list:
            altitudes = h
        else:
            altitudes = [h]
        
        plt.figure("Thrust Required Diagram")
        for h in altitudes:
            h = float(h)
            
            flightChars = self.Calculate_FlightCharacteristics(h, n, W, atm)
            M_range = flightChars['MachRange']
            D_range = flightChars['DragRange'] 

            plt.plot(M_range, D_range, label="h={} ft".format(h))
        plt.title('{} Thrust Required\n W = {:.1f} lbs   n = {}g'.format(self.Name,W,n))
        plt.legend()
        plt.xlabel('Mach Number')
        plt.ylabel('Thrust Required / Drag (lbs)')
        plt.grid()
    
    def Plot_LD_vs_Mach(self, h=0, W=None, n=1):
        '''
        Plots the Lift-to-Drag Ratio (CL/CD) versus Mach number for specified altitudes.
        
        Parameters
        ----------
        h : float or list of floats, optional
            Altitude(s) (in feet) at which to calculate and plot the Lift-to-Drag ratio. The default is 0 ft.
        W : float, optional
            Aircraft weight (in lbs) used for calculations. If not specified, the default is the aircraft's performance weight.
        n : float, optional
            Load factor in g's. The default is 1.
        
        Returns
        -------
        None
            Displays a plot of Lift-to-Drag Ratio versus Mach number for the specified altitudes.
        '''

        self.Check()
        
        W = self.W if W == None else W
        
        atm = Atmosphere()
        
        # Setup to handle a vector of altitudes:
        if type(h) == np.ndarray or type(h) == list:
            altitudes = h
        else:
            altitudes = [h]
        
        plt.figure("Lift to Drag Ratio Diagram")
        for h in altitudes:
            h = float(h)
            # At Altitude values
            h = float(h)
            
            flightChars = self.Calculate_FlightCharacteristics(h, n, W, atm)
            M_range = flightChars['MachRange']
            LD_range = flightChars['LDRatioRange']
            
            plt.plot(M_range, LD_range, label="h={} ft".format(h))
        plt.title('{} Lift to Drag Ratio\n W = {:.1f} lbs  n = {}g'.format(self.Name,W,n))
        plt.legend()
        plt.xlabel('Mach Number')
        plt.ylabel('CL/CD')
        plt.grid()
        
    
    def Plot_sqrtCL_CD_vs_Mach(self, h=0, W=None, n=1):
        '''
        Plots the square root of the Lift Coefficient (sqrt(CL)) divided by the Drag Coefficient (CD) versus Mach number.
        
        Parameters
        ----------
        h : float or list of floats, optional
            Altitude(s) (in feet) at which to calculate and plot sqrt(CL)/CD. The default is 0 ft.
        W : float, optional
            Aircraft weight (in lbs) used for calculations. If not specified, the default is the aircraft's performance weight.
        n : float, optional
            Load factor in g's. The default is 1.
        
        Returns
        -------
        None
            Displays a plot of sqrt(CL)/CD versus Mach number for the specified altitudes.
        '''

        self.Check()
        
        W = self.W if W == None else W
        
        atm = Atmosphere()
        
        # Setup to handle a vector of altitudes:
        if type(h) == np.ndarray or type(h) == list:
            altitudes = h
        else:
            altitudes = [h]
        
        plt.figure("Sqrt(Lift) to Drag Ratio Diagram")
        for h in altitudes:
            h = float(h)
            # At Altitude values
            h = float(h)
            
            flightChars = self.Calculate_FlightCharacteristics(h, n, W, atm)
            M_range = flightChars['MachRange']
            CL_range = flightChars['LiftCoefficientRange']
            CD_range = flightChars['DragCoefficientRange']
            
            plt.plot(M_range, np.divide(np.sqrt(CL_range), CD_range), label="h={} ft".format(h))
        plt.title('{} \n Square root of Lift Coef by Drag Coef \n W = {:.1f} lbs  n = {}g'.format(self.Name,W,n))
        plt.legend()
        plt.xlabel('Mach Number')
        plt.ylabel(r"$\sqrt{C_L}$/C_D")
        plt.grid()
        
        
    def Plot_EndFact_vs_Mach(self, h=0, W=None, n=1):
        '''
        Plots the endurance factor versus Mach number for specified altitudes.
        
        Parameters
        ----------
        h : float or list of floats, optional
            Altitude(s) (in feet) at which to calculate and plot the endurance factor. The default is 0 ft.
        W : float, optional
            Aircraft weight (in lbs) used for calculations. If not specified, the default is the aircraft's performance weight.
        n : float, optional
            Load factor in g's. The default is 1.
        
        Returns
        -------
        None
            Displays a plot of the endurance factor versus Mach number for the specified altitudes.
        '''

        self.Check()
        W = self.W if W == None else W
       
        atm = Atmosphere()
        
        # Setup to handle a vector of altitudes:
        if type(h) == np.ndarray or type(h) == list:
            altitudes = h
        else:
            altitudes = [h]
        
        plt.figure("Endurance Factor Diagram")
        for h in altitudes:
            h = float(h)
            
            flightChars = self.Calculate_FlightCharacteristics(h, n, W, atm)
            M_range = flightChars['MachRange']
            ef_range = flightChars['EnduranceFactorRange']
            
            plt.plot(M_range, ef_range, label="h={} ft".format(h))
        plt.title('{} Endurace Factor\n W = {:.1f} lbs  n = {}g'.format(self.Name,W,n))
        plt.legend()
        plt.xlabel('Mach Number')
        plt.ylabel('Endurance Factor (hr)')
        plt.grid()
    
    
    def Plot_RangeFact_vs_Mach(self, h=0, W=None, n=1):
        '''
        Plots the range factor versus Mach number for specified altitudes.
        
        Parameters
        ----------
        h : float or list of floats, optional
            Altitude(s) (in feet) at which to calculate and plot the range factor. The default is 0 ft.
        W : float, optional
            Aircraft weight (in lbs) used for calculations. If not specified, the default is the aircraft's performance weight.
        n : float, optional
            Load factor in g's. The default is 1.
        
        Returns
        -------
        None
            Displays a plot of the range factor versus Mach number for the specified altitudes.
        '''

        self.Check()
        W = self.W if W == None else W
       
        atm = Atmosphere()
        
        # Setup to handle a vector of altitudes:
        if type(h) == np.ndarray or type(h) == list:
            altitudes = h
        else:
            altitudes = [h]
        
        plt.figure("Range Factor Diagram")
        for h in altitudes:
            h = float(h)
            
            flightChars = self.Calculate_FlightCharacteristics(h, n, W, atm)
            M_range = flightChars['MachRange']
            rf_range = flightChars['RangeFactorRange']
            
            plt.plot(M_range, rf_range, label="h={} ft".format(h))
        plt.title('{} Range Factor\n W = {:.1f} lbs  n = {}g'.format(self.Name,W,n))
        plt.legend()
        plt.xlabel('Mach Number')
        plt.ylabel('Range Factor (NM)')
        plt.grid()
     
        
     
    def Plot_Thrust_vs_Mach(self, h=0, W=None, n=1):
        '''
        Plots the range factor versus Mach number for specified altitudes.
        
        Parameters
        ----------
        h : float or list of floats, optional
            Altitude(s) (in feet) at which to calculate and plot the range factor. The default is 0 ft.
        W : float, optional
            Aircraft weight (in lbs) used for calculations. If not specified, the default is the aircraft's performance weight.
        n : float, optional
            Load factor in g's. The default is 1.
        
        Returns
        -------
        None
            Displays a plot of the Thrust versus Mach number for some specified altitudes.
        '''

        self.Check()
        W = self.W if W == None else W
       
        atm = Atmosphere()
        
        # Setup to handle a vector of altitudes:
        if type(h) == np.ndarray or type(h) == list:
            altitudes = h
        else:
            altitudes = [h]
        
        fig = plt.figure("Thrust(Drag) vs Mach Diagram")
        plt.title('{} Thrust vs Mach Number\n W = {:.1f} lbs  n = {}g'.format(self.Name,W,n))
        for i, h in enumerate(altitudes):
            h = float(h)
            
            flightChars = self.Calculate_FlightCharacteristics(h, n, W, atm)
            M_range = flightChars['MachRange']
            Thrust_req = flightChars['DragRange']
            Thrust_avail = flightChars['ThrustAvailable']
            Mach_min = flightChars['StallMach']
            Mach_max = flightChars['MaxMachAtAltitude']
            print(Mach_max)
            
            T_ax_range = [0, 1.5*Thrust_avail]
            M_ax_range = [0,np.ceil(Mach_max*10)/10]
            
            sbplt_int = str(len(altitudes))+'1'+str(i+1)
            print(sbplt_int)
            plt.subplot(len(altitudes),1,i+1)
            plt.title('h={} ft'.format(h))
            plt.plot(M_ax_range, [Thrust_avail,Thrust_avail],'k',label='_nolegend_')
            plt.plot([Mach_max,Mach_max], T_ax_range, 'k',label='_nolegend_')
            plt.plot([Mach_min,Mach_min], T_ax_range, 'k',label='_nolegend_')
            plt.plot(M_range, Thrust_req, label="Thrust Req")
            
            plt.legend()
            plt.margins(x=0, y=0)
            plt.xlabel('Mach Number')
            plt.ylabel('Thrust (Drag) [lbf]')
            plt.grid()  

        
        plt.tight_layout()
        
         
            
            
        
    #def Calc_V_Stall(self, h=0):
        
        
        
        
if __name__ == '__main__':
    # Main Testing
    atm = Atmosphere()
    # print(atm.linterp_h(0, atm.QMS))

    # _____ Wing Planiform Test ______
    pl1 = Plane()
    pl1.S_ref = 950 # ft^2 - Wing Reference area
    pl1.b = 74     # ft   - Wing Span
    pl1.LE_Sweep = 32 # deg - Sweep of leading edge - Λ
    pl1.c_t = 6.125   # ft - Tip Chord
    pl1.c_r = 20.15   # ft - Root Chord
    # pl1.WingPlaniform()

    # _____ BD-5J Plane ________
    BD_5J = Plane()
    BD_5J.Name = "BD-5J"
    BD_5J.b = 17
    BD_5J.S_ref = 37.8

    BD_5J.max_n = 6
    BD_5J.min_n = -1.5

    BD_5J.max_M = 0.55
    BD_5J.max_V_KEAS = 260
    BD_5J.max_V_fts = 300*3600/5280
    BD_5J.max_h = 26000

    BD_5J.W = 900
    BD_5J.GrossTOW = 960
    BD_5J.FuelCap = 370

    BD_5J.Ta_SL = 202
    BD_5J.sfc_SL = 1.3

    BD_5J.CL_max = 1.35
    BD_5J.CD0 = 0.02
    BD_5J.K = 0.062


    # Tests
    # BD_5J.PlotFlightEnvelope()
    # BD_5J.Plot_V_n_Diagram(0)
    # BD_5J.Plot_V_n_Diagram(np.arange(0,26000,1000))
    # BD_5J.PlotThrustRequired(np.arange(0,26000,1000))
    # print('\n\n\n')
    # BD_5J.Plot_LD_vs_Mach(np.arange(0,26000,1000))
    BD_5J.Plot_EndFact_vs_Mach(np.array([0,10000]))#(np.arange(0,26000,1000))
    BD_5J.Plot_RangeFact_vs_Mach(np.array([0,10000]))#(np.arange(0,26000,1000))
    BD_5J.Plot_Thrust_vs_Mach(h=[0,10000])

    
    