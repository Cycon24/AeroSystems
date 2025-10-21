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
import _tools.writeDictTable as wDict 
import numpy as np


# =============================================================================
# Main Class that can generate any engine 
# =============================================================================

class Engine():
    '''
    20251012 - Building this to encompase all engines to have inherited functions
                to reduce repetition and make building new engines easier. Will
                serve as the framework for developing the engine creator class.
            -  The tough thing will be enabling the ability to have multiple bypass
                ducts and handling airflow for the calculation loop in addition to 
                creating the more generalized engine performance functions.
            -  In order to make this easier, it may be necessay to only define
                certain parameters once as some performance values are tied to a 
                specific component, such as:
                    TPR: Total Pressure Recovery is tied to a single inlet, so multiple
                    would make calculation difficult.
    20251014 - Verified changes made reproduce the same values for the performance parameters 
                so all output equation are correct for Turbojet_Afterburner and Turbofan_Mixedflow_Afterburner
                    
    '''
    def __init__(self, **kwargs):
        
        self.MAIN_CALCULATION_UPDATED = False 
        
    def calculate(self, printVals: bool = True) -> None:
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
        # for i in range(0,len(self.AllStages)):
        #     # Calculate each row and print outputs
        #     self.AllStages[i][0].calculate()
        #     if printVals: self.AllStages[i][0].printOutputs()
        #     # Check if current stage has a parallel (ie, prev stage passes air to 2 stages)
        #     if self.AllStages[i][1] != None:
        #         self.AllStages[i][1].calculate()
        #         if printVals: self.AllStages[i][1].printOutputs()
                
        #     # Move forward/propogate
        #     if i != len(self.AllStages)-1: # It is not at the end, so forward
        #         if self.AllStages[i+1][1] != None: 
        #             # Means that this stage delivers to two stages: fan -> HPC & BP Noz
        #             self.AllStages[i][0].forward(self.AllStages[i+1][0],self.AllStages[i+1][1])
        #         else:
        #             # Stage delivers to one stage
        #             self.AllStages[i][0].forward(self.AllStages[i+1][0])
        
        
        return self._base_calculate(printVals)
    
    
        
    def calculateEnginePerformance(self) -> dict:
        '''
        Calculates the performance parameters for any engine infrastructure. 
        
        NOTE: 202051013 - In the future may need to update thermal efficiency
                        to include different h_PR's if different fuels are utilized.

        Returns
        -------
        Dict.
            'F_mdot': Specific Thrust (uninstalled)
            'T_mdot': Specific Thrust (isntalled)
            'F'     : Thrust (uninstalled)
            'T'     : Thrust (installed)
            'SFC'   : Thrust Specific Fuel Consumption (uninstalled)
            'TSFC'  : Thrust Specific Fuel Consumption (installed)
            'f_n'   : Fuel-Air Ratio of Combustor "n" 
            'f_tot' : FUel-Air Ratio of total engine
            'eta_T' : Thermal Efficiency
            'eta_P' : Propulsion Efficiency
            'eta_O' : Overall Efficiency
            'TPR_n' : Total Pressure Recovery of Inlet "n"
        '''
        # Get lists of each type of stage:
        Intakes     = self._getStagesOfType("Intake")
        Compressors = self._getStagesOfType("Compressor")
        Combustors  = self._getStagesOfType("Combustor")
        Turbines    = self._getStagesOfType("Turbine")
        Nozzles     = self._getStagesOfType("Nozzle")
        Mixers      = self._getStagesOfType("Mixer")
        Ducts       = self._getStagesOfType("Duct")
        
        
        # Get general props
        gam_a = self.inputs.get('Gamma_c')
        Ta = self.inputs.get('Ta')
        Pa = self.inputs.get('Pa')
        
        # from inlet stage (first in list is core inlet)
        gc = Intakes[0].gc
        gf = Intakes[0].gf
        R =  Intakes[0].R_i
        gam_a = Intakes[0].gam_i
        
        # Claculate flow velocities
        Minf = self.inputs.get('Minf')
        h_PR = self.inputs.get('h_PR')
        a0 = np.sqrt(gam_a*R*Ta*gc)
        V0 = Minf*a0
        
        # Calculate fuel-air-ratios
        fs = {}
        f_tot = 0
        for comb in Combustors:
            fs[f'f_{comb.StageID}'] = comb.f
            
            f_tot += comb.f*comb.mdot_ratio / (1 + comb.f) if not comb.IS_IDEAL else comb.f*comb.mdot_ratio
        
        # Calculate thrust chars
        noz_mom_sum = 0 
        noz_pres_sum = 0
        D_noz  = 0
        KE_rat_noz = 0 
        for noz in Nozzles:
            noz_mom_sum += noz.mdot_ratio*noz.Ve/a0
            noz_pres_sum += noz.Ae_mdota*(noz.Pe - Pa)
            
            D_noz  += 0 if noz.Drag == None else noz.Drag 
            KE_rat_noz  += noz.mdot_ratio*(noz.Ve)**2
          
        F_mdot = (a0/gc)*(noz_mom_sum - Minf) + noz_pres_sum 
        SFC = f_tot / F_mdot 
        
        # Calculate thrust if we have mass flow rate
        # Get mass flow rate
        mdot0 = Intakes[0].mdot # May need to update to incorperate more intakes
        D_inlet = 0
        TPRs = {}  # Dict of total pressure recoverie
        for inlet in Intakes: 
            D_inlet += 0 if inlet.D_additive == None else inlet.D_additive
            TPRs[f'TPR_{inlet.StageID}'] = inlet.eta_r * inlet.pi * inlet.pi_r
            # mdot0 += 0 if inlet.mdot == None else inlet.mdot # Add all massflowrate
        
        if mdot0 == None:
            F = None # Uninstalled Thrust
            T = None # Installed Thrust
            T_mdot = None # UPDATE (20251012) Can maybe calculate this if inlet/nozzle drag is calculated relative to mdot0
            TSFC = None 
        else: 
            F = F_mdot * mdot0 
            T = F - D_inlet + D_noz
            T_mdot = T / mdot0 
            TSFC = f_tot / T_mdot
            
        # Recalc SFC for units
        if self.inputs["Units"] == "SI":
            SFC *= 1e6 # go from kg/n*s * (1000g/1kg) * (1000 N / 1 kN)
            if TSFC != None: TSFC *= 1e6 
        else:
            SFC *= 3600 # go from lbm/lbf*s * (3600s/1hr)
            if TSFC != None: TSFC *= 3600
        
        
        # Calculate efficies
        eta_T = 0.5*(KE_rat_noz - V0**2) / (f_tot*h_PR*gf*gc)
        eta_P = (F_mdot)*V0 / (0.5*gc*(KE_rat_noz - V0**2))
        eta_P_inst = None if T_mdot == None else eta_P * (T_mdot/F_mdot)
        eta_O = eta_T * eta_P 
        
        
        
        # STILL INCORRECT FOR SOME REASON
        # tau_lambda = (self.Turbine.cp_e/self.Compressor.cp_e)*(self.Combustor.Toe/self.inputs['Ta'])
        # alpha_2 = ((tau_lambda/(self.Inlet.tau_r*self.Fan.tau))*(1+f_b)*(self.Turbine.tau-1) + self.Fan.tau - self.Compressor.tau/self.Fan.tau)/(1-self.Fan.tau)
        outputs = {
            'F_mdot': F_mdot,
            'T_mdot': T_mdot,
            'F': F,
            'T': T,
            'S':SFC,
            'TSFC': TSFC, 
            **fs, # unpack combustor f's
            'f_tot':f_tot,
            'eta_T':eta_T,
            'eta_P':eta_P,
            'eta_O':eta_O,
            **TPRs
            }
       
        return outputs
    
    def runParameterSweep(self, paramKey, paramList, perfFunctions=None, printVals=False) -> dict:
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
            finalOutputs[func.__name__] = self._OutputsToEmptyOutputsArray(arrayFormat, func)
        
        
        
        
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
        

    
    
    def printInputs(self, form: str ='{:9.3f}') -> None:
        '''
        Prints out all of the kwargs entered on intialization

        Returns
        -------
        None.

        '''
        print_str = "\t {}" + f"  =  {form}"
        print('Engine Inputs')
        for key,val in self.inputs.items():
            print(print_str.format(key,val))
            
        return None
    
    
    def printOutputs(self, form: str ='{:9.3f}') -> None:
        # Loop through length of all stages
        for i in range(0,len(self.AllStages)):   
            # Get stage vals for core stage
             self.AllStages[i][0].printOutputs(form)
             # Check if bypass stages exit
             if isinstance(self.AllStages[i][1], list):
                # Loop through list
                for stag in self.AllStages[i][1]:
                    stag.printOutputs(form)
                    
             # elif self.AllStages[i][1] != None:
             #     # Calculate if not list or none
             #    self.AllStages[i][1].printOutputs(form)
                
        return None
    
    def writePerformanceToFile(self, performanceOutput: dict, filepath: str = "performance.txt") -> None:
        # utilize writeDictTable to write the data to a table, this function acts like a wrapper
        # to ensure that the formatting is correct since the performance output can be any length
        # and to include correct units
        fmts = {
            # Possible variable params
            "Tt4": ".0f",
            "Tt7": ".0f",
            "Minf": "0.3f",
            
            # Thrust
            "F_mdot": ".2f",
            "T_mdot": ".2f",
            "T": ".2f",
            "F": ".2f",
            "S": ".5f",
            "TSFC": ".5f",
            # "f_tot": ".5f", # Should be handled later
            "Me": ".3f", 
            "Ve": ".2f",
            
            # Efficiencies
            "eta_T": ".2f", 
            "eta_P": ".2f", 
            "eta_O": ".2f",
            
            # Multiple Stage params (TPR_n and f_n handled later)
        }
        
        perfunits = self.getPerformanceUnits()
        
        units_map = {
            # Possible variable params
            "Tt4": self.AllStages[0][0].units_labels["T"],
            "Tt7": self.AllStages[0][0].units_labels["T"],
            
            "F_mdot": perfunits["F_mdot"],  
            "T_mdot": perfunits["F_mdot"],
            "F": perfunits["F"],
            "T": perfunits["F"],
            "S": perfunits["SFC"],
            "TSFC": perfunits["SFC"],
            "f_tot": "", # Dont need to repeat for other f's
            "eta_T": "%",
            "eta_P": "%",
            "eta_O": "%",
            "TPR": "",
            "Me": "",
            "Ve": self.AllStages[0][0].units_labels["V"],
        }
        
        # Update all efficiencies mult by 100 to convert to percentage
        for key, val in performanceOutput.items():
            # First check if its a list, if it isnt, make it one
            if not (isinstance(val, list) or isinstance(val, np.ndarray)):
                performanceOutput[key] = [val]
            
            # Value shaping
            if key.startswith("eta"):
                #I ts a list, iterate through and multiple all by 100
                for i, v in enumerate(performanceOutput[key]):
                    if v != None:
                        performanceOutput[key][i] = v*100 

            # Formatting        
            if key.startswith("f_"):
                fmts[key] = ".5f"
            if key.startswith("TPR"):
                fmts[key] = ".4f", 
  
                    

        wDict.write_dict_table(
            performanceOutput,
            file_path=filepath,
            formats=fmts,
            spacing=4,
            units=units_map,
            units_style="brackets",  # or "plain"
        )
        
        return None
    
    def getPerformanceUnits(self) -> dict:
        units = {
            "F_mdot": 'N*s\kg' if self.inputs["Units"] == "SI" else "lbf*s/lbm",
            "F": 'N' if self.inputs["Units"] == "SI" else "lbf",
            "SFC": "g/kN*s" if self.inputs["Units"] == "SI" else "lbm/lbf*hr"}
        return units
    
    def getStageStagnationVals(self) -> dict:
        '''
        Outputs a dictionary of the stagnation temperature and pressure at the 
        exit of all stages.

        Returns
        -------
        StageStagnationProps : dict
            Contains .

        '''
        StageStagnationProps = {}
        # Loop through length of all stages
        for i in range(0,len(self.AllStages)):   
            # Get stage vals for core stage
             stageOuts = self.AllStages[i][0].StageValues()['outputs']
             StageStagnationProps[self.AllStages[i][0].StageName] = {'Poe': stageOuts['Poe'], 'Toe': stageOuts['Toe']}
             # Check if bypass stages exit
             if isinstance(self.AllStages[i][1], list):
                # Loop through list
                for stag in self.AllStages[i][1]:
                    stageOuts = stag.StageValues()['outputs']
                    StageStagnationProps[stag.StageName] = {'Poe': stageOuts['Poe'], 'Toe': stageOuts['Toe']}
             # elif self.AllStages[i][1] != None:
             #     # Calculate if not list or none
             #    stageOuts = self.AllStages[i][1].StageValues()['outputs']
             #    StageStagnationProps[self.AllStages[i][1].StageName] = {'Poe': stageOuts['Poe'], 'Toe': stageOuts['Toe']}
                
        return StageStagnationProps
    
    def getStageVals(self) -> dict:
        '''
        Outputs a dictionary of all stage outputs of all stages

        Returns
        -------
        StageProps : dict
            DESCRIPTION.

        '''
        StageProps = {} 
        
        for i in range(0,len(self.AllStages)):   
             stageOuts = self.AllStages[i][0].StageValues()
             StageProps[self.AllStages[i][0].StageName] = stageOuts
             
             if isinstance(self.AllStages[i][1], list):
                 for stag in self.AllStages[i][1]:
                     stageOuts = stag.StageValues()
                     StageProps[stag.StageName] = stageOuts
             # elif self.AllStages[i][1] != None:
             #     stageOuts = self.AllStages[i][1].StageValues()
             #     StageProps[self.AllStages[i][1].StageName] = stageOuts
        return StageProps
    
    def _base_calculate(self, printVals: bool = False) -> None:
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
            if isinstance(self.AllStages[i][1], list):
                for stag in self.AllStages[i][1]: 
                    stag.calculate()
                    if printVals: stag.printOutputs()
                
            # Move forward/propogate
            if i != len(self.AllStages)-1: # It is not at the end, so forward
                # Always pass both entries, bypass handled at stage level
                self.AllStages[i][0].forward(self.AllStages[i+1][0], self.AllStages[i+1][1])
                
                    
        self.MAIN_CALCULATION_UPDATED = True 
        return None
    
    def _getStagesOfType(self, stageType: str) -> list: 
        '''
        Returns a list of all stages of type stageType

        Parameters
        ----------
        stageType : str
            Must be a string of only one of the following:
                - "Intake"
                - "Compressor"
                - "Combustor"
                - "Turbine"
                - "Nozzle"
                - "Mixer"
                - "Duct"

        Returns
        -------
        list
            The list of stages of type stageType. Can be empty if no matches found.
        '''
        stages = [] 
        
        numCoreStages = len(self.AllStages)
        # Loop through all stages list
        for i in range(0,numCoreStages):
            # Loop through items of each row
            for j in range(0, len(self.AllStages[i])):
                # Check if its a stage
                if isinstance(self.AllStages[i][j], SS.Stage):
                    # It is a stage, check if nmatches stagetype
                    if str(self.AllStages[i][j]) == stageType:
                        # Its a match, append it
                        stages.append(self.AllStages[i][j])
                    # Dont need to do anything if its not a stage
                
                # Check if its a list
                elif isinstance(self.AllStages[i][j], list):
                    # It is a list, open it up
                    for k in range(0, len(self.AllStages[i][j])):
                        # Check if its a stage
                        if isinstance(self.AllStages[i][j][k], SS.Stage):
                            # It is a stage, check if nmatches stagetype
                            if str(self.AllStages[i][j][k]) == stageType:
                                # Its a match, append it
                                stages.append(self.AllStages[i][j][k])
                        else:
                            # Not a stage, should be a stage here
                            raise ValueError('[Error]\t List within all stages contains non-stage objects')
                
                # Dont need to do anything if it is not a stage or list
        return stages 
    
    
    def _OutputsToEmptyOutputsArray(self, arrayFormat, func):
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
    

        
            
# =============================================================================
# Example Classes
# =============================================================================
class Turbofan_SingleSpool(Engine):
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
        Engine.__init__(self)
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
                          [self.HP_comp,  [self.BP_nozzle]],
                          [self.combustor,None],
                          [self.HP_turb, None],
                          [self.nozzle, None]]
        
        # ___________________ End ______________________________
        
                    
        



class Turbofan_DoubleSpool(Engine):
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
        Engine.__init__(self)
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
                          [self.HP_comp,  [self.BP_nozzle]],
                          [self.combustor,None],
                          [self.HP_turb, None],
                          [self.LP_turb, None],
                          [self.nozzle, None]]
        
        # ___________________ End ______________________________   
    
                    
    
                

class Turbofan_MixedFlow_AfterBurner(Engine):
    def __init__(self, **kwargs):
        '''
        A single spool turbofan engine with a fan with bypass, HPC, burner, HPT,
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
        Engine.__init__(self)
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
        eta_c = self.inputs.get('eta_c') # Compressor
        eta_t = self.inputs.get('eta_t') # Turbine
        eta_f = self.inputs.get('eta_f')
        
        npf = self.inputs.get('e_f') # Fan - Polytropic
        npc = self.inputs.get('e_c') # Compressor - Polytropic
        npt = self.inputs.get('e_t') # Turbine - Polytropic
        # Pressure Ratios/Relations
        pi_d = self.inputs.get('pi_d') # Diffuser total pressure ratio
        pi_b = self.inputs.get('pi_b') # Combustor total pressure ratio
        pi_f = self.inputs.get('pi_f') # Fan PR
        pi_overall = self.inputs.get('pi_c')   # Compressor PR
        self.pi_c_overall = pi_overall
        pi_c  = pi_overall if pi_f == None else pi_overall / pi_f
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
            'IS_IDEAL': self.inputs.get('IS_IDEAL'),
            # 'BPR': self.inputs.get('BPR', None)  # Must be passed as None so we can 
            # handle BPR vs pi_f calcs
            }
        BPR = self.inputs.get('BPR', None)
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
            self.Fan       = SS.Compressor(**self.gen_kwargs, pi=pi_f, np=npf, ni=eta_f, BPR=BPR)
            self.BP_duct   = SS.Duct(**self.gen_kwargs, pi=1)
            self.Compressor= SS.Compressor(**self.gen_kwargs, pi=pi_c, np=npc,ni=eta_c, pi_overall=pi_overall)
            self.Combustor = SS.Combustor(**self.gen_kwargs, Toe=To_ti, pi=pi_b, ni=eta_b, cp_e=cp_b, Gamma_e=gam_b)
            self.Turbine   = SS.Turbine([self.Compressor, self.Fan], **self.gen_kwargs, nm=eta_m, np=npt)
            self.Mixer     = SS.Mixer(self.Turbine, self.BP_duct, **self.gen_kwargs, pi=pi_M, MIX_GAS_PROPERTIES=self.inputs.get('MIX_GAS_PROPERTIES'))
            self.Afterburner = SS.Combustor(**self.gen_kwargs, pi=pi_AB, ni=eta_ab, dTb=dTo_ab, cp_e=cp_ab, Gamma_e=gam_ab)
            self.Nozzle    = SS.Nozzle(air_type='hot',nozzle_type='CD',**self.gen_kwargs, pi=pi_n) 
            
            # Set names for easier readout checks
            self.Fan.UpdateInputs(StageName = 'Fan', StageID="f")
            self.Afterburner.UpdateInputs(StageName = 'Afterburner', StageID="ab")
            
            # Set other conditions
            self.Combustor.Toe = To_ti # Set combustor outlet temperature
            self.Turbine.Me = M6
            self.Afterburner.Toe = To_ab_e 
            
            
            # Define all stages in engine to iterate through
            # Two dimensional since there is a bypass, ie one stage
            # passes params to two different stages
            self.AllStages = [[self.Inlet, None ],
                              [self.Fan, None], 
                              [self.Compressor,  [self.BP_duct]],
                              [self.Combustor,None],
                              [self.Turbine, None],
                              [self.Mixer, None],
                              [self.Afterburner, None],
                              [self.Nozzle, None]]
            
        else:
            # Only pass values into each stage rather than all
             self.Inlet      .UpdateInputs(**self.gen_kwargs, m_dot=mdot, pi=pi_d)
             self.Fan        .UpdateInputs(**self.gen_kwargs, pi=pi_f, np=npf, BPR=BPR)
             self.BP_duct    .UpdateInputs(pi=1)
             self.Compressor .UpdateInputs(**self.gen_kwargs, pi=pi_c, np=npc,pi_overall=pi_overall)
             self.Combustor  .UpdateInputs(**self.gen_kwargs, Toe=To_ti, pi=pi_b, ni=eta_b)
             self.Turbine    .UpdateInputs([self.Compressor, self.Fan], **self.gen_kwargs, nm=eta_m, np=npt)
             self.Mixer      .UpdateInputs(self.Turbine, self.BP_duct, **self.gen_kwargs, pi=pi_M)
             self.Afterburner.UpdateInputs(**self.gen_kwargs, pi=pi_AB, ni=eta_ab, dTo=dTo_ab)
             self.Nozzle     .UpdateInputs(**self.gen_kwargs, pi=pi_n) 
             
             # Set other conditions
             self.Combustor.Toe = To_ti # Set combustor outlet temperature
             self.Turbine.Me = M6
             self.Afterburner.Toe = To_ab_e 
        
             
    def calculate(self, printVals: bool = True) -> None:
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
        iterations = 0
        if hasattr(self, 'solverType'):
            self.inputs[self.solverType] = None
        
        if self.inputs["BPR"] == None and self.inputs["pi_f"] != None:
            print("[Note]\t Fan Pressure Ratio given, iterative solving for BPR")
            self.solverType = "BPR"
            # Need to iterate to solve for BPR
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
           
            error = 0
            dPt = 1
            tol = 1e-4
            while abs(dPt) > tol:
                iterations += 1
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
                
                self._base_calculate()
                
                # Check error
                if abs(alpha) < tol:
                    dPt = 0 
                    error = 0 
                    # print('Warning: Bypass Ratio reached 0...')
                else:
                    dPt = self.Turbine.Poe - self.BP_duct.Poe 
                    error = dPt/(self.Turbine.Poe + self.BP_duct.Poe) #/self.Turbine.Poe 
                
                
                
                
                
        elif self.inputs["BPR"] != None and self.inputs["pi_f"] == None:
            print("[Note]\t BPR given, iterative solving for Fan Pressure Ratio") 
            self.solverType = "pi_f"
            # Start by assuming a pressure ratio of 1
            pi_f = 1 
            
            # Set up error for loop
            error = 0
            dPt = 1
            tol = 1e-4

            # Enter loop
            while abs(dPt) > tol and dPt != None:
                # print(f"Loop iteration {i}")
                pi_f += error 
                iterations +=1
                # print(f"pi_f = {pi_f:.4f}, error = {error:.4f}")
                # check if alpha went < 0
                if pi_f < 1:
                    print('[Info]\t pi_f set to 1')
                    pi_f = 1

                # print(' f = {}\n alpha = {}\n alpha_f = {}'.format(f, alpha,alpha_f))
                
                # Update
                # pi_c  = self.pi_c_overall if pi_f == None else self.pi_c_overall / pi_f
                self.UpdateInputs(pi_f=pi_f)
                
                
                # Run Calculations
                self._base_calculate()
                
   
                # Check error
                if abs(pi_f) < 1-tol:
                    dPt = 0 
                    error = 0 
                    print('[Warning]\t Fan Pressure Ratio reached 1...')
                else:
                    dPt = self.Turbine.Poe - self.BP_duct.Poe 
                    error = dPt/(self.Turbine.Poe + self.BP_duct.Poe) #/self.Turbine.Poe 
                    # print(f'dPt = {dPt:.4f}, P_t4 {self.Turbine.Poe:.2f} P_tbp {self.BP_duct.Poe:.2f}')
                
    
        print(f"[Info]\t Iterative solver solution for {self.solverType}  = {self.inputs[self.solverType]:.4f} in {iterations} iterations")  

        if printVals:
            for i in range(0,len(self.AllStages)):   
                 self.AllStages[i][0].printOutputs()
                 if isinstance(self.AllStages[i][1], list):
                     for stag in self.AllStages[i][1]:
                         stag.printOutputs()
        return None
            
    
                                  
               
    
    
   
    # def CalculatePerformanceParams(self):
    #     '''
    #     Need to get: 
    #     F/mdot_0
    #     SFC = f_tot / F/mdot0
    #     f
    #     f_AB
    #     f_0 = f_b / (1+ alpha) + f_ab 
    #     eta_T = 1 - 1 / (tau_r*tau_c)
    #     eta_P = 2*Minf * (V9/a0 - Minf) / ((V9/a0)^2 - M0^2) 
    #     eta_O = T*P
    #     alpha

    #     Returns
    #     -------
    #     Dict.
    #         'F_mdot': Specific Thrust
    #         'SFC': Thrust Specific Fuel Consumption
    #         'f': Fuel-Air Ratio of Combustor
    #         'f_AB': Fuel-Air Ratio of Afterburner
    #         'f_tot': FUel-Air Ratio of total engine
    #         'eta_T': Thermal Efficiency
    #         'eta_P': Propulsion Efficiency
    #         'eta_O': Overall Efficiency
    #         'alpha': Bypass Ratio (mdot_bypass / mdot_core)

    #     '''
    #     # Thrust, since Pe = Pa no need to incorperate 
    #     gam_a = self.inputs.get('Gamma_c')
    #     gc = self.Inlet.gc
    #     gf = self.Inlet.gf
    #     R = self.Inlet.R_i
    #     Ta = self.inputs.get('Ta')
    #     # print(f"{gam_a}*{R}*{Ta}*{gc}")
    #     a0 = np.sqrt(gam_a*R*Ta*gc)
    #     Minf = self.inputs.get('Minf')
    #     V0 = Minf*a0
    #     h_PR = self.inputs.get('h_PR')
        
    #     F_mdot = (a0/gc)*(self.Nozzle.mdot_ratio*self.Nozzle.Ve/a0 - Minf)
    #     f_b = self.Combustor.f
    #     f_ab = self.Afterburner.f
    #     alpha = self.Fan.BPR 
    #     f_tot = f_b/(1 + alpha) + f_ab
    #     SFC = 3600*f_tot / F_mdot
        

    #     mdot_N = self.Inlet.m_dot
    #     Thrust = None if self.Nozzle.m_dot == None else F_mdot * mdot_N
    #     Installed_Thrust =  None if self.Inlet.D_additive == None or Thrust == None else Thrust - self.Inlet.D_additive
    #     T_mdot = None if Installed_Thrust == None else Installed_Thrust / self.Inlet.m_dot
    #     mdot_f_tot =  None if self.Inlet.m_dot == None else f_tot * self.Inlet.m_dot
    #     TSFC =  None if Installed_Thrust == None else 3600*mdot_f_tot / Installed_Thrust
        
    #     eta_T = 0.5*(self.Nozzle.mdot_ratio*self.Nozzle.Ve**2 - V0**2) / (f_tot*h_PR*gf*gc)
    #     #1 - 1 / (self.Inlet.tau_r * self.Compressor.tau)
    #     eta_P = 2*Minf * (self.Nozzle.Ve/a0 - Minf) / ((self.Nozzle.Ve/a0)**2 - Minf**2)
    #     eta_O = eta_T * eta_P 
        
    #     TPR = self.Inlet.eta_r * self.Inlet.pi * self.Inlet.pi_r
        
    #     # STILL INCORRECT FOR SOME REASON
    #     # tau_lambda = (self.Turbine.cp_e/self.Compressor.cp_e)*(self.Combustor.Toe/self.inputs['Ta'])
    #     # alpha_2 = ((tau_lambda/(self.Inlet.tau_r*self.Fan.tau))*(1+f_b)*(self.Turbine.tau-1) + self.Fan.tau - self.Compressor.tau/self.Fan.tau)/(1-self.Fan.tau)
    #     outputs = {
    #         'F_mdot': F_mdot,
    #         'T_mdot': T_mdot,
    #         'F': Thrust,
    #         'T': Installed_Thrust,
    #         'S':SFC,
    #         'TSFC': TSFC, 
    #         'f_b':f_b,
    #         'f_ab':f_ab,
    #         'f_tot':f_tot,
    #         'eta_T':eta_T,
    #         'eta_P':eta_P,
    #         'eta_O':eta_O,
    #         'alpha':alpha,
    #         'TPR': TPR}
    #     # print('M0 = {:.2f}\tα1 = {:.3f}\tα2 = {:.3f}'.format(Minf, alpha, alpha_2))
    #     return outputs
    
    # def getStageStagnationVals(self):
    #     StageStagnationProps = {}
    #     for i in range(0,len(self.AllStages)):   
    #          stageOuts = self.AllStages[i][0].StageValues()['outputs']
    #          StageStagnationProps[self.AllStages[i][0].StageName] = {'Poe': stageOuts['Poe'], 'Toe': stageOuts['Toe']}
    #          if self.AllStages[i][1] != None:
    #              stageOuts = self.AllStages[i][1].StageValues()['outputs']
    #              StageStagnationProps[self.AllStages[i][1].StageName] = {'Poe': stageOuts['Poe'], 'Toe': stageOuts['Toe']}
    #     return StageStagnationProps
    
    # def getStageVals(self):
    #     StageProps = {} 
        
    #     for i in range(0,len(self.AllStages)):   
    #          stageOuts = self.AllStages[i][0].StageValues()
    #          StageProps[self.AllStages[i][0].StageName] = stageOuts
             
    #          if self.AllStages[i][1] != None:
    #              stageOuts = self.AllStages[i][1].StageValues()
    #              StageProps[self.AllStages[i][1].StageName] = stageOuts
    #     return StageProps
    
    
    
class Turbojet_Afterburner(Engine):
    def __init__(self, **kwargs):
        '''
        A single spool turbojet engine with diffuser, HPC, burner, HPT,
        afterburner and nozzle. 
        Parameters
        ----------
        **kwargs : Dictionary
            Contains all needed and optional parameters with the keys listed below.
            Required:
            'Ta': Atmospheric static temperature
            'Pa': Atmospheric static pressure
            'pi_c':   Compressor Pressure Ratio
            'T_turb_in': Turbine inlet temp (HP Turbine)
            
            Optional
            'Vinf':  Or Minf
            'Minf':  Or Vinf, if none its assumed stationary
            'mdot_a': Mass flow rate of air into engine (kg/s) (lbm/s)
            'Q_fuel': Heat energy of fuel, pass if real to calculate f and include fuel flow in power/velecity calcs
            
            Efficiencies (Assumed to be 1 if not passed)
            'ni': Inlet Isentropic Efficiency
            'nj': Nozzle Isentropic Efficiency
            'nc': Compressor Isentropic Efficiency
            'nt': Turbine Isentropic Efficiency
            'nb': Cobustor Efficincy
            'nm': Mechanical Efficiency
            'npc': Compressor Polytropic Efficiency (overrides isentropic)
            'npt': Turbine Polytropic Efficiency (overrides isentropic)
            
            'pi_b': Total pressure ratio through the burner
            'pi_ab': Total pressure ratio through after burner

        Returns
        -------
        None.

        '''
        # Stages
        # Atm moving
        # Inlet
        # HP Compressor
        # Combustor
        # HP Turbine
        # Afterburner
        # Nozzle
        Engine.__init__(self)
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
        eta_c = self.inputs.get('eta_c') # Compressor
        eta_t = self.inputs.get('eta_t') # Turbine
      
        npc = self.inputs.get('e_c') # Compressor - Polytropic
        npt = self.inputs.get('e_t') # Turbine - Polytropic
        # Pressure Ratios/Relations
        pi_d = self.inputs.get('pi_d') # Diffuser total pressure ratio
        pi_b = self.inputs.get('pi_b') # Combustor total pressure ratio

        pi_overall = self.inputs.get('pi_c')   # Compressor PR
        pi_c  = pi_overall
        
        pi_n  = self.inputs.get('pi_n') # Nozzle total pressure ratio
        pi_AB = self.inputs.get('pi_ab') # Afterburner total ressure raito
        pi_AB_off = self.inputs.get('pi_ab_off', pi_AB) # AB total pressure ratio without operation
        
        # Turbine Inlet / Combustor Outlet
        To_ti = self.inputs.get('Tt4') # K - Turbine inlet temp
        To_ab_e = self.inputs.get('Tt7')
        dTo_ab  = self.inputs.get('dTt_ab', None) # Change in total temp of Afterburner
        
        # Air Mass flow
        mdot = self.inputs.get('mdot_a') # kg/s or lbm/s
        
        
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
            self.Inlet     = SS.Intake(**self.gen_kwargs, mdot=mdot, pi=pi_d, cp_i=cp_air, Gamma_i=gam_air)
            self.Compressor   = SS.Compressor(**self.gen_kwargs, pi=pi_c, ni=eta_c,np=npc, pi_overall=pi_overall)
            self.Combustor = SS.Combustor(**self.gen_kwargs, Toe=To_ti, pi=pi_b, ni=eta_b, cp_e=cp_b, Gamma_e=gam_b)
            self.Turbine   = SS.Turbine(self.Compressor, **self.gen_kwargs, nm=eta_m, ni=eta_t, np=npt)
            self.Afterburner = SS.Combustor(**self.gen_kwargs, pi=pi_AB, pi_off=pi_AB_off, ni=eta_ab, dTb=dTo_ab, cp_e=cp_ab, Gamma_e=gam_ab)
            self.Nozzle    = SS.Nozzle(air_type='hot',nozzle_type='CD',**self.gen_kwargs, pi=pi_n) 
            
            # Set names for easier readout checks
            self.Afterburner.UpdateInputs(StageName = 'Afterburner', StageID="ab")
            
            # Set other conditions
            self.Combustor.Toe = To_ti # Set combustor outlet temperature 
            self.Afterburner.Toe = To_ab_e # if delta is used, its already passed into afterburner definition
            
            
            # Define all stages in engine to iterate through
            # Two dimensional since there is a bypass, ie one stage
            # passes params to two different stages
            self.AllStages = [[self.Inlet, None ],
                              [self.Compressor,  None],
                              [self.Combustor,None],
                              [self.Turbine, None],
                              [self.Afterburner, None],
                              [self.Nozzle, None]]
            
        else:
            # Only pass values into each stage rather than all
             self.Inlet      .UpdateInputs(**self.gen_kwargs, m_dot=mdot, pi=pi_d)
             self.Compressor .UpdateInputs(**self.gen_kwargs, pi=pi_c, np=npc,pi_overall=pi_overall)
             self.Combustor  .UpdateInputs(**self.gen_kwargs, Toe=To_ti, pi=pi_b, ni=eta_b)
             self.Turbine    .UpdateInputs(self.Compressor, **self.gen_kwargs, nm=eta_m, np=npt)
             self.Afterburner.UpdateInputs(**self.gen_kwargs, pi=pi_AB, ni=eta_ab)
             self.Nozzle     .UpdateInputs(**self.gen_kwargs, pi=pi_n) 
             
             # Set other conditions
             self.Combustor.Toe = To_ti # Set combustor outlet temperature
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
        # Since there isnt a need to iterate for this engine type, we can just calculate outright
        # Get gas props
                   
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
            

    
        if printVals:
            for i in range(0,len(self.AllStages)):   
                 self.AllStages[i][0].printOutputs()
                 if self.AllStages[i][1] != None:
                     self.AllStages[i][1].printOutputs()
            
    # def RunParameterSweep(self, paramKey, paramList, perfFunctions=None, printVals=False):
    #     # Setup the output dictionary
    #     # Needs to contain an item for each perfFunction
    #     # Then the values within each perfFunction needs to be lists. 
    #     # Lengths of lists based on param list, this needs to be defined
    #     #   prior to running sweep
        
    #     arrayFormat = np.zeros((len(paramList),))
    #     # Make output dic
    #     finalOutputs = {paramKey: paramList}
        
    #     # If perf functions arent a list, make them one
    #     if type(perfFunctions) != list:
    #         perfFunctions = [perfFunctions]
    #     for func in perfFunctions:
    #         finalOutputs[func.__name__] = self._OutputsToEmptyOutputsArray(arrayFormat, func)
        
        
        
    #     for i, p in enumerate(paramList):
    #         # Update inputs with parameter
    #         # self.inputs[paramKey] = p
    #         self.UpdateInputs(**{paramKey: p})
            
    #         # Run Calculation Loop
    #         self.calculate(printVals=printVals)
            
    #         # Run output functions
    #         if perfFunctions != None:
    #             # Iterate through functions                  
    #             for func in perfFunctions:
    #                 # Get function outputs
    #                 func_outs = func()
    #                 # Run through each thing in outputs
    #                 for key in func_outs.keys():
    #                     if type(func_outs[key]) == dict:
    #                         # func outputs nested dict (assuming only 3 layers)
    #                         for key2 in func_outs[key].keys():
    #                             if type(func_outs[key][key2]) == dict:
    #                                 for key3 in func_outs[key][key2].keys():
    #                                     if type(func_outs[key][key2][key3]) != str:
    #                                         finalOutputs[func.__name__][key][key2][key3][i] = func_outs[key][key2][key3]
    #                                     else:
    #                                         finalOutputs[func.__name__][key][key2][key3] = func_outs[key][key2][key3]
                                
    #                             elif type(func_outs[key][key2]) != str:
    #                                 # Set an empty list as the value
    #                                 finalOutputs[func.__name__][key][key2][i] = func_outs[key][key2]
    #                             else:
    #                                 finalOutputs[func.__name__][key][key2] = func_outs[key][key2]
    #                     elif type(func_outs[key]) != str:
    #                         finalOutputs[func.__name__][key][i] = func_outs[key]
    #                     else:
    #                         # Is a string
    #                         finalOutputs[func.__name__][key] = func_outs[key]
                        
                                   
    #     return finalOutputs
        

    # def _OutputsToEmptyOutputsArray(self, arrayFormat, func):
    #     '''
    #     Returns a dictionary containing lists in the same format as
    #     the function 'func' outputs (which outputs a dict with single values)

    #     Parameters
    #     ----------
    #     arrayFormat : np Array
    #         An empty numpy array the same size as needed for param sweep.
    #     func : Function
    #         The performance/output function.

    #     Returns
    #     -------
    #     funcOuts : Dict
    #         Contains dict with values as lists.

    #     '''
    #     # Run the function
    #     self.calculate(False)
    #     funcOuts = func()
    #     for key in funcOuts.keys():
    #         if type(funcOuts[key]) == dict:
    #             # func outputs nested dict (assuming only 3 layers)
    #             for key2 in funcOuts[key].keys():
    #                 if type(funcOuts[key][key2]) == dict:
    #                     for key3 in funcOuts[key][key2].keys():
    #                         if type(funcOuts[key][key2][key3]) != str:
    #                             funcOuts[key][key2][key3] = arrayFormat.copy()
                            
    #                 elif type(funcOuts[key][key2]) != str:
    #                     # Set an empty list as the value
    #                     funcOuts[key][key2] = arrayFormat.copy()
                        
    #         elif type(funcOuts[key]) != str:
    #             funcOuts[key] = arrayFormat.copy()
                
    #     return funcOuts
                                  
               
    
    
            
   
    # def CalculatePerformanceParams(self):
    #     return self.gen_CalculatePerformanceParams()
    #     '''
    #     Need to get: 
    #     F/mdot_0
    #     SFC = f_tot / F/mdot0
    #     f
    #     f_AB
    #     f_0 = f_b / (1+ alpha) + f_ab 
    #     eta_T = 1 - 1 / (tau_r*tau_c)
    #     eta_P = 2*Minf * (V9/a0 - Minf) / ((V9/a0)^2 - M0^2) 
    #     eta_O = T*P
    #     alpha

    #     Returns
    #     -------
    #     Dict.
    #         'F_mdot': Specific Thrust
    #         'SFC': Thrust Specific Fuel Consumption
    #         'f': Fuel-Air Ratio of Combustor
    #         'f_AB': Fuel-Air Ratio of Afterburner
    #         'f_tot': FUel-Air Ratio of total engine
    #         'eta_T': Thermal Efficiency
    #         'eta_P': Propulsion Efficiency
    #         'eta_O': Overall Efficiency
    #         'alpha': Bypass Ratio (mdot_bypass / mdot_core)

    #     '''
    #     # Thrust, since Pe = Pa no need to incorperate 
    #     gam_a = self.inputs.get('Gamma_c')
    #     gc = self.Inlet.gc
    #     gf = self.Inlet.gf
    #     R = self.Inlet.R_i
    #     Ta = self.inputs.get('Ta')
    #     a0 = np.sqrt(gam_a*R*Ta*gc)
    #     Minf = self.inputs.get('Minf')
    #     V0 = Minf*a0
    #     h_PR = self.inputs.get('h_PR')
        
    #     F_mdot = (a0/gc)*(self.Nozzle.mdot_ratio*self.Nozzle.Ve/a0 - Minf)
    #     f_b = self.Combustor.f
    #     f_ab = self.Afterburner.f
    #     alpha = 0 
    #     f_tot = f_b/(1 + alpha) + f_ab
    #     SFC = 3600*f_tot / F_mdot # Metric: kg / (N*s)
        

    #     mdot_N = self.Inlet.m_dot
    #     Thrust = None if self.Nozzle.m_dot == None else F_mdot * mdot_N
    #     Installed_Thrust =  None if self.Inlet.D_additive == None or Thrust == None else Thrust - self.Inlet.D_additive
    #     T_mdot = None if Installed_Thrust == None else Installed_Thrust / self.Inlet.m_dot
    #     mdot_f_tot =  None if self.Inlet.m_dot == None else f_tot * self.Inlet.m_dot
    #     TSFC =  None if Installed_Thrust == None else 3600*mdot_f_tot / Installed_Thrust
        
    #     eta_T = 0.5*(self.Nozzle.mdot_ratio*self.Nozzle.Ve**2 - V0**2) / (f_tot*h_PR*gf*gc)
    #     #1 - 1 / (self.Inlet.tau_r * self.Compressor.tau)
    #     eta_P = 2*Minf * (self.Nozzle.Ve/a0 - Minf) / ((self.Nozzle.Ve/a0)**2 - Minf**2)
    #     eta_O = eta_T * eta_P 
        
    #     TPR = self.Inlet.eta_r * self.Inlet.pi * self.Inlet.pi_r
        
    #     # STILL INCORRECT FOR SOME REASON
    #     # tau_lambda = (self.Turbine.cp_e/self.Compressor.cp_e)*(self.Combustor.Toe/self.inputs['Ta'])
    #     # alpha_2 = ((tau_lambda/(self.Inlet.tau_r*self.Fan.tau))*(1+f_b)*(self.Turbine.tau-1) + self.Fan.tau - self.Compressor.tau/self.Fan.tau)/(1-self.Fan.tau)
    #     outputs = {
    #         'F_mdot': F_mdot,
    #         'T_mdot': T_mdot,
    #         'F': Thrust,
    #         'T': Installed_Thrust,
    #         'S':SFC,
    #         'TSFC': TSFC, 
    #         'f_b':f_b,
    #         'f_ab':f_ab,
    #         'f_tot':f_tot,
    #         'eta_T':eta_T,
    #         'eta_P':eta_P,
    #         'eta_O':eta_O,
    #         'alpha':alpha,
    #         'TPR': TPR}
    #     # print('M0 = {:.2f}\tα1 = {:.3f}\tα2 = {:.3f}'.format(Minf, alpha, alpha_2))
    #     return outputs
    
    
    
    
# =============================================================================
#  Main Testing
# =============================================================================
if __name__=="__main__":
    eng = Turbojet_Afterburner()
    print(eng._getStagesOfType("Nozzle"))
    # Works
    

