# -*- coding: utf-8 -*-
"""
Created on Fri Nov 17 21:18:14 2023

@author: cycon
"""
def get_Fuel_Table():
    '''
    Returns a dictionary of fuel data, including their n number, formula, common name,
    formation enthalpy, and LHV.

    Returns
    -------
    Fuel_Table : dict
        Dictionary with array values according to the keys. Index of 'n' array corresponds to that molecule's properties.

    '''
    Fuel_Table = {
        'n': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
        'Formula': ['CH4', 'C2H6', 'C3H8', 'C4H10', 'C5H12', 'C6H14', 'C7H16', 'C8H18', 'C9H20', 'C10H22', 'C11H24', 'C12H26'],
        'Common name': ['Methane', 'Ethane', 'Propane', 'Butane', 'Pentane', 'Hexane', 'Heptane', 'Octane', 'Nonane', 'Decane', 'Undecane', 'Dodecane'],
        'hf_o (kJ/mol)': [-74850, -84680, -103850, -126150, -147100, -167200, -187900, -208450, -229300, -249400, -270300, -291010],
        'LHV (kJ/kg)': [50050, 47520, 46340, 45370, 45350, 44752, 44566, 44430, 44311, 44240, 44194, 44147]
    }
    return Fuel_Table


def get_Gas_Enthalpies():
    '''
    Returns a dictionary of the formation enthalpies and enthalpy at 298K for CO2, H2O, N2, and O2

    Returns
    -------
    Gas_Enthalpies : Dict
        Contains the items Molecule, hf_o (kJ/mol), and h_298K (kJ/mol) with list data

    '''
    Gas_Enthalpies = {
        'Molecule': ['CO2', 'H2O', 'N2', 'O2'],
        'hf_o (kJ/mol)': [-393520, -241820, 0, 0],
        'h_298K (kJ/mol)': [9364, 9904, 8669, 8682]
    }
    return Gas_Enthalpies
