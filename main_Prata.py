"""
File: main.py
Authors: Bruno Dias and Jeremie Meurisse
Date: 06/27/24
Description: main script to run the 0D reactor with different reaction rates.
"""

# This is a sample Python script.
import mutationpp
# Press ⌃R to execute it or replace it with your code.
# Press Double ⇧ to search everywhere for classes, files, tool windows, actions, and settings.
import mutationpp as mpp

import OxidationModel
from OxidationModel import OxidationRateModel
from OxidationModel import packages as pk
import matplotlib.pyplot as plt
import numpy as np

def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
    print(f'Hi, {name}')  # Press ⌘F8 to toggle the breakpoint.
    # help(mpp.Mixture)
    
    Tw = np.linspace(800, 2400, 100)
    p_CO = []
    p_O = []
    p_O2 = []
    p_CO2 = []
    
    f_CO = []
    f_O = []
    f_O2 = []
    f_CO2 = []
    
    for T in Tw:
    
        T_beam = 1000  # beam temperature [K]
        p_beam = 2.4e-2  # beam pressure [Pa]
        
        oxidation = OxidationRateModel()
        oxidation.pressure_in = p_beam
        oxidation.temperature_in = T_beam
        
        x_in = {"O": 1, "O2": 0}  # mole fraction of inflow gas: x_O or x_O2 [-]
        
        oxidation.composition_in = x_in

        # oxidation_bruno = oxidation.oxidation_model("Bruno_oxidation")
        oxidation_prata = oxidation.oxidation_model("prata_oxidation")
        oxidation_prata.wall_temperature =T
        p_i, f_i = oxidation_prata.solve_ODEs(T)
        p_CO.append(p_i["CO"])
        p_O.append(p_i["O"])
        p_O2.append(p_i["O2"])
        p_CO2.append(p_i["CO2"])
        
        f_CO.append(f_i["CO"])
        f_O.append(f_i["O"])
        f_O2.append(f_i["O2"])
        f_CO2.append(f_i["CO2"])

        # oxidation_prata.plot_reaction_rates()
        # oxidation_prata.plot_model_prediction()
        # oxidation_prata.plot_surface_coverage()
        
        print(oxidation_prata.wall_temperature)
        # oxidation_za = OxidationRateModel("ZA_oxidation")
    p_CO_Prata = np.loadtxt("data/CO_probability.dat", delimiter=',')
    p_CO2_Prata = np.loadtxt("data/CO2_probability.dat", delimiter=',')
    p_O2_Prata = np.loadtxt("data/O2_probability.dat", delimiter=',')
    p_O_Prata = np.loadtxt("data/O_probability.dat", delimiter=',')

    plt.figure(1)
    plt.plot(p_CO_Prata[:, 0], p_CO_Prata[:, 1], 'k*', label='CO Prata')
    plt.plot(p_CO2_Prata[:, 0], p_CO2_Prata[:, 1], 'r*', label='CO2 Prata')
    plt.plot(p_O2_Prata[:, 0], p_O2_Prata[:, 1], 'g*', label='O2 Prata')
    plt.plot(p_O_Prata[:, 0], p_O_Prata[:, 1], 'b*', label='O Prata')

    plt.plot(Tw, np.array(p_CO), 'k', label='CO')
    plt.plot(Tw, np.array(p_CO2), 'r', label='CO2')
    plt.plot(Tw, np.array(p_O2), 'g', label='O2')
    plt.plot(Tw, np.array(p_O), 'b', label='O')
    plt.legend(loc='best')
    plt.xlabel('T [K]')
    plt.xlim(800, 2400)
    plt.ylim(-0.1, 1)
    plt.grid()

    
    plt.figure(2)
    plt.plot(Tw, np.array(f_CO), 'k', label='CO')
    plt.plot(Tw, np.array(f_CO2), 'r', label='CO2')
    plt.plot(Tw, np.array(f_O2), 'g', label='O2')
    plt.plot(Tw, np.array(f_O), 'b', label='O')
    plt.legend(loc='best')
    plt.xlabel('T [K]')
    plt.xlim(800, 2400)
    plt.grid()
    
    plt.show()    


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    print_hi('PyCharm')

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
