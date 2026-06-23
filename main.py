"""
File: main.py
Authors: Bruno Dias and Jeremie Meurisse
Date: 06/27/24
Description: main script to run the 0D reactor with different reaction rates.
"""

# Press ⌃R to execute it or replace it with your code.
# Press Double ⇧ to search everywhere for classes, files, tool windows, actions, and settings.
import mutationpp as mpp

import OxidationModel
from OxidationModel import OxidationRateModel
from OxidationModel import packages as pk
import matplotlib.pyplot as plt
import numpy as np
from labellines import labelLines


SMALL_SIZE = 12
MEDIUM_SIZE = 14
BIGGER_SIZE = 18


plt.rc('axes', titlesize=MEDIUM_SIZE)  # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)  # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)  # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)  # fontsize of the tick labels
plt.rcParams['lines.linewidth'] = 2
# cc = cycler('linestyle', ['-', '--', ':', '-.', '-']) * \
#     cycler('color', ['r', 'g', 'b', 'm', 'k'])
# plt.rc('axes', prop_cycle=cc)

# plt.rc('text', usetex=True)
plt.rc('font', family='Helvetica')

def zeroD_reactor():
    
    
    experimental_presure = np.array(
                    [2172,
                     2922,
                     3546,
                     3673,
                     3668,
                     4079,
                     5325,
                     6180])

    experimental_temperature = np.array(
                    [518,
                     718,
                     818,
                     928,
                     1085,
                     1116,
                     1314,
                     1502])

    labels =[   "A",
                "B",
                "C",
                "D",
                "E",
                "F",
                "G",
                "H"]

    labels_O2 = ['O2_'+i for i in labels]   
    
    p_CO = []
    p_O = []
    p_O2 = []
    p_CO2 = []
    
    f_CO = []
    f_O = []
    f_O2 = []
    f_CO2 = []
    
    f_CO_Prata = []
    f_O_Prata = []
    f_O2_Prata = []
    f_CO2_Prata = []
    
    f_CO_ZA = []
    f_O_ZA = []
    f_O2_ZA = []
    f_CO2_ZA = []
    
    
    for T,P in zip(experimental_temperature, experimental_presure):
    
        T_beam = T  # beam temperature [K]
        p_beam = P  # beam pressure [Pa]
        
        print("Temp", T, "Press",P)
        
        oxidation = OxidationRateModel()
        oxidation.pressure_in = p_beam
        oxidation.temperature_in = T_beam
        
        x_in = {"O": 0, "O2": 10}  # mole fraction of inflow gas: x_O or x_O2 [-]
        
        oxidation.composition_in = x_in

        oxidation_bruno = oxidation.oxidation_model("Bruno_oxidation")
        oxidation_prata = oxidation.oxidation_model("prata_oxidation")
        oxidation_ZA= oxidation.oxidation_model("za_oxidation")
        
        oxidation_bruno.wall_temperature =T
        f_i = oxidation_bruno.solve_ODEs(T)
        
        f_CO.append(f_i["CO"])
        f_O2.append(f_i["O2"])
        f_CO2.append(f_i["CO2"])

        oxidation_prata.wall_temperature =T
        p_i, f_i_Prata = oxidation_prata.solve_ODEs(T)
        
        f_CO_Prata.append(f_i_Prata["CO"])
        f_O2_Prata.append(f_i_Prata["O2"])
        f_CO2_Prata.append(f_i_Prata["CO2"])
        
        oxidation_ZA.wall_temperature = T
        f_i_ZA = oxidation_ZA.solve_ODEs(T)
        
        f_CO_ZA.append(f_i_ZA["CO"])
        f_O2_ZA.append(f_i_ZA["O2"])
        f_CO2_ZA.append(f_i_ZA["CO2"])
        print("f_i_ZA",f_i_ZA)
        
 

    
    plt.figure(2)
    plt.plot(np.array(f_CO), '--k', label='CO')
    plt.plot(np.array(f_CO2), '-k', label=r'CO$_2$')
    
    plt.plot(np.array(f_CO_Prata), '--r', label='CO')
    plt.plot(np.array(f_CO2_Prata), '-r', label=r'CO$_2$')
    xvals = [6,2,6,2]
    labelLines(plt.gca().get_lines(), align=True, fontsize=10, ha='left',xvals=xvals)
    # plt.plot(np.array(f_CO_ZA), '--b', label='CO ZA')
    # plt.plot(np.array(f_CO2_ZA), '-b', label='CO2 ZA')

    plt.ylabel(r'$\dot{\omega}$, kg m$^{-2}$  s$^{-1}$',  rotation=0, horizontalalignment='left', labelpad=0)
    ax = plt.gca()
    ax.yaxis.set_label_coords(0., 1.01)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_yscale('log')
    xticks = np.arange(0, 8, 1)
    
    ax.set_xticks(xticks,labels=labels_O2, fontsize=MEDIUM_SIZE)
    # plt.legend(loc='best')
    # plt.xlabel('T [K]')
    # plt.legend(loc='best',frameon=False)
    plt.savefig('rates_oxidation.eps', format='eps',transparent=True,bbox_inches='tight', pad_inches=0.2)
    plt.savefig('rates_oxidation.pdf', format='pdf',transparent=True,bbox_inches='tight', pad_inches=0.2)
    plt.tight_layout()
    plt.show()    


def ZA_reactor():
    
    experimental_presure = np.array(
                    [2172])

    experimental_temperature = np.array(
                    [1518])
    
    labels =[   "A",
                "B",
                "C",
                "D",
                "E",
                "F",
                "G",
                "H"]

    labels_O2 = ['O2_'+i for i in labels]   
    
    f_CO = []
    f_O = []
    f_O2 = []
    f_CO2 = []
    
    for T,P in zip(experimental_temperature, experimental_presure):
    
        T_beam = T  # beam temperature [K]
        p_beam = P  # beam pressure [Pa]
        
        print("Temp", T, "Press",P)
        
        oxidation = OxidationRateModel()
        oxidation.pressure_in = p_beam
        oxidation.temperature_in = T_beam
        
        x_in = {"O": 0, "O2": 10}  # mole fraction of inflow gas: x_O or x_O2 [-]
        
        oxidation.composition_in = x_in

        oxidation_ZA= oxidation.oxidation_model("za_oxidation")
        
        oxidation_ZA.wall_temperature = T
        # oxidation_ZA.plot_surface_coverage()
        f_i = oxidation_ZA.solve_ODEs(T)
        print(f_i)
        f_CO.append(f_i["CO"])
        f_O2.append(f_i["O2"])
        f_CO2.append(f_i["CO2"])

# Press the green button in the gutter to run the script.
if __name__ == '__main__':
     #ZA_reactor()
     zeroD_reactor()

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
