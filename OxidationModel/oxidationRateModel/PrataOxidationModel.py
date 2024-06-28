"""
File: main.py
Authors: Bruno Dias and Jeremie Meurisse
Date: 06/27/24
Description: 0D reactor with different reaction rates using Prata model.
"""

import OxidationModel.oxidationRateModel.OxidationRateModel as Ox
from OxidationModel.packages import ABC
import math
from scipy.integrate import odeint
from scipy.optimize import fsolve
import numpy as np
import sys

class PrataOxidationModel(Ox.OxidationRateModelSelector, ABC):

    # constant values
    T_beam = 1000 # beam temperature [K]
    p_beam = 2.4e-2 # beam pressure [Pa]
    B = 1e-5 # total active site density [mol/m2]
    Mm_O = 0.0159994 # molar mass of oxygen [kg/mol]
    Av =  6.022e23 # Avogadro number [1/mol]
    m_O = Mm_O / Av # mass of oxygen atom [kg]
    k_B = 1.380649e-23 # Boltzmann constant [J/K]
    h = 6.62607015e-34 # Planck constant [J.s]
    R = 8.3144598 # Universal gas constant [J/mol/K]
    
    def __init__(self):
        print("Initialize PrataOxidationModel class")
        self.solve_ODEs()
        return

    def compute_rates(self, Tw):
        """
        compute the reaction rates

        :param Tw: surface temperature [K]
        :return: w_O oxygen concentration [mol/m3], k reaction rates [mol/m2/s]
        """ 
        # constant variables
        F_O = 1/4 * math.sqrt(8 * self.k_B * self.T_beam / (math.pi * self.m_O)) # flux of gas species to the surface = F_O [O] [mol/m2/s] 
        F_O_2D = math.sqrt(math.pi * self.k_B *  Tw / (2 * self.m_O)) # mean thermal speed of the mobile adsorbed species on the surface [m/s] 
        w_O = self.p_beam / (self.R * self.T_beam) # oxygen concentration [mol/m3]

        # reaction rates
        k={}
        k["kO1"] = F_O * 0.3 / self.B
        k["kO2"] = 2 * math.pi * self.m_O * self.k_B**2 * Tw**2 / (self.Av * self.B * self.h**3) * math.exp(-44277 /  Tw)
        k["kO3"] = F_O / self.B * 100 * math.exp(-4000 /  Tw)
        k["kO4"] = F_O / self.B * math.exp(-500 /  Tw)
        k["kO5"] = F_O / self.B * 0.7
        k["kO6"] = 2 * math.pi * self.m_O * self.k_B**2 *  Tw**2 / (self.Av * self.B * self.h**3) * math.exp(-96500 /  Tw)
        k["kO7"] = F_O / self.B * 1000 * math.exp(-4000 /  Tw)
        k["kO8"] = math.sqrt(self.Av / self.B) * F_O_2D * 1e-3 * math.exp(-15000/ Tw)
        k["kO9"] = math.sqrt(self.Av / self.B) * F_O_2D * 5e-5 * math.exp(-15000/ Tw)

        return w_O, k

    def surface_reaction_rates(self, w_O, k):
        """
        compute the steady-state surface concentrations

        :param w_O: oxygen concentration [mol/m3]
        :param k: reaction rates [mol/m2/s]
        :return: steady-state surface density of w_s empty sites [mol/m2], w_Os absorbed oxygen with weakly bound [mol/m2], and w_Oss absorbed oxygen with relatively strong bound [mol/m2]
        """ 
        A1 = 0
        B1 = k["kO1"] * w_O
        C1 = 2 * k["kO9"]
        D1 = k["kO2"]+ (k["kO3"] + k["kO4"]) * w_O + 0
        A2 = 0
        B2 = k["kO5"] * w_O
        C2 = 2 * k["kO8"]
        D2 = k["kO6"] + k["kO7"] * w_O + 0
        A3 = 0
        B3 = 0
        C3 = 0
        D3 = 0

        def func_s(w_s):
            f = - w_s + self.B\
                   - 2 * (A1 * w_s**2 + B1 * w_s) / (D1 + math.sqrt(D1**2 + 4 * C1 * (A1 * w_s**2 + B1 * w_s)))\
                   - 2 * (A2 * w_s**2 + B2 * w_s) / (D2 + math.sqrt(D2**2 + 4 * C2 * (A2 * w_s**2 + B2 * w_s)))
                   #- 2 * (B3 * w_s) / (D3 + math.sqrt(D3**2 + 4 * C3 * B2 * w_s))
            return f
        
        w_s_root = fsolve(func_s, self.B)
        w_s = w_s_root[0]
        tol = 1e-12
        if abs(func_s(w_s)) > tol:
            print("Error: f > 0, f =",str(w_s_root[1]))
            sys.exit()
        w_Os = 2 * (A1 * w_s**2 + B1 * w_s) / (D1 + math.sqrt(D1**2 + 4 * C1 * (A1 * w_s**2 + B1 * w_s)))
        w_Oss = 2 * (A2 * w_s**2 + B2 * w_s) / (D2 + math.sqrt(D2**2 + 4 * C2 * (A2 * w_s**2 + B2 * w_s)))

        return w_s, w_Os, w_Oss
    
    # def f_CO(self, Tw, t):
    #     w_O, kO1, kO2, kO3, kO4, kO5, kO6, kO7, kO8, kO9 = self.compute_rates(Tw)
    #     w_Os, w_Oss = self.surface_reaction_rates(w_O)
    #     dwdt = kO3 * w_O * w_Os + kO7 * w_O * w_Oss
    #     return dwdt
    
    def solve_ODEs(self):
        # w0 = [1, 0] # Oxygen concentration [mol/m3]
        # t = np.linspace(0, 1, int(1e6))
        # sol = odeint(self.f_CO, w_O_0, t, args=())
        Tw = np.linspace(800,2000,100)
        w_s = []
        w_Os = []
        w_Oss = []
        for Tw_i in Tw:
            w_O, k = self.compute_rates(Tw_i)
            w_s_i, w_Os_i, w_Oss_i = self.surface_reaction_rates(w_O, k)
            w_s.append(w_s_i)
            w_Os.append(w_Os_i)
            w_Oss.append(w_Oss_i)
        import matplotlib.pyplot as plt

        plt.plot(Tw, np.array(w_Os)/self.B, 'r', label='w_Os')
        plt.plot(Tw, np.array(w_Oss)/self.B, 'b', label='w_Oss')
        plt.plot(Tw, np.array(w_s)/self.B, 'k', label='w_s')
        plt.plot(Tw, (np.array(w_Os)+np.array(w_Oss))/self.B, 'g', label='w_Os')

        # plt.ylim(1e-4,1)
        plt.legend(loc='best')
        plt.xlabel('t')
        plt.grid()
        plt.yscale("log")
        plt.show()

        print("")       





