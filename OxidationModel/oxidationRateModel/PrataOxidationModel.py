"""
File: main.py
Authors: Bruno Dias and Jeremie Meurisse
Date: 06/27/24
Description: 0D reactor with different reaction rates using Prata model.
"""

import OxidationModel.oxidationRateModel.OxidationRateModel as Ox
from OxidationModel.packages import ABC
import matplotlib.pyplot as plt
import math
from scipy.integrate import odeint
from scipy.optimize import fsolve
import numpy as np
import sys, os

class PrataOxidationModel(Ox.OxidationRateModelSelector, ABC):
    """Oxidation class using Prata model

    Compute and plot the propability of products for O-atom
    Plot surface coverage for O-atom

    Attributes:
        None
    """

    # constant values
    T_beam = 1000 # beam temperature [K]
    p_beam = 2.4e-2 # beam pressure [Pa]
    B = 1e-5 # total active site density [mol/m2]
    Mm_O = 0.0159994 # molar mass of oxygen [kg/mol]
    Av =  6.022e23 # Avogadro number [1/mol]
    m_O = Mm_O / Av # mass of oxygen atom [kg]
    m_O2 = 2 * m_O # mass of oxygen [kg]
    k_B = 1.380649e-23 # Boltzmann constant [J/K]
    h = 6.62607015e-34 # Planck constant [J.s]
    R = 8.3144598 # Universal gas constant [J/mol/K]
    F_O = 1/4 * math.sqrt(8 * k_B * T_beam / (math.pi * m_O)) #  mean thermal speed of O-atom to the surface [m/s] 
    w_O = p_beam / (R * T_beam) # [O] oxygen concentration [mol/m3]
    f_Oin = F_O * w_O # flux of O-atom to the surface = F_O [O] [mol/m2/s]
    # F_O2 = 1/4 * math.sqrt(8 * k_B * T_beam / (math.pi * m_O2)) # flux of O2 to the surface = F_O2 [O2] [mol/m2/s]

    def __init__(self):
        print("Initialize PrataOxidationModel class")
        print("flux =",str(self.f_Oin*self.Av),"[atoms/m2/s]")
        self.plot_model_prediction()
        self.plot_surface_coverage()
        return

    def compute_rates(self, Tw):
        """
        compute the reaction rates

        :param Tw: surface temperature [K]
        :return: w_O oxygen concentration [mol/m3], k reaction rates [mol/m2/s]
        """ 
        # thermal dependent variables
        F_O_2D = math.sqrt(math.pi * self.k_B * Tw / (2 * self.m_O)) # mean thermal speed of the mobile adsorbed O-atom on the surface [m/s] 
        
        # reaction rates
        k={}
        k["kO1"] = self.F_O * 0.3 / self.B
        k["kO2"] = 2 * math.pi * self.m_O * self.k_B**2 * Tw**2 / (self.Av * self.B * self.h**3) * math.exp(-44277 /  Tw)
        k["kO3"] = self.F_O / self.B * 100 * math.exp(-4000 /  Tw)
        k["kO4"] = self.F_O / self.B * math.exp(-500 /  Tw)
        k["kO5"] = self.F_O / self.B * 0.7
        k["kO6"] = 2 * math.pi * self.m_O * self.k_B**2 *  Tw**2 / (self.Av * self.B * self.h**3) * math.exp(-96500 /  Tw)
        k["kO7"] = self.F_O / self.B * 1000 * math.exp(-4000 /  Tw)
        k["kO8"] = math.sqrt(self.Av / self.B) * F_O_2D * 1e-3 * math.exp(-15000/ Tw)
        k["kO9"] = math.sqrt(self.Av / self.B) * F_O_2D * 5e-5 * math.exp(-15000/ Tw)

        # k["kox1"] =  self.F_O2 / self.B**2 * math.exp(-8000 /  Tw)
        # k["kox2"] =  self.F_O2 / self.B * 100 * math.exp(-4000 /  Tw)
        # k["kox3"] =  self.F_O2 / self.B * math.exp(-500 /  Tw)
        # k["kox4"] =  self.F_O2 / self.B**2 * math.exp(-8000 /  Tw)
        # k["kox5"] =  self.F_O2 / self.B * 1000 * math.exp(-4000 /  Tw)

        return k

    def surface_coverage(self, k):
        """
        compute the steady-state surface concentrations

        :param k: reaction rates [mol/m2/s]
        :return: steady-state surface density of w_s empty sites [mol/m2], of w_Os absorbed oxygen with weakly bound [mol/m2], and of w_Oss absorbed oxygen with relatively strong bound [mol/m2]
        """ 
        A1 = 0 # 2 * k["kox1"] *  self.w_O2
        B1 = k["kO1"] *  self.w_O
        C1 = 2 * k["kO9"]
        D1 = k["kO2"]+ (k["kO3"] + k["kO4"]) *  self.w_O + 0
        A2 = 0
        B2 = k["kO5"] *  self.w_O
        C2 = 2 * k["kO8"]
        D2 = k["kO6"] + k["kO7"] *  self.w_O + 0
        A3 = 0
        B3 = 0
        C3 = 0
        D3 = 0

        def func_s(w_s):
            """
            compute the function for the steady-state surface density of empty sites

            :param w_s: steady-state surface density of empty sites [mol/m2]
            :return: f = 0, function for computing the steady-state surface density of empty sites
            """
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

    def fun_f_CO(self, Tw, t, k, w_s, w_Os, w_Oss):
        """
        compute the flux of CO product to the surface

        :param Tw: surface temperature [K]
        :return: f_CO = d[CO]/dt flux of CO product to the surface [mol/m2/s]
        """
        f_CO = k["kO3"] * self.w_O * w_Os + k["kO7"] * self.w_O * w_Oss
        return f_CO
    
    def fun_f_O(self, Tw, t, k, w_s, w_Os, w_Oss):
        """
        compute the flux of O product to the surface

        :param Tw: surface temperature [K]
        :return: f_O = d[O]/dt flux of O product to the surface [mol/m2/s]
        """
        press_O = self.p_beam # O-atom pressure [Pa] assume x_O = 1
        f_O = press_O/(self.Av * math.sqrt(2 * math.pi * self.m_O * self.k_B * self.T_beam)) - k["kO1"] * self.w_O * w_s + k["kO2"] * w_Os - k["kO4"] * self.w_O * w_Os - k["kO5"] * self.w_O * w_s + k["kO6"] * w_Oss + 0
        return f_O
    
    def fun_f_O2(self, Tw, t, k, w_s, w_Os, w_Oss):
        """
        compute the flux of O2 product to the surface

        :param Tw: surface temperature [K]
        :return: f_O2 = d[O2]/dt flux of O2 product to the surface [mol/m2/s]
        """
        f_O2 = 0 + k["kO8"] * w_Oss**2 + k["kO9"] * w_Os**2 + 0
        return f_O2
    
    def fun_f_CO2(self, Tw, t, k, w_s, w_Os, w_Oss):
        """
        compute the flux of CO2 product to the surface

        :param Tw: surface temperature [K]
        :return: f_CO2 = d[CO2]/dt flux of CO2 product to the surface [mol/m2/s]
        """
        f_CO2 = k["kO4"] * self.w_O * w_Os 
        return f_CO2
    
    def solve_ODEs(self, Tw):
        """
        solve ODEs for probability of products

        :param Tw: surface temperature [K]
        :return: [p_CO, p_O]: probability of products after one time step of 1 sec [-]
        """
        # compute reaction rates and surface coverage
        k = self.compute_rates(Tw)
        w_s, w_Os, w_Oss = self.surface_coverage(k)

        # flux of products
        t_span = (0,1)
        t_eval = np.linspace(t_span[0], t_span[1],2)
        sol_CO = odeint(self.fun_f_CO, Tw, t_eval, args=(k,w_s,w_Os,w_Oss))
        f_CO=sol_CO[1]-sol_CO[0]
        sol_O = odeint(self.fun_f_O, Tw, t_eval, args=(k,w_s,w_Os,w_Oss))
        f_O=sol_O[1]-sol_O[0]
        sol_O2 = odeint(self.fun_f_O2, Tw, t_eval, args=(k,w_s,w_Os,w_Oss))
        f_O2=sol_O2[1]-sol_O2[0]
        sol_CO2 = odeint(self.fun_f_CO2, Tw, t_eval, args=(k,w_s,w_Os,w_Oss))
        f_CO2=sol_CO2[1]-sol_CO2[0]

        # probabilities
        p_CO = f_CO / self.f_Oin
        p_O = f_O / self.f_Oin
        p_O2 = 2 * f_O2 / self.f_Oin # if x_O2 > 0: (2 * f_CO2 + f_CO) / (2 * self.f_O2in)
        p_CO2 = 2 * f_CO2 / self.f_Oin

        return p_CO, p_O, p_O2, p_CO2
    
    def plot_model_prediction(self):
        """
        plot the model prediction for O-atom
        """ 
        Tw = np.linspace(800,2400,100)
        p_CO = []
        p_O = []
        p_O2 = []
        p_CO2 = []
        for Tw_i in Tw:
            p_CO_i, p_O_i, p_O2_i, p_CO2_i  = self.solve_ODEs(Tw_i)
            p_CO.append(p_CO_i)
            p_O.append(p_O_i)
            p_O2.append(p_O2_i)
            p_CO2.append(p_CO2_i)

        p_CO_Prata = np.loadtxt("data/CO_probability.dat", delimiter=',')
        p_CO2_Prata = np.loadtxt("data/CO2_probability.dat", delimiter=',')
        p_O2_Prata = np.loadtxt("data/O2_probability.dat", delimiter=',')
        p_O_Prata = np.loadtxt("data/O_probability.dat", delimiter=',')

        plt.plot(p_CO_Prata[:,0], p_CO_Prata[:,1], 'k*', label='CO Prata')
        plt.plot(p_CO2_Prata[:,0], p_CO2_Prata[:,1], 'r*', label='CO2 Prata')
        plt.plot(p_O2_Prata[:,0], p_O2_Prata[:,1], 'g*', label='O2 Prata')
        plt.plot(p_O_Prata[:,0], p_O_Prata[:,1], 'b*', label='O Prata')


        plt.plot(Tw, np.array(p_CO), 'k', label='CO')
        plt.plot(Tw, np.array(p_CO2), 'r', label='CO2')
        plt.plot(Tw, np.array(p_O2), 'g', label='O2')
        plt.plot(Tw, np.array(p_O), 'b', label='O')
        plt.legend(loc='best')
        plt.xlabel('T [K]')
        plt.xlim(800,2400)
        plt.ylim(-0.1,1)
        plt.grid()
        home_dir=os.getenv("HOME")
        plt.savefig(home_dir+'/Desktop/fig1.png', transparent=True)
        plt.show()     

    def plot_surface_coverage(self):
        """
        plot the surface coverage for O-atom
        """ 
        Tw = np.linspace(800,2000,100)
        w_s = []
        w_Os = []
        w_Oss = []
        for Tw_i in Tw:
            k = self.compute_rates(Tw_i)
            w_s_i, w_Os_i, w_Oss_i = self.surface_coverage(k)
            w_s.append(w_s_i)
            w_Os.append(w_Os_i)
            w_Oss.append(w_Oss_i)

        print("coverage sum",(np.array(w_s) + np.array(w_Os) + np.array(w_Oss))/self.B )

        O_coverage_Prata =  np.loadtxt("data/O_coverage.dat",delimiter=',')
        O_star_coverage_Prata =  np.loadtxt("data/O_star_coverage.dat",delimiter=',')
        s_coverage_Prata =  np.loadtxt("data/s_coverage.dat",delimiter=',')
        total_O_coverage_Prata =  np.loadtxt("data/total_O_coverage.dat",delimiter=',')


        plt.plot(Tw, np.array(w_Os)/self.B, 'r', label='w_Os')
        plt.plot(O_coverage_Prata[:,0], O_coverage_Prata[:,1], 'r*', label='w_Os Prata')

        plt.plot(Tw, np.array(w_Oss)/self.B, 'b', label='w_Oss')
        plt.plot(O_star_coverage_Prata[:, 0], O_star_coverage_Prata[:, 1], 'b*', label='w_Oss Prata')

        plt.plot(Tw, np.array(w_s)/self.B, 'k', label='w_s')
        plt.plot(s_coverage_Prata[:, 0], s_coverage_Prata[:, 1], 'k*', label='w_s Prata')

        plt.plot(Tw, (np.array(w_Os)+np.array(w_Oss))/self.B, 'g', label='w_Os+w_Oss')
        plt.plot(total_O_coverage_Prata[:, 0], total_O_coverage_Prata[:, 1], 'g*', label='w_Os+w_Oss Prata')

        plt.legend(loc='best')
        plt.xlabel('T [K]')
        plt.grid()
        plt.yscale("log")
        plt.xlim(800,2000)
        plt.ylim(1e-4,1)
        home_dir=os.getenv("HOME")
        plt.savefig(home_dir+'/Desktop/fig2.png', transparent=True)
        plt.show()



