import OxidationModel.oxidationRateModel.OxidationRateModel as Ox
from OxidationModel.packages import ABC
import math
from scipy.integrate import odeint
from scipy.optimize import fsolve
import numpy as np
import sys
import mutationpp as mpp

class BrunoOxidationModel(Ox.OxidationRateModelSelector, ABC):
    
    speciesList = "O2 O CO CO2"
    mixtureOption = mpp.MixtureOptions()
    mixtureOption.setSpeciesDescriptor(speciesList)
    mixtureOption.setThermodynamicDatabase('RRHO')
    mixture = mpp.Mixture(mixtureOption)

    Mm_O = mixture.speciesMw(mixture.speciesIndex('O'))  # 0.0159994 # molar mass of oxygen [kg/mol]
    Mm_C = mixture.speciesMw(mixture.speciesIndex('C'))  # 0.012011 # molar mass of carbon [kg/mol]
    Mm_O2 = mixture.speciesMw(mixture.speciesIndex('O2'))  # 2 * Mm_O # molar mass of O2 [kg/mol]
    Mm_CO = mixture.speciesMw(mixture.speciesIndex('CO'))  # Mm_C + Mm_O # molar mass of CO [kg/mol]
    Mm_CO2 = mixture.speciesMw(mixture.speciesIndex('CO2'))  # Mm_C + 2 * Mm_O # molar mass of CO2 [kg/mol]
    Av = mixture.NA()  # Avogadro number [1/mol]
    m_O = Mm_O / Av  # mass of oxygen atom [kg]
    m_O2 = 2 * m_O  # mass of oxygen [kg]
    k_B = mixture.KB()  # Boltzmann constant [J/K]
    h = mixture.HP()  # Planck constant [J.s]
    R = mixture.RU()  # Universal gas constant [J/mol/K]
    
    def __init__(self, T_beam, p_beam,x_in):
        super().__init__()
        
        self.T_beam = T_beam
        self.p_beam = p_beam
        self.x_in = x_in
        self._wall_temperature = []

        print("Initialize brunoOxidation Model class")
        self.F_O = 1 / 4 * math.sqrt(8 * self.k_B * T_beam / (math.pi * self.m_O))  # mean thermal speed of O-atom to the surface [m/s]
        self.w_O = self.x_in["O"] * p_beam / (self.R * T_beam)  # [O] oxygen concentration [mol/m3]
        self.f_Oin = self.F_O * self.w_O  # flux of O-atom to the surface = F_O [O] [mol/m2/s]
        self.F_O2 = 1 / 4 * math.sqrt(
            8 * self.k_B * T_beam / (math.pi * self.m_O2))  # mean thermal speed of O-atom to the surface [m/s]
        self.w_O2 = self.x_in["O2"] * p_beam / (self.R * T_beam)  # [O2] oxygen concentration [mol/m3]
        self.f_O2in = self.F_O2 * self.w_O2  # flux of O2 to the surface = F_O2 [O2] [mol/m2/s]
        print("flux =", str(self.f_O2in * self.Av), "[atoms/m2/s]")
        return


    @Ox.OxidationRateModelSelector.wall_temperature.setter
    def wall_temperature(self, Twall):
        self._wall_temperature = Twall
    
    
    def compute_rates(self, Tw):
        """
        compute the reaction rates

        :param Tw: surface temperature [K]
        :k reaction rates [mol/m2/s]
        """

        # reaction rates
        k = {}

        # reaction 1 $C(s)+O_2(g) \rightarrow CO_2(g)$
        
        k["kr1"] =  0.4*np.exp(-1000/Tw)
        
        # reaction 2 $2 \ C(s)+  O_2(g) \rightarrow 2 \ CO(g)$
        k["kr2"] =0.3*np.power(Tw,0.5)*np.exp(-5105/Tw)
 
        return k    


    def fun_f_CO(self, Tw, t, k, w_s, w_Os, w_Oss):
        """
        compute the flux of CO product to the surface

        :param Tw: surface temperature [K]
        :return: f_CO = d[CO]/dt flux of CO product to the surface [mol/m2/s]
        """
        f_CO = 2*k["kr2"] * self.w_O2
        return f_CO
    
    def fun_f_O2(self, Tw, t, k, w_s, w_Os, w_Oss):
        """
        compute the flux of O2 product to the surface

        :param Tw: surface temperature [K]
        :return: f_O2 = d[O2]/dt flux of O2 product to the surface [mol/m2/s]
        """
        f_O2 = -k["kr1"] * self.w_O2 -k["kr2"] * self.w_O2

        return f_O2

    def fun_f_CO2(self, Tw, t, k, w_s, w_Os, w_Oss):
        """
        compute the flux of CO2 product to the surface

        :param Tw: surface temperature [K]
        :return: f_CO2 = d[CO2]/dt flux of CO2 product to the surface [mol/m2/s]
        """
        f_CO2 = k["kr1"] * self.w_O2

        return f_CO2
    
    def plot_reaction_rates(self):
        pass

    def plot_model_prediction(self):
        pass

    def plot_surface_coverage(self):
        pass
    
    def solve_ODEs(self, Tw):
        
        """
        solve ODEs for probability of products

        :param Tw: surface temperature [K]
        :return: [p_CO, p_O]: probability of products after one time step of 1 sec [-]
        """
        # compute reaction rates and surface coverage
        k = self.compute_rates(Tw)

        w_s = w_Os = w_Oss = 0.
        
        # flux of products
        t_span = (0, 100)
        t_eval = np.linspace(t_span[0], t_span[1], 1000)
        sol_CO = odeint(self.fun_f_CO, Tw, t_eval, args=(k, w_s, w_Os, w_Oss))
        f_CO = (sol_CO[-1] - sol_CO[0])/(t_span[1]-t_span[0])
        sol_O2 = odeint(self.fun_f_O2, Tw, t_eval, args=(k, w_s, w_Os, w_Oss))
        f_O2 = (sol_O2[-1] - sol_O2[0])/(t_span[1]-t_span[0])
        sol_CO2 = odeint(self.fun_f_CO2, Tw, t_eval, args=(k, w_s, w_Os, w_Oss))
        f_CO2 = (sol_CO2[-1] - sol_CO2[0])/(t_span[1]-t_span[0])

        f = {"CO": f_CO,"O2": f_O2, "CO2": f_CO2}

        return f
