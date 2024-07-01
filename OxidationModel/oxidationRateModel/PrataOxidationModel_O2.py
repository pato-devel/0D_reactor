import OxidationModel.oxidationRateModel.OxidationRateModel as Ox
from OxidationModel.packages import ABC
import math
from scipy.integrate import odeint

class PrataOxidationModel(Ox.OxidationRateModelSelector, ABC):

    def __init__(self,Tw=1000):
        print("Inside PrataOxidationModel class")
        self.Tw = Tw # surface temperature [K]
        return

    def compute_rates(self):
        # constant variables
        T_beam = 1000 # beam temperature [K]
        B = 1e-5 # total active site density [mol/m2]
        Mm_O = 0.0159994 # molar mass of oxygen [kg/mol]
        Av =  6.022e23 # Avogadro number [1/mol]
        m_O = Mm_O / Av # mass of oxygen atom [kg]
        m_02 = m_O*2  # mass of oxygen molecular [kg]
        k_B = 1.380649e-23 # Boltzmann constant [J/K]
        h = 6.62607015e-34 # Planck constant [J.s]
        F_O = 1/4 * math.sqrt(8 * k_B * T_beam / (math.pi * m_O)) # flux of gas species to the surface = F_O [O] [mol/m2/s] 
        F_O2 = 1/4 * math.sqrt(8 * k_B * T_beam / (math.pi * m_O2)) # flux of gas species to the surface = F_O2 [O2] [mol/m2/s]
        F_O_2D = math.sqrt(math.pi * k_B * self.Tw / (2 * m_O)) # mean thermal speed of the mobile adsobed species on the surface [m/s]

        # reaction rates
        self.kO1 = F_O * 0.3 / B
        self.kO2 = 2 * math.pi * m_O * k_B**2 * self.Tw**2 / (Av * B * h**3) * math.exp(-44277 / self.Tw)
        self.kO3 = F_O / B * 100 * math.exp(-4000 / self.Tw)
        self.kO4 = F_O / B * math.exp(-500 / self.Tw)
        self.k05 = F_O / B * 0.7
        self.k06 = 2 * math.pi * m_O * k_B**2 * self.Tw**2 / (Av * B * h**3) * math.exp(-96500 / self.Tw)
        self.k07 = F_O / B * 1000 * math.exp(-4000 / self.Tw)
        self.k08 = math.sqrt(Av / B) * F_O_2D * 1e-3 * math.exp(-15000/self.Tw)
        self.k09 = math.sqrt(Av / B) * F_O_2D * 5e-5 * math.exp(-15000/self.Tw)
        self.k16 = F_O2/B**2*math.exp(-8000/self.Tw)
        self.k17 = F_O2 / B *100 * math.exp(-4000 / self.Tw)
        self.k18 = F_O2 / B  * math.exp(-500 / self.Tw)
        self.k19 = F_O2 / (B*B) * math.exp(-8000 / self.Tw)
        self.k20 = F_O2 / B * 1000 * math.exp(-4000 / self.Tw)
        
    # steady-state surface concentrations
    def surface_reaction_rates(self, w):
        w_O, w_Os, w_Oss, f_CO = w
        output={}
        output["A1"] = 0
        output["B1"] = self.kO1 * w_O
        output["C1"] = 2 * self.k09
        output["D1"] = self.kO2 + (self.kO3 + self.kO4) * w_O + 0
        output["A2"] = 0
        output["B2"] = self.kO5 * w_O
        output["C2"] = 2 * self.kO8
        output["D2"] = self.kO6




    
    def f_CO(self, w, t):
        w_O, w_Os, w_Oss, f_CO = w
        dwdt = [f_CO, self.kO3 * w_O * w_Os + self.kO7 * w_O * w_Oss]
        return dwdt
    
    def solve_ODEs(self):
        w0 = [1, 0, 0, 0]
        sol = odeint(f_CO, y0, t, args=(b, c))
        





