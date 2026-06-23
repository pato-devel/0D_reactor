import OxidationModel.oxidationRateModel.OxidationRateModel as Ox
from OxidationModel.packages import ABC
import math
from scipy.integrate import odeint
from scipy.optimize import fsolve
import numpy as np
import sys
import mutationpp as mpp

class ZAOxidationModel(Ox.OxidationRateModelSelector, ABC):

    speciesList = "O2 O CO CO2"
    mixtureOption = mpp.MixtureOptions()
    mixtureOption.setSpeciesDescriptor(speciesList)
    mixtureOption.setThermodynamicDatabase('RRHO')
    mixture = mpp.Mixture(mixtureOption)

    B = 3.5e19 # units [1/m2] from Alba, Z&A and Chen
    phi = 5.8E-5  # total active site density [mol/m2] page 63 Alba thesis
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

        print("Initialize ZAOxidationModel Model class")
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
    
    
    def compute_forward_constants(self, Tw):
        """
        compute the reaction rates

        :param Tw: surface temperature [K]
        :return: w_O oxygen concentration [mol/m3], k reaction rates [mol/m2/s]
        """
        # thermal dependent variables
        thermal_speed_o2 = math.sqrt(8*self.R*self.T_beam/(math.pi*self.Mm_O2))
        
        k = {}
        
        # 2O(s) ↔ O2 + 2(s) Arrh type page 61
        #k["k03"] = 3.58E13*Tw*math.exp(-256.07E3/(self.R*Tw))
        
        #  O2 + 2(s) ↔ 2O(s) from Z&A paper and Chen2012 reaction 2
        k["k03"] = 0.0008*self.B*(self.k_B*Tw/self.h)*np.exp(-30800./self.T_beam)
        
        # O2 +(s)↔O+O(s) ER type page 61
        k["k04"] = thermal_speed_o2/(4.*self.phi)*1.*math.exp(-118.06e3/(self.R*Tw))
        
        # O(s)+C(b)↔CO+(s)  Arrh type page 61
        k["k06"] = 2.08E9*Tw*math.exp(-332.56E3/(self.R*Tw))
        
        # 2O(s) + C(b) ↔ CO2 + 2(s) Arrh type page 61
        k["k08"] = 3.58e17*math.exp(-332.56E3/(self.R*Tw))
    
        return k
    
    def compute_rates(self, k, w_Os):
        """
        compute the reaction rates

        :param Tw: surface temperature [K], oxygen surface cover
        :return rates
        """
        
        w_es = self.phi - w_Os
        Po = 1.01325E5
        over_K_one = self.k_B*self.T_beam/Po*math.pow(2*math.pi*self.m_O*self.k_B*self.T_beam/self.h**2,1.5)*math.exp(-45000/self.T_beam)
        K_one = 1./over_K_one
        
        Gibbs_energy = self.mixture.getSTGibbsMass(self.T_beam)
        
        delta_gibbs_O2 = Gibbs_energy[self.mixture.speciesIndex('O2')] - 2*Gibbs_energy[self.mixture.speciesIndex('O')]
    
        Kc_O2 = np.exp(-delta_gibbs_O2/(self.R*self.T_beam))
        print(Gibbs_energy)
        print(delta_gibbs_O2)
        print("Kc_O2",delta_gibbs_O2/(self.R))
        
        Kp_O2 = np.power(Po/(self.R*self.T_beam),-1)*Kc_O2

        K3 = Kp_O2*K_one*K_one
        print(K3)
        kb3 = k["k03"]/K3       # backward rate 
        
        r= {}
        
        # 2O(s) ↔ O2 + 2(s) Arrh type page 61
       # r["r03"] = k["k03"]*w_Os**2 - kb3*self.w_O2*w_es**2
       
        #  O2 + 2(s) ↔ 2O(s) from Z&A paper and Chen2012 reaction 2
        r["r03"] = k["k03"]*(K3*self.p_beam*self.x_in["O2"]*(w_es/self.phi)**2 - (w_Os/self.phi)**2 ) 
        
        # O2 +(s)↔O+O(s) ER type page 61
        r["r04"] = k["k04"]*self.w_O2*w_es
        
        # O(s)+C(b)↔CO+(s)  Arrh type page 61
        r["r06"] = k["k06"]*w_Os
        
        # 2O(s) + C(b) ↔ CO2 + 2(s) Arrh type page 61
        r["r08"] = k["k08"]*w_Os**2
    
        return r
    
    def surface_coverage(self, k):
        """
        compute the steady-state surface concentrations

        :param k: reaction rates [mol/m2/s] :return: steady-state surface density of w_s empty sites [mol/m2],
        of w_Os absorbed oxygen with weakly bound [mol/m2], and of w_Oss absorbed oxygen with relatively strong bound
        [mol/m2]
        """

        def func_s(w_s):
            """
            compute the function for the steady-state surface density of empty sites

            :param w_s: steady-state surface density of empty sites [mol/m2]
            :return: f = 0, function for computing the steady-state surface density of empty sites
            """
            r = self.compute_rates(k, w_s)
            
            f = 2*r["r03"] + r["r04"] - r["r06"] -2*r["r08"]
    
            return f

        w_Os_root = fsolve(func_s, 0.)
        w_Os = w_Os_root[0]
        tol = 1e-12
        if abs(func_s(w_Os)) > tol:
            print("Error: f >", tol, "f =", str(func_s(w_Os)))
            sys.exit()
        w_s = self.phi - w_Os

        return w_s, w_Os
    

    def fun_f_CO(self, Tw, t, r, w_s, w_Os, w_Oss):
        """
        compute the flux of CO product to the surface

        :param Tw: surface temperature [K]
        :return: f_CO = d[CO]/dt flux of CO product to the surface [mol/m2/s]
        """
        f_CO = r["r06"]
        return f_CO
    
    def fun_f_O2(self, Tw, t, r, w_s, w_Os, w_Oss):
        """
        compute the flux of O2 product to the surface

        :param Tw: surface temperature [K]
        :return: f_O2 = d[O2]/dt flux of O2 product to the surface [mol/m2/s]
        """
        f_O2 = -r["r03"]-r["r04"]

        return f_O2

    def fun_f_CO2(self, Tw, t, r, w_s, w_Os, w_Oss):
        """
        compute the flux of CO2 product to the surface

        :param Tw: surface temperature [K]
        :return: f_CO2 = d[CO2]/dt flux of CO2 product to the surface [mol/m2/s]
        """
        f_CO2 = r["r08"]

        return f_CO2
    
    def fun_f_O(self, Tw, t, r, w_s, w_Os, w_Oss):
        """
        compute the flux of O product to the surface

        :param Tw: surface temperature [K]
        :return: f_O = d[O]/dt flux of O product to the surface [mol/m2/s]
        """
        f_O = r["r04"]
        return f_O
    
    def plot_reaction_rates(self):
        pass

    def plot_model_prediction(self):
        pass

    def plot_surface_coverage(self):
        
        Tw = np.array(self._wall_temperature)
        
        if Tw.size==0:
            raise TypeError("wall temperatures not set")
        
        w_s = []
        w_Os = []
        print("Temper",Tw)
        
        for Tw_i in Tw:
    
            k = self.compute_forward_constants(Tw)
            w_s_i, w_Os_i = self.surface_coverage(k)
            w_s.append(w_s_i)
            w_Os.append(w_Os_i)
        
        print("w_s", w_s, "w_Os", w_Os)
        
        pass
    
    def solve_ODEs(self, Tw):
        
        """
        solve ODEs for probability of products

        :param Tw: surface temperature [K]
        :return: [p_CO, p_O]: probability of products after one time step of 1 sec [-]
        """
        # compute reaction rates and surface coverage
    
        k = self.compute_forward_constants(Tw)
        w_s, w_Os = self.surface_coverage(k)
        w_Oss = 0.
        r = self.compute_rates(k,w_Os)
        
        K3 = self.w_O2*w_s**2/max(w_Os**2,1e-16)
        print(K3)
        kb3 = k["k03"]/K3
        print("kb3",kb3)
        print("k[O3]*w_Os**2",k["k03"]*w_Os**2)
        print("kb3*self.w_O2*w_s**2",kb3*self.w_O2*w_s**2)
        print(r)
        
        # flux of products
        t_span = (0, 1000)
        t_eval = np.linspace(t_span[0], t_span[1], 10000)
        sol_CO = odeint(self.fun_f_CO, Tw, t_eval, args=(r, w_s, w_Os, w_Oss))
        f_CO = (sol_CO[-1] - sol_CO[0])/(t_span[1]-t_span[0])
        sol_O = odeint(self.fun_f_O, Tw, t_eval, args=(r, w_s, w_Os, w_Oss))
        f_O = (sol_O[-1] - sol_O[0])/(t_span[1]-t_span[0])
        sol_O2 = odeint(self.fun_f_O2, Tw, t_eval, args=(r, w_s, w_Os, w_Oss))
        f_O2 = (sol_O2[-1] - sol_O2[0])/(t_span[1]-t_span[0])
        sol_CO2 = odeint(self.fun_f_CO2, Tw, t_eval, args=(r, w_s, w_Os, w_Oss))
        f_CO2 = (sol_CO2[-1] - sol_CO2[0])/(t_span[1]-t_span[0])

        f = {"CO": f_CO, "O": f_O, "O2": f_O2, "CO2": f_CO2}

        return f

