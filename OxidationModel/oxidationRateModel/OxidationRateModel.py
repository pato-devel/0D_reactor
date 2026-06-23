import OxidationModel
from OxidationModel.packages import ABC,abstractmethod
import OxidationModel.oxidationRateModel as bOx

class OxidationRateModel:
    def __init__(self):

        self._Temperature_in = 0.
        self._pressure_in = 0.
        self._composition={}

        pass

    @property
    def temperature_in(self):
        return self._Temperature_in

    @temperature_in.setter
    def temperature_in(self, Temp):
        self._Temperature_in = Temp

    @property
    def pressure_in(self):
        return self._pressure_in

    @pressure_in.setter
    def pressure_in(self, Press):
        self._pressure_in = Press
        
    @property
    def composition_in(self):
        return self._composition_in

    @composition_in.setter
    def composition_in(self, comp):
        
        comp = {k: v / total for total in (sum(comp.values()),) for k, v in comp.items()}
        
        self._composition_in = comp    

    def oxidation_model(self, oxidation_model: str = "Bruno")-> object:
        return self.__oxidation_model_selector(oxidation_model)

    def __oxidation_model_selector(self,name: str):
        try:
            if name.lower() == "bruno_oxidation":
                return bOx.BrunoOxidationModel(self._Temperature_in,self.pressure_in, self.composition_in)
            elif name.lower() == "prata_oxidation":
                return bOx.PrataOxidationModel(self.temperature_in,self.pressure_in, self.composition_in)
            elif name.lower() == "za_oxidation":
                return bOx.ZAOxidationModel(self._Temperature_in,self.pressure_in, self.composition_in)
            raise TypeError("Oxidation type is not valid.")
        except TypeError as e:
            print(e)


class OxidationRateModelSelector(ABC):

    def __init__(self):
        self._wall_temperature = []
        return

    @abstractmethod
    def plot_reaction_rates(self):
        pass

    @abstractmethod
    def plot_model_prediction(self):
        pass

    @abstractmethod
    def plot_surface_coverage(self):
        pass
    
    @abstractmethod
    def solve_ODEs(self, Twall):
        pass
    
    
    @property
    def wall_temperature(self):
        return self._wall_temperature
    
    @wall_temperature.setter
    @abstractmethod
    def wall_temperature(self, Twall):
        pass