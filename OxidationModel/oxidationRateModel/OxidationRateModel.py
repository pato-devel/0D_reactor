import OxidationModel
from OxidationModel.packages import ABC
import OxidationModel.oxidationRateModel as bOx

class OxidationRateModel:
    def __init__(self, oxidation_model: str = "Bruno") -> object:
        self.__oxidation_model = self.__oxidation_model_selector(oxidation_model)

    @property
    def oxidation_model(self):
        return self.__oxidation_model

    @staticmethod
    def __oxidation_model_selector(name: str):
        try:
            if name.lower() == "bruno_oxidation":
                return bOx.BrunoOxidationModel()
            elif name.lower() == "prata_oxidation":
                return bOx.PrataOxidationModel()
            elif name.lower() == "za_oxidation":
                return bOx.ZAOxidationModel()
            raise TypeError("Oxidation type is not valid.")
        except TypeError as e:
            print(e)


class OxidationRateModelSelector(ABC):

    def __init__(self):
        print("Inside master class")
        return
