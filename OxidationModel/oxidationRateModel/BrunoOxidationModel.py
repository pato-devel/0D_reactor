import OxidationModel.oxidationRateModel.OxidationRateModel as Ox
from OxidationModel.packages import ABC

class BrunoOxidationModel(Ox.OxidationRateModelSelector, ABC):

    def __init__(self):
        print("Inside brunoOxidationModel class")
        return

