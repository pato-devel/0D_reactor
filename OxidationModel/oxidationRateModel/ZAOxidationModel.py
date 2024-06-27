import OxidationModel.oxidationRateModel.OxidationRateModel as Ox
from OxidationModel.packages import ABC

class ZAOxidationModel(Ox.OxidationRateModelSelector, ABC):

    def __init__(self):
        print("Inside ZAOxidationModel class")
        return

