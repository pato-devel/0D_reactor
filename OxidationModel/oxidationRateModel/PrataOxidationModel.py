import OxidationModel.oxidationRateModel.OxidationRateModel as Ox
from OxidationModel.packages import ABC

class PrataOxidationModel(Ox.OxidationRateModelSelector, ABC):

    def __init__(self):
        print("Inside PrataOxidationModel class")
        return

