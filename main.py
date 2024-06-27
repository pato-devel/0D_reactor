# This is a sample Python script.
import mutationpp
# Press ⌃R to execute it or replace it with your code.
# Press Double ⇧ to search everywhere for classes, files, tool windows, actions, and settings.
import mutationpp as mpp

import OxidationModel
from OxidationModel import OxidationRateModel
from OxidationModel import packages as pk


def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
    print(f'Hi, {name}')  # Press ⌘F8 to toggle the breakpoint.
    # help(mpp.Mixture)
    oxidation = OxidationRateModel("Bruno_oxidation")
    oxidation_prata = OxidationRateModel("prata_oxidation")
    oxidation_za = OxidationRateModel("ZA_oxidation")
    variable = pk.np.linspace(1,100,10)
    print(variable)

# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    print_hi('PyCharm')

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
