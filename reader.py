import pandas as pd
import numpy as np
# import SPMSys

file_name = "PSCompPars_2024.08.05_04.42.23.csv"

# Čtení ze souboru:
# (Pandas dataframe)
df = pd.read_csv(file_name, skiprows=54)
names = df['pl_name'].values
radii  = df['pl_rade'].values      # poloměry planet (v poloměrech Země)
masses = df['pl_bmasse'].values    # hmotnosti planet (v hmotnostech Země)
eorbs  = df['pl_orbeccen'].values  # výstřednosti dráhy
aorbs  = df['pl_orbsmax'].values   # hlavní poloosy
mstars = df['st_mass'].values      # hmotnosti mateřských hvězd

# class Exoplanet():
#     def __init__(self, name):
#         self.name = name
#         self.regime = regime

#         self.RPlanet = radii[names.index(self.name)] * SPMSys.REarth
#         self.MPlanet = masses[names.index(self.name)] * SPMSys.MEarth
#         self.eccPlanet = eorbs[names.index(self.name)]
#         self.aPlanet = aorbs[names.index(self.name)]

#         self.MStar = mstars[names.index(self.name)] * SPMSys.MSun

# Hodnoty uložené v dataframu lze převést na numpy pole pomocí ".values":
#již numpy.ndarray
# print(radii.values)
# print(names)
# print(min(eorbs), max(eorbs))
# print(min(masses), max(masses))