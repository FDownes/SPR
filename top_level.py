"""
Top Level Script

@author: Fionn Downes
"""


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from concatenate_dispersion import CombineDispersion as cd
from concatenate_dispersion import PlotDispersion as plot_cd


#%%
# Define a wavelengt range in (nm):
wavelength_array = np.linspace(450, 2000, 50)

# Combine Dispersion into a DataFrame:
cd_mod = cd(wavelength_array=wavelength_array)
dis_df = cd_mod.get_df()

# Plot each figure:
plot_dis = plot_cd(dis_df)
plot_dis.dielectrics()
plot_dis.metals()
plot_dis.ta2o5()










