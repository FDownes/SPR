"""
This is a script for testing dispersion for various materials.



"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from material_dispersion import n_Dielectrics as ndie
from material_dispersion import n_Metals as nmet
from material_dispersion import n_Other as noth


#%%

tst = ndie(wavelength_array=900)
print(tst)
print(tst.glass_bk7())
print(tst.glass_schott())
print(tst.sio2_py())
print(tst.sio2_m())

tst2 = nmet(wavelength_array=900)
print(tst2)
print(tst2.ag_bb())
print(tst2.pd_bb())

tst3 = noth(wavelength_array=900)
print(tst3)
print(tst3.ta2o5_ld())

#%%
wa = np.linspace(450, 2000, 50)
tst = ndie(wavelength_array=wa)
print(tst)
print(tst.glass_bk7())
print(tst.glass_schott())
print(tst.sio2_py())
print(tst.sio2_m())

tst2 = nmet(wavelength_array=wa)
print(tst2)
print(tst2.ag_bb())
print(tst2.pd_bb())

tst3 = noth(wavelength_array=wa)
print(tst3)
print(tst3.ta2o5_ld())
