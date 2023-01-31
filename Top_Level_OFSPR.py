"""
Top Level Script for modelling Optical Fibre Surface Plasmon Resonance (OFSPR).
"""

import numpy as np
import matplotlib.pyplot as plt
import DispDielectrics as dispD
import DispMetals as dispM
import DispOther as dispO

plt.close('all')


#%%
# Top Level Script:
#------------------------------------------------------------------------------
# Define wavelength range and step-size; all in (nm)
wavelength_array = np.linspace(450, 2000, 50)
plot=True

# Define all empty lists
n_BK7, n_schott, n_sio2M, n_sio2PY, e_Pd, e_Ag, n_ta2o5= ([] for i in range(7))

for i in range(len(wavelength_array)):
    n_BK7.append(dispD.BK7(wavelength=wavelength_array[i]))
    n_schott.append(dispD.schott(wavelength=wavelength_array[i]))
    n_sio2M.append(dispD.sio2M(wavelength=wavelength_array[i]))
    n_sio2PY.append(dispD.sio2PY(wavelength=wavelength_array[i]))
    e_Pd.append(dispM.BB_Pd(wavelength=wavelength_array[i]))
    e_Ag.append(dispM.BB_Ag(wavelength=wavelength_array[i]))
    n_ta2o5.append(dispO.ta2o5(wavelength=wavelength_array[i]))


#%%

if plot==True:
    # Plotting Dielectrics
    #------------------------------------------------------------------------------
    fig = plt.figure()
    ax = plt.subplot(111)
    ax.plot(wavelength_array, n_BK7)
    ax.plot(wavelength_array, n_schott)
    ax.plot(wavelength_array, n_sio2M)
    ax.plot(wavelength_array, n_sio2PY)
    ax.set_xlim([wavelength_array[0], wavelength_array[len(wavelength_array)-1]])
    plt.xlabel("Wavelength (nm)")
    plt.ylabel("Refractive Index")
    plt.title("Material Dispersion: Dielectrics")
    fig.legend((
            'BK7 Glass',
            'Schott Glass',
            'SiO2 Fibre Core',
            'SiO2 Thin Film',
            ), loc='upper right')
    plt.show()
    
    
    #%%
    # Plotting Metals
    fig = plt.figure(figsize=[12,8])
    ax = plt.subplot(111)
    ax.plot(wavelength_array, np.array(e_Ag).real, 'b')
    ax.plot(wavelength_array, np.array(e_Ag).imag, '--b')
    ax.plot(wavelength_array, np.array(e_Pd).real, 'r')
    ax.plot(wavelength_array, np.array(e_Pd).imag, '--r')
    ax.plot(wavelength_array, np.linspace(0, 0, len(wavelength_array)), '--k')
    ax.set_xlim([wavelength_array[0], wavelength_array[len(wavelength_array)-1]])
    plt.xlabel("Wavelength (nm)")
    plt.ylabel("Dielectric Permittivity")
    plt.title("Material Dispersion: Metals")
    fig.legend((
            'Ag: Real',
            'Ag: Imag',
            'Pd: Real',
            'Pd: Imag',
            'y = 0',
            ), loc='upper right')
    plt.show()
    
    #%%
    # Plotting Othe: Ta2O5
    fig = plt.figure(figsize=[12,8])
    ax = plt.subplot(111)
    ax.plot(wavelength_array, n_ta2o5)
    ax.set_xlim([wavelength_array[0], wavelength_array[len(wavelength_array)-1]])
    plt.xlabel("Wavelength (nm)")
    plt.ylabel("Dielectric Permittivity")
    plt.title("Material Dispersion: Ta2O5")
    plt.show()
    
    
    
