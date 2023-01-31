# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 10:00:04 2022

@author: fionn
"""

def plot_dielectrics(df):
    """Plotting refractive index of Dielectrics."""
    # Set up figure:
    fig = plt.figure()
    ax = plt.subplot(111)
    ax.set_xlim([df.wavelength.iloc[0], df.wavelength.iloc[-1]])
    # Plot all refractive indices from the df:
    ax.plot(df.wavelength, df.glass_bk7)
    ax.plot(df.wavelength, df.glass_schott)
    ax.plot(df.wavelength, df.sio2_m)
    ax.plot(df.wavelength, df.sio2_py)
    # Assign labels:
    plt.xlabel("Wavelength (nm)")
    plt.ylabel("Refractive Index")
    plt.title("Material Dispersion: Dielectrics")
    fig.legend((
            'BK7 Glass',
            'Schott Glass',
            '$\mathregular{SiO_2}$ Fibre Core',
            '$\mathregular{SiO_2}$ Thin Film',
            ), loc='upper right')
    plt.show()


#%%

def plot_metals(df):
    """Plotting dielectric permittivity of Metals."""
    # Set up figure:
    fig = plt.figure()
    ax = plt.subplot(111)
    ax.set_xlim([df.wavelength.iloc[0], df.wavelength.iloc[-1]])
    # Plot the real and imagffinart component of all dielectric permittivities:
    ax.plot(df.wavelength, (np.array(df.ag)**2).real, 'b-')
    ax.plot(df.wavelength, (np.array(df.ag)**2).imag, 'b--')
    ax.plot(df.wavelength, (np.array(df.pd)**2).real, 'r-')
    ax.plot(df.wavelength, (np.array(df.pd)**2).imag, 'r--')
    # Include a zero crossing:
    ax.plot(df.wavelength, np.linspace(0, 0, len(df.wavelength)), '--k')
    # Assign labels:
    plt.xlabel("Wavelength (nm)")
    plt.ylabel("Dielectric Permittivity")
    plt.title("Material Dispersion: Metals")
    fig.legend((
            'Ag: real',
            'Ag: imag',
            'Pd: real',
            'Pd: imag',
            ), loc='upper right')
    plt.show()


#%%

def plot_ta2o5(df):
    """Plotting refractive index of Ta2O5."""
    # Set up figure:
    fig = plt.figure()
    ax = plt.subplot(111)
    ax.set_xlim([df.wavelength.iloc[0], df.wavelength.iloc[-1]])
    # Plot all refractive indices from the df:
    ax.plot(df.wavelength, df.ta2o5)
    # Assign labels:
    plt.xlabel("Wavelength (nm)")
    plt.ylabel("Refractive Index")
    plt.title("Material Dispersion: $\mathregular{Ta_2O_5}$")
    plt.show()
