"""
Retrive and Combine Material Dispersion Information

Classes
-------
 - CombineDispersion()
 - PlotDispersion()

Makes Use Of
------------
 - n_Dielectrics()
 - n_Metals()
 - n_Other()

"""

__all__ = ['CombineDispersion', 'PlotDispersion',]

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from material_dispersion import n_Dielectrics as ndie
from material_dispersion import n_Metals as nmet
from material_dispersion import n_Other as noth


class CombineDispersion():
    """
    Generate the refractive index of various dielectric materials

    Attributes
    ----------
    wavelength_array : numpy.array
        Wavelength / colour of light [nm].

    Methods
    -------
    - get_df
    """

    def __init__(self, wavelength_array,):
        """Initialize; see class doc."""
        self.wavelength_array = wavelength_array


    def get_df(self,):
        """
        Generates the material dispersion for several materials. 
        
        Parameters
        ----------
        wavelength_array : numpy_array
            Wavelength of incident light (nm).
    
        Returns
        -------
        pd : Pandas DataFramedispersion_list.
            Dataframe containing wavelength (nm), and material dispersion from each
            material in material_list.
        """
    
        # Convery wavelength_array to numpy array if it is not the case already:
        wavelength_array = np.array(self.wavelength_array)
    
        # Call each class in material_dispersion:
        materials_die = ndie(wavelength_array=wavelength_array)
        materials_met = nmet(wavelength_array=wavelength_array)
        materials_oth = noth(wavelength_array=wavelength_array)
    
        # Define empty lists for material dispersion:
        n_glass_bk7, n_glass_schott, n_sio2_m, n_sio2_py, = ([] for i in range(4))
        n_pd, n_ag, = ([] for i in range(2))
        n_ta2o5 = []
    
        # Dielectrics:
        n_glass_bk7 = materials_die.glass_bk7()
        n_glass_schott = materials_die.glass_schott()
        n_sio2_m = materials_die.sio2_py()
        n_sio2_py = materials_die.sio2_m()
        # Metals:
        n_ag = np.sqrt(materials_met.ag_bb())
        n_pd = np.sqrt(materials_met.pd_bb())
        # Other:
        n_ta2o5 = materials_oth.ta2o5_ld()
    
        # Add all materials to a Pandas dataframe:
        if wavelength_array.size==1:
            # Avoiding ValueError: If using all scalar values, you must pass an index
            self.df = pd.DataFrame({
                'wavelength' : [wavelength_array],
                'glass_bk7' :  [n_glass_bk7],
                'glass_schott' : [n_glass_schott],
                'sio2_m' : [n_sio2_m],
                'sio2_py' : [n_sio2_py],
                'pd' : [n_pd],
                'ag' : [n_ag],
                'ta2o5' : [n_ta2o5],
                })
    
        else:
            self.df = pd.DataFrame({
                'wavelength' : wavelength_array,
                'glass_bk7' :  n_glass_bk7,
                'glass_schott' : n_glass_schott,
                'sio2_m' : n_sio2_m,
                'sio2_py' : n_sio2_py,
                'pd' : n_pd,
                'ag' : n_ag,
                'ta2o5' : n_ta2o5,
                })

        return self.df


class PlotDispersion():
    """
    Plot the Material Dispersion of components in 

    Attributes
    ----------
    df : Pandas Dataaframe
        Wavelength / colour of light [nm].

    Methods
    -------
    - 
    """

    def __init__(self, df,):
        """Initialize; see class doc."""
        self.df = df


    def dielectrics(self):
        """Plotting refractive index of Dielectrics."""
        # Set up figure:
        fig = plt.figure()
        ax = plt.subplot(111)
        ax.set_xlim([self.df.wavelength.iloc[0], self.df.wavelength.iloc[-1]])
        # Plot all refractive indices from the df:
        ax.plot(self.df.wavelength, self.df.glass_bk7)
        ax.plot(self.df.wavelength, self.df.glass_schott)
        ax.plot(self.df.wavelength, self.df.sio2_m)
        ax.plot(self.df.wavelength, self.df.sio2_py)
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


    def metals(self):
        """Plotting dielectric permittivity of Metals."""
        # Set up figure:
        fig = plt.figure()
        ax = plt.subplot(111)
        ax.set_xlim([self.df.wavelength.iloc[0], self.df.wavelength.iloc[-1]])
        # Plot the real and imagffinart component of all dielectric permittivities:
        ax.plot(self.df.wavelength, (np.array(self.df.ag)**2).real, 'b-')
        ax.plot(self.df.wavelength, (np.array(self.df.ag)**2).imag, 'b--')
        ax.plot(self.df.wavelength, (np.array(self.df.pd)**2).real, 'r-')
        ax.plot(self.df.wavelength, (np.array(self.df.pd)**2).imag, 'r--')
        # Include a zero crossing:
        ax.plot(self.df.wavelength, np.linspace(0, 0, len(self.df.wavelength)), '--k')
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


    def ta2o5(self):
        """Plotting refractive index of Ta2O5."""
        # Set up figure:
        fig = plt.figure()
        ax = plt.subplot(111)
        ax.set_xlim([self.df.wavelength.iloc[0], self.df.wavelength.iloc[-1]])
        # Plot all refractive indices from the df:
        ax.plot(self.df.wavelength, self.df.ta2o5)
        # Assign labels:
        plt.xlabel("Wavelength (nm)")
        plt.ylabel("Refractive Index")
        plt.title("Material Dispersion: $\mathregular{Ta_2O_5}$")
        plt.show()
