"""
Calculate Material Dispersion Materials outside of simple Dielectrics & Metals.

List of Materials:
 - Ta2O5

"""

import numpy as np

#%%
# Ta2O5
#------------------------------------------------------------------------------
def ta2o5(wavelength=900):

    """
    Calculates the refractive index of Ta2O5.
    
    Parameters
    ----------
    wavelegnth : numpy.array
        Wavelegnth of light [nm].
    
    Returns
    -------
    n : numpy.array
        Refractive index. 
    
    Notes
    -----
    Model referenced from (Bright 2013).

    """

    wavelength_um = wavelength * 1e-9 * 1e6  # wavelength (nm), convert to (m), then convert to (um)
    wavelength_cm = wavelength * 1e-9 * 1e2
    
    A = 2.06
    B = 0.025
    e_inf = (A + B/(wavelength_um)**2)**2
    
    w_p0 = 649              # in cm-1
    gamma_0 = 6.5 * 1e5
    w_cm = 1 / wavelength_cm
    t2 = w_p0**2 / (w_cm**2 - (1j * gamma_0 * w_cm))
    
    w_j = np.array([266, 500, 609, 672, 868, 3020])
    w_pj = np.array([1040, 573, 634, 408, 277, 373])
    gamma_j = np.array([188, 112, 88, 43, 113, 652])
    t3=[]
    for i in range(len(w_j)):
        t3.append(w_pj[i]**2 / (w_j[i]**2 - w_cm**2 - 1j * gamma_j[i] * w_cm))
    tt3 = t3[0] + t3[1] + t3[2] + t3[3] + t3[4] + t3[5]
    ε = e_inf - t2 + tt3
    n = (ε**.5).real

    return n

