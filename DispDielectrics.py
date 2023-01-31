"""
Calculates Material Dispersion of various Dielectric Materials.
 
List of Materials:   
 - BK7 Glass.
 - Schott Glass.
 - SiO2, thin film.
 - SiO2, fusde slica fibre core.

"""

import numpy as np


#%%
# BK7 Glass
#------------------------------------------------------------------------------
def BK7(wavelength=900):
    """
    Calculates the refractive index of BK7 glass.

    Parameters
    ----------
    wavelength : numpy.array
        Wavelengthof light [nm].

    Returns
    -------
    n : numpy.array
        Refractive index. 

    Notes
    -----
    Model referenced from (Refractiveinedx.info).

    """


    # wavelength (nm), convert to (m), then convert to (um)
    wavelength_um = wavelength * 1e-9 * 1e6
    a1 = 1.03961212
    a2 = 0.231792344
    a3 = 1.0146945
    b1 = 6.0006987 * 1e-3
    b2 = 2.00179144 * 1e-2
    b3 = 1.013560653 * 1e2
    n = np.sqrt(
            1 + 
            (a1 * wavelength_um**2)/(wavelength_um**2 - b1) +
            (a2 * wavelength_um**2)/(wavelength_um**2 - b2) +
            (a3 * wavelength_um**2)/(wavelength_um**2 - b3)
            )

    return(n)


#%%
# Schott Glass
#------------------------------------------------------------------------------
def schott(wavelength=900):

    """
    Calculates the refractive index of Schott glass.

    Parameters
    ----------
    wavelength : numpy.array
        Wavelength of light [nm].

    Returns
    -------
    n : numpy.array
        Refractive index. 

    Notes
    -----
    Model referenced from (Refractiveinedx.info).

    """
    # wavelength (nm), convert to (m), then convert to (um)
    wavelength_um = wavelength * 1e-9 * 1e6
    p1 = 1.4579
    p2 = -0.00365
    p3 = -0.00243
    p4 = 0.00012
    p5 = -4.606 *1e-6
    p6 = 9.635 * 1e-8
    n = (p1 + 
         p2 * wavelength_um**2 + 
         p3 * wavelength_um**-2 +
         p4 * wavelength_um**-2 + 
         p5 * wavelength_um**-6 + 
         p6 * wavelength_um **-8
         )

    return(n)


#%%
# SiO2, grown on a multilayer structure
#------------------------------------------------------------------------------
def sio2PY(wavelength=900):

    """
    Calculates the refractive index of thin film silica grown on a multilayer
    structure

    Parameters
    ----------
    wavelength : numpy.array
        Wavelength of light [nm].

    Returns
    -------
    n : numpy.array
        Refractive index. 

    Notes
    -----
    Model referenced from (Postava and Yamaguchi 2001).

    """

    # wavelength (nm), convert to (m), then convert to (um)
    wavelength_um = wavelength * 1e-9 * 1e6
    a1 = 1.1336
    b1 = 9.261*1e-2
    n = np.sqrt(1 + (a1 * wavelength_um**2) / (wavelength_um**2 - b1**2))

    return(n)


#%%
# Fused Silica
#------------------------------------------------------------------------------
def sio2M(wavelength=900):

    """
    Calculates the refractive index of fused silica

    Parameters
    ----------
    wavelength : numpy.array
        Wavelength of light [nm].

    Returns
    -------
    n : numpy.array
        Refractive index. 

    Notes
    -----
    Model referenced from (Malitson 1965).

    """

    # wavelength (nm), convert to (m), then convert to (um)
    wavelength_um = wavelength * 1e-9 * 1e6
    a1 = 0.6961663
    a2 = 0.407942
    a3 = 0.8974794
    b1 = 0.004679148
    b2 = 0.01351206
    b3 = 97.934000
    n = np.sqrt(
            1 + 
            (a1 * wavelength_um**2)/(wavelength_um**2 - b1) +
            (a2 * wavelength_um**2)/(wavelength_um**2 - b2) +
            (a3 * wavelength_um**2)/(wavelength_um**2 - b3)
            )

    return(n)

