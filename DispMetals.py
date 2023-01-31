"""
Calculate Material Dispersion for Metals using the Brendel Bormann Model
Original data: Rakić et al. 1998, https://doi.org/10.1364/AO.37.005271

List of Materials:
 - Pd
 - Ag

"""

import numpy as np
from scipy.special import wofz as w
π = np.pi

#%%
# Pd
#------------------------------------------------------------------------------
def BB_Pd(wavelength=900):

    """
    Calculates the refractive index of Pd.

    Parameters
    ----------
    wavelegnth : numpy.array
        Wavelegnth of light [nm].

    Returns
    -------
    ε : numpy.array
        Dielectric permittivity (complex representation) 

    Notes
    -----
    Model referenced from:
        - Rakić 1998.
        - Brendel and Bormann 1992.

    """

    wavelength_m = wavelength * 1e-9       
    h = 4.135667516 * 1e-15     # (eV*s)
    c = 299792458               # (m/s)
    ω = h * c / wavelength_m             # (eV)

    ωp = 9.72  #eV
    f0 = 0.330
    Γ0 = 0.009 #eV
    
    f1 = 0.769
    Γ1 = 2.343 #eV
    ω1 = 0.066 #eV
    σ1 = 0.694 #eV
    
    f2 = 0.093
    Γ2 = 0.497 #eV
    ω2 = 0.502 #eV
    σ2 = 0.027 #eV
    
    f3 = 0.309
    Γ3 = 2.022 #eV
    ω3 = 2.432 #eV
    σ3 = 1.167 #eV
    
    f4 = 0.409
    Γ4 = 0.119 #eV
    ω4 = 5.987 #eV
    σ4 = 1.331 #eV
    
    Ωp = f0**.5 * ωp  #eV
    ε = 1-Ωp**2/(ω*(ω+1j*Γ0))

    α = (ω**2+1j*ω*Γ1)**.5
    za = (α-ω1)/(2**.5*σ1)
    zb = (α+ω1)/(2**.5*σ1)
    ε += 1j*π**.5*f1*ωp**2 / (2**1.5*α*σ1) * (w(za)+w(zb)) #χ1
    
    α = (ω**2+1j*ω*Γ2)**.5
    za = (α-ω2)/(2**.5*σ2)
    zb = (α+ω2)/(2**.5*σ2)
    ε += 1j*π**.5*f2*ωp**2 / (2**1.5*α*σ2) * (w(za)+w(zb)) #χ2
    
    α = (ω**2+1j*ω*Γ3)**.5
    za = (α-ω3)/(2**.5*σ3)
    zb = (α+ω3)/(2**.5*σ3)
    ε += 1j*π**.5*f3*ωp**2 / (2**1.5*α*σ3) * (w(za)+w(zb)) #χ3
    
    α = (ω**2+1j*ω*Γ4)**.5
    za = (α-ω4)/(2**.5*σ4)
    zb = (α+ω4)/(2**.5*σ4)
    ε += 1j*π**.5*f4*ωp**2 / (2**1.5*α*σ4) * (w(za)+w(zb)) #χ4
    
    return ε


#%%
# Ag
#------------------------------------------------------------------------------
def BB_Ag(wavelength=900):

    """
    Calculates the refractive index of Ag.

    Parameters
    ----------
    wavelegnth : numpy.array
        Wavelegnth of light [nm].

    Returns
    -------
    ε : numpy.array
        Dielectric permittivity (complex representation) 

    Notes
    -----
    Model referenced from:
        - Rakić 1998.
        - Brendel and Bormann 1992.

    """

    wavelength_m = wavelength * 1e-9       
    h = 4.135667516 * 1e-15     # (eV*s)
    c = 299792458               # (m/s)
    ω = h * c / wavelength_m             # (eV)

    ωp = 9.01  #eV
    f0 = 0.821
    Γ0 = 0.049 #eV
    
    f1 = 0.050
    Γ1 = 0.189 #eV
    ω1 = 2.025 #eV
    σ1 = 1.894 #eV
    
    f2 = 0.133
    Γ2 = 0.067 #eV
    ω2 = 5.185 #eV
    σ2 = 0.665 #eV
    
    f3 = 0.051
    Γ3 = 0.019 #eV
    ω3 = 4.343 #eV
    σ3 = 0.189 #eV
    
    f4 = 0.467
    Γ4 = 0.117 #eV
    ω4 = 9.809 #eV
    σ4 = 1.170 #eV
    
    f5 = 4.000
    Γ5 = 0.052 #eV
    ω5 = 18.56 #eV
    σ5 = 0.516 #eV
    
    Ωp = f0**.5 * ωp  #eV
    ε = 1-Ωp**2/(ω*(ω+1j*Γ0))

    α = (ω**2+1j*ω*Γ1)**.5
    za = (α-ω1)/(2**.5*σ1)
    zb = (α+ω1)/(2**.5*σ1)
    ε += 1j*π**.5*f1*ωp**2 / (2**1.5*α*σ1) * (w(za)+w(zb)) #χ1
    
    α = (ω**2+1j*ω*Γ2)**.5
    za = (α-ω2)/(2**.5*σ2)
    zb = (α+ω2)/(2**.5*σ2)
    ε += 1j*π**.5*f2*ωp**2 / (2**1.5*α*σ2) * (w(za)+w(zb)) #χ2
    
    α = (ω**2+1j*ω*Γ3)**.5
    za = (α-ω3)/(2**.5*σ3)
    zb = (α+ω3)/(2**.5*σ3)
    ε += 1j*π**.5*f3*ωp**2 / (2**1.5*α*σ3) * (w(za)+w(zb)) #χ3
    
    α = (ω**2+1j*ω*Γ4)**.5
    za = (α-ω4)/(2**.5*σ4)
    zb = (α+ω4)/(2**.5*σ4)
    ε += 1j*π**.5*f4*ωp**2 / (2**1.5*α*σ4) * (w(za)+w(zb)) #χ4
    
    α = (ω**2+1j*ω*Γ5)**.5
    za = (α-ω5)/(2**.5*σ5)
    zb = (α+ω5)/(2**.5*σ5)
    ε += 1j*π**.5*f5*ωp**2 / (2**1.5*α*σ5) * (w(za)+w(zb)) #χ5 
    
    return ε

