"""
Calculates Material Dispersion of various Dielectric and Metallic Materials.

@author: Fionn Downes

Classes
-------
 - n_Dielectrics
 - n_Metals
 - n_Other


 List of Materials
 -----------------
 - BK7 Glass.
 - Schott Glass.
 - SiO2, thin film.
 - SiO2, fused slica fibre core.
 - Ag, Brendel Bormann Model.
 - Pd, Brendel Bormann Model.
 - Ta2O5, Lorentz-Drude Model.
"""

__all__ = ['n_Dielectrics', 'n_Metals', 'n_Other',]


import numpy as np
from scipy.special import wofz


class n_Dielectrics():
    """
    Generate the refractive index of various dielectric materials

    Attributes
    ----------
    wavelength_array : numpy.array
        Wavelength / colour of light [nm].

    Methods
    -------
    - glass_bk7
    - glass_schott    
    - sio2_py
    - sio2_m
    """

    def __init__(self, wavelength_array,):
        """Initialize; see class doc."""
        self.wavelength_array = wavelength_array


    def glass_bk7(self):
        """
        Calculates the refractive index of BK7 glass.
    
        Parameters
        ----------
        wavelength_array : numpy.array
            Wavelengthof light [nm].
    
        Returns
        -------
        n : numpy.array
            Refractive index. 
    
        Info
        ----
        - Schott Zemax catalog (2017-01-20b) http://www.schott.com/
        - https://refractiveindex.info/
        """

        # wavelength (nm), convert to (m), then convert to (um)
        self.wavelength_um = self.wavelength_array * 1e-9 * 1e6
        self.a1 = 1.03961212
        self.a2 = 0.231792344
        self.a3 = 1.0146945
        self.b1 = 6.0006987 * 1e-3
        self.b2 = 2.00179144 * 1e-2
        self.b3 = 1.013560653 * 1e2
        self.n = np.sqrt(
                1 + 
                (self.a1 * self.wavelength_um**2)/(self.wavelength_um**2 - self.b1) +
                (self.a2 * self.wavelength_um**2)/(self.wavelength_um**2 - self.b2) +
                (self.a3 * self.wavelength_um**2)/(self.wavelength_um**2 - self.b3)
                )
    
        return self.n


    def glass_schott(self):
        """
        Calculates the refractive index of Schott glass.
    
        Parameters
        ----------
        wavelength_array : numpy.array
            Wavelength of light [nm].
    
        Returns
        -------
        n : numpy.array
            Refractive index. 
    
        Info
        ----
        - Schott Zemax catalog (2017-01-20b) http://www.schott.com/.
        - https://refractiveindex.info/.
        """

        # wavelength (nm), convert to (m), then convert to (um)
        self.wavelength_um = self.wavelength_array * 1e-9 * 1e6
        self.p1 = 1.4579
        self.p2 = -0.00365
        self.p3 = -0.00243
        self.p4 = 0.00012
        self.p5 = -4.606 *1e-6
        self.p6 = 9.635 * 1e-8
        self.n = (self.p1 + 
             self.p2 * self.wavelength_um**2 + 
             self.p3 * self.wavelength_um**-2 +
             self.p4 * self.wavelength_um**-2 + 
             self.p5 * self.wavelength_um**-6 + 
             self.p6 * self.wavelength_um **-8
             )
    
        return self.n


    def sio2_py(self):
        """
        Calculates the refractive index of thin film silica grown on a multilayer
        structure.

        Parameters
        ----------
        wavelength_array  : numpy.array
            Wavelength of light [nm].

        Returns
        -------
        n : numpy.array
            Refractive index. 

        Info
        ----
        - Postava and Yamaguchi (2001) Optical functions of low-k materials for
          interlayer dielectrics.
        """

        # wavelength (nm), convert to (m), then convert to (um)
        self.wavelength_um = self.wavelength_array * 1e-9 * 1e6
        self.a1 = 1.1336
        self.b1 = 9.261 * 1e-2
        self.n = np.sqrt(1 + (self.a1 * self.wavelength_um**2) / (self.wavelength_um**2 - self.b1**2))
    
        return self.n


    def sio2_m(self):
        """
        Calculates the refractive index of fused silica.

        Parameters
        ----------
        wavelength_array  : numpy.array
            Wavelength of light [nm].

        Returns
        -------
        n : numpy.array
            Refractive index. 

        Info
        ----
        - Malitson (1965) Interspecimen Comparison of the Refractive Index of
          Fused Silica.
        """

        # wavelength (nm), convert to (m), then convert to (um)
        self.wavelength_um = self.wavelength_array * 1e-9 * 1e6
        self.a1 = 0.6961663
        self.a2 = 0.407942
        self.a3 = 0.8974794
        self.b1 = 0.004679148
        self.b2 = 0.01351206
        self.b3 = 97.934000
        self.n = np.sqrt(
                1 + 
                (self.a1 * self.wavelength_um**2) / (self.wavelength_um**2 - self.b1) +
                (self.a2 * self.wavelength_um**2) / (self.wavelength_um**2 - self.b2) +
                (self.a3 * self.wavelength_um**2) / (self.wavelength_um**2 - self.b3)
                )

        return self.n


    def __repr__(self):
        return repr(f'Refractive index of Dielectrics at wavelength: ' + str(self.wavelength_array) + ' nm')



class n_Metals():
    """
    Generate the refractive index of various metallic materials.

    Attributes
    ----------
    wavelength_array : numpy.array
        Wavelength / colour of light [nm].

    Methods
    -------
    - ag_bb
    - pd_bb
    """

    def __init__(self, wavelength_array,):
        """Initialize; see class doc."""
        self.wavelength_array = wavelength_array


    def ag_bb(self):
    
        """
        Calculates the refractive index of Ag.
    
        Parameters
        ----------
        wavelegnth : numpy.array
            Wavelegnth of light [nm].
    
        Returns
        -------
        eps : numpy.array
            Dielectric permittivity (complex representation).
    
        Info
        ----
        - Brendel and Bormann (1992) An infrared dielectric function model for
          amorphous solids.
        - Rakic (1998) Optical properties of metallic films for
          vertical-cavity optoelectronic devices.

        Parameters
        ----------
        - h (eVs),
        - c (m/s),
        - w (eV),
        - wp (eV),
        - gamma0 (eV),
        - gamma1 (eV),
        - gamma2 (eV),
        - gamma3 (eV),
        - gamma4 (eV),
        - sigma0 (eV),
        - sigma1 (eV),
        - sigma2 (eV),
        - sigma3 (eV),
        - sigma4 (eV),
        - omegap (eV),
        - w is the Faddeeva function,
        """

        # Constants:
        self.wavelength_m = self.wavelength_array * 1e-9       
        self.h = 4.135667516 * 1e-15
        self.c = 299792458
        self.w = self.h * self.c / self.wavelength_m
        self.wp = 9.01
        self.f0 = 0.821
        self.gamma0 = 0.049
        self.f1 = 0.050
        self.gamma1 = 0.189
        self.w1 = 2.025
        self.sigma1 = 1.894
        self.f2 = 0.133
        self.gamma2 = 0.067
        self.w2 = 5.185
        self.sigma2 = 0.665
        self.f3 = 0.051
        self.gamma3 = 0.019
        self.w3 = 4.343
        self.sigma3 = 0.189
        self.f4 = 0.467
        self.gamma4 = 0.117
        self.w4 = 9.809
        self.sigma4 = 1.170
        self.f5 = 4.000
        self.gamma5 = 0.052
        self.w5 = 18.56
        self.sigma5 = 0.516

        # Equations:
        self.omegap = self.f0**.5 * self.wp
        self.eps = 1 - self.omegap**2 / (self.w * (self.w + 1j * self.gamma0))
        self.alpha = (self.w**2 + 1j * self.w * self.gamma1)**.5
        self.za = (self.alpha - self.w1) / (2**.5 * self.sigma1)
        self.zb = (self.alpha + self.w1) / (2**.5 * self.sigma1)
        self.eps += 1j * (np.pi)**.5 * self.f1 * self.wp**2 / (2**1.5 * self.alpha * self.sigma1) * (wofz(self.za) + wofz(self.zb))
        self.alpha = (self.w**2 + 1j * self.w * self.gamma2)**.5
        self.za = (self.alpha - self.w2) / (2**.5 * self.sigma2)
        self.zb = (self.alpha + self.w2) / (2**.5 * self.sigma2)
        self.eps += 1j * (np.pi)**.5 * self.f2 * self.wp**2 / (2**1.5 * self.alpha * self.sigma2) * (wofz(self.za) + wofz(self.zb))
        self.alpha = (self.w**2 + 1j * self.w * self.gamma3)**.5
        self.za = (self.alpha - self.w3) / (2**.5 * self.sigma3)
        self.zb = (self.alpha + self.w3) / (2**.5 * self.sigma3)
        self.eps += 1j * (np.pi)**.5 * self.f3 * self.wp**2 / (2**1.5 * self.alpha * self.sigma3) * (wofz(self.za) + wofz(self.zb))
        self.alpha = (self.w**2 + 1j * self.w * self.gamma4)**.5
        self.za = (self.alpha - self.w4) / (2**.5 * self.sigma4)
        self.zb = (self.alpha + self.w4) / (2**.5 * self.sigma4)
        self.eps += 1j * (np.pi)**.5 * self.f4 * self.wp**2 / (2**1.5 * self.alpha * self.sigma4) * (wofz(self.za) + wofz(self.zb))
        self.alpha = (self.w**2 + 1j * self.w * self.gamma5)**.5
        self.za = (self.alpha - self.w5) / (2**.5 * self.sigma5)
        self.zb = (self.alpha + self.w5) / (2**.5 * self.sigma5)
        self.eps += 1j * (np.pi)**.5 * self.f5 * self.wp**2 / (2**1.5 * self.alpha * self.sigma5) * (wofz(self.za) + wofz(self.zb))

        return self.eps


    def pd_bb(self):
    
        """
        Calculates the refractive index of Pd.
    
        Parameters
        ----------
        wavelegnth : numpy.array
            Wavelegnth of light [nm].
    
        Returns
        -------
        eps : numpy.array
            Dielectric permittivity (complex representation).
    
        Info
        ----
        - Brendel and Bormann (1992) An infrared dielectric function model for
          amorphous solids.
        - Rakic (1998) Optical properties of metallic films for
          vertical-cavity optoelectronic devices.
          
        Parameters
        ----------
        - h (eVs),
        - c (m/s),
        - w (eV),
        - wp (eV),
        - gamma0 (eV),
        - gamma1 (eV),
        - gamma2 (eV),
        - gamma3 (eV),
        - gamma4 (eV),
        - sigma0 (eV),
        - sigma1 (eV),
        - sigma2 (eV),
        - sigma3 (eV),
        - sigma4 (eV),
        - omegap (eV),
        - w is the Faddeeva function,
        """

        # Constants:
        self.wavelength_m = self.wavelength_array * 1e-9       
        self.h = 4.135667516 * 1e-15
        self.c = 299792458
        self.w = self.h * self.c / self.wavelength_m
        self.wp = 9.72
        self.f0 = 0.330
        self.gamma0 = 0.009
        self.f1 = 0.769
        self.gamma1 = 2.343
        self.w1 = 0.066
        self.sigma1 = 0.694
        self.f2 = 0.093
        self.gamma2 = 0.497
        self.w2 = 0.502
        self.sigma2 = 0.027
        self.f3 = 0.309
        self.gamma3 = 2.022
        self.w3 = 2.432
        self.sigma3 = 1.167
        self.f4 = 0.409
        self.gamma4 = 0.119
        self.w4 = 5.987
        self.sigma4 = 1.331

        # Equations:
        self.omegap = self.f0**.5 * self.wp
        self.eps = 1 - self.omegap**2 / (self.w * (self.w + 1j * self.gamma0))
        self.alpha = (self.w**2 + 1j * self.w * self.gamma1)**.5
        self.za = (self.alpha - self.w1) / (2**.5 * self.sigma1)
        self.zb = (self.alpha + self.w1) / (2**.5 * self.sigma1)
        self.eps += 1j * (np.pi)**.5 * self.f1 * self.wp**2 / (2**1.5 * self.alpha * self.sigma1) * (wofz(self.za) + wofz(self.zb))
        self.alpha = (self.w**2 + 1j * self.w * self.gamma2)**.5
        self.za = (self.alpha - self.w2) / (2**.5 * self.sigma2)
        self.zb = (self.alpha + self.w2) / (2**.5 * self.sigma2)
        self.eps += 1j * (np.pi)**.5 * self.f2 * self.wp**2 / (2**1.5 * self.alpha * self.sigma2) * (wofz(self.za) + wofz(self.zb))
        self.alpha = (self.w**2 + 1j * self.w * self.gamma3)**.5
        self.za = (self.alpha - self.w3) / (2**.5 * self.sigma3)
        self.zb = (self.alpha + self.w3) / (2**.5 * self.sigma3)
        self.eps += 1j * (np.pi)**.5 * self.f3 * self.wp**2 / (2**1.5 * self.alpha * self.sigma3) * (wofz(self.za) + wofz(self.zb))
        self.alpha = (self.w**2 + 1j * self.w * self.gamma4)**.5
        self.za = (self.alpha - self.w4) / (2**.5 * self.sigma4)
        self.zb = (self.alpha + self.w4) / (2**.5 * self.sigma4)
        self.eps += 1j * (np.pi)**.5 * self.f4 * self.wp**2 / (2**1.5 * self.alpha * self.sigma4) * (wofz(self.za) + wofz(self.zb))

        return self.eps


    def __repr__(self):
        return repr(f'Dielectric Permittivity of Metals at wavelength: ' + str(self.wavelength_array) + ' nm')



class n_Other():
    """
    Generate the refractive index of unusual materials.

    Attributes
    ----------
    wavelength_array : numpy.array
        Wavelength / colour of light [nm].

    Methods
    -------
    - ta2o5_ld
    """

    def __init__(self, wavelength_array,):
        """Initialize; see class doc."""
        self.wavelength_array = wavelength_array


    def ta2o5_ld(self):
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
        
        Info
        ----
        - Bright et al. (2013), Infrared optical properties of amorphous and
          nanocrystalline Ta2O5 thin films.
        """
    
        self.wavelength_um = self.wavelength_array * 1e-9 * 1e6  # wavelength (nm), convert to (m), then convert to (um)
        self.wavelength_cm = self.wavelength_array * 1e-9 * 1e2

        self.a = 2.06
        self.b = 0.025
        self.e_inf = (self.a + self.b / (self.wavelength_um)**2)**2

        self.w_p0 = 649              # in cm-1
        self.gamma_0 = 6.5 * 1e5
        self.w_per_cm = 1 / self.wavelength_cm
        self.t2 = self.w_p0**2 / (self.w_per_cm**2 - (1j * self.gamma_0 * self.w_per_cm))

        self.w_j = np.array([266, 500, 609, 672, 868, 3020])
        self.w_pj = np.array([1040, 573, 634, 408, 277, 373])
        self.gamma_j = np.array([188, 112, 88, 43, 113, 652])

        self.t3 = []
        for i in range(len(self.w_j)):
            self.t3.append(self.w_pj[i]**2 / (self.w_j[i]**2 - self.w_per_cm**2 - 1j * self.gamma_j[i] * self.w_per_cm))
        self.t_sum = self.t3[0] + self.t3[1] + self.t3[2] + self.t3[3] + self.t3[4] + self.t3[5]
        self.eps = self.e_inf - self.t2 + self.t_sum
        self.n = (self.eps**.5).real

        return self.n


    def __repr__(self):
        return repr(f'Refractive Index of unususal materials at wavelength: ' + str(self.wavelength_array) + ' nm')

