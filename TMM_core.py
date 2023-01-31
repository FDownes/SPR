"""
TMM Model - Kretschmann configuration
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from material_dispersion import n_Dielectrics as ndie
from material_dispersion import n_Metals as nmet
from material_dispersion import n_Other as noth


#%%

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


#%%


def get_material_dispersion(wavelength_array, plot=False):
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
    wavelength_array = np.array(wavelength_array)

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
        df = pd.DataFrame({
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
        df = pd.DataFrame({
            'wavelength' : wavelength_array,
            'glass_bk7' :  n_glass_bk7,
            'glass_schott' : n_glass_schott,
            'sio2_m' : n_sio2_m,
            'sio2_py' : n_sio2_py,
            'pd' : n_pd,
            'ag' : n_ag,
            'ta2o5' : n_ta2o5,
            })

    if plot==True:
        plot_dielectrics(df=df)
        plot_metals(df=df)
        plot_ta2o5(df=df)

    return df


#%%
# testing is material dispersion works:
wa = np.linspace(450, 2000, 50)
df = get_material_dispersion(wavelength_array=wa, plot=True)


#%%
wa1 = 633
m1 = get_material_dispersion(wavelength_array=wa1)
wa2 = np.array(633)
m2 = get_material_dispersion(wavelength_array=wa2)


#%%
n_BK7, e_Ag, e_Pd, n_ta2o5, n_sio2M, n_sio2PY, rs, Rs, rp, Rp, l = ([] for i in range(11))

#wavelength_array = np.linspace(450, 2000, 50)
wavelength_array = [633]                                 # list or  linspace??
# t_array = np.arange(1, 90, 1)                      # list or linspace?

d = [0, 30, 100, 3, 0]            # Thickness of layers (in order from core / prism)

f_l = 0.01      # fibre length
f_d = 200e-6    # fibre diameter
NA = 0.22
num_int_points = 106  # number of sample points used to compute the Gaus-Legendre quadrature

def tmm(wavelength_array=633,
        layer_thickness_list=[0, 30, 0],
        materials_list=['bk7', 'ag', 'air'],
        fibre_length=0.01,
        fibre_diameter=200e-6,
        NA=0.22,
        num_int_points=106,
        
        ):

    """
    
    """
    

N = len(d)                       # Number of layers in the Multilayer

# Include all wavelengths:
for i in range(len(wavelength_array)):
    l.append(wavelength_array[i])

    # Generate Material Properties :
    n_BK7.append(dispD.BK7(l[i]))       # Glass prism
    e_Ag.append(dispM.BB_Ag(l[i]))
    e_Pd.append(dispM.BB_Pd(l[i]))
    n_ta2o5.append(dispO.ta2o5(l[i]))
    n_sio2M.append(dispD.sio2M(l[i]))   # SiO2 Fibre Core
    n_sio2PY.append(dispD.sio2PY(l[i])) # SiO2 This Film


    # Will be moved, need to be able to initialize material properties
    # Material in each layer, same order as thicknesses list
    n = [n_sio2M[i], np.sqrt(e_Ag[i]), n_sio2PY[i], np.sqrt(e_Pd[i]), 1]    

    # Optical fibre properties:
    n_core = n_sio2M[i]
    n_cladding = np.sqrt(n_core**2 - NA**2)
    t_critical = np.arcsin(n_cladding / n_core)
    t_upper = 90
    t_lower = 0     # t_critical

    # Gauss Legendre Quadrature:
    sample_points, weights = np.polynomial.legendre.leggauss(num_int_points)
    norm_int =((t_upper - t_lower) / 2) * np.sum(weights)
    t_gauss = ((t_upper - t_lower) / 2) * sample_points + ((t_lower + t_upper) / 2);

    # Include all angles of incidence:
    for t in range(len(t_gauss)):
        theta = t_gauss[t]
        t_rad = theta * np.pi /180

        n_ref = f_l / (f_d * np.tan(t_rad)) # number of fibre reflections
        
        Ms = np.matrix([[1, 0], [0, 1]])
        Mp = np.matrix([[1, 0], [0, 1]])
        
        e, qs, qp, beta_j = ([] for i in range(4))
        for per in range(N):
            # Inner components of the multilayer:
            e.append((np.array(n[per]).real**2) - (np.array(n[per]).imag**2) + 1j * (2 * np.array(n[per]).real * np.array(n[per]).imag))

        # Transfer matrix multilayer setup:
        for k in range(1, N-1):
            # e.append((np.array(n[k]).real**2) - (np.array(n[k]).imag**2) + 1j * (2 * np.array(n[k]).real * np.array(n[k]).imag))
            qs.append(np.sqrt(e[k] - n[0]**2 * np.sin(t_rad)**2))
            qp.append(np.sqrt(e[k] - (n[0]**2 * np.sin(t_rad)**2)) / e[k])
            beta_j.append((2 * np.pi / l[i]) * d[k]  * np.sqrt(e[k] - n[0]**2 * np.sin(t_rad)**2))

            # S-polarization:
            m11s = np.cos(np.array(beta_j[k-1]))
            m12s = -1 * (1j / np.array(qs[k-1])) * np.sin(np.array(beta_j[k-1]))
            m21s = -1 * (1j * np.array(qs[k-1])) * np.sin(np.array(beta_j[k-1]))
            m22s = np.cos(np.array(beta_j[k-1]))
            Ms = Ms * np.matrix([[m11s, m12s], [m21s, m22s]])

            # P-polarization:
            m11p = np.cos(np.array(beta_j[k-1]))
            m12p = -1 * (1j / np.array(qp[k-1])) * np.sin(np.array(beta_j[k-1]))
            m21p = -1 * (1j * np.array(qp[k-1])) * np.sin(np.array(beta_j[k-1]))
            m22p = np.cos(np.array(beta_j[k-1]))
            Mp = Mp * np.matrix([[m11p, m12p], [m21p, m22p]])

        # S-polarization:
        qs0 = np.array(np.sqrt(e[0] - n[0]**2 * np.sin(t_rad)**2))
        qsN = np.array(np.sqrt(e[-1] - n[0]**2 * np.sin(t_rad)**2))
        a_s = (qs0 * (Ms[0,0] + (Ms[0,1] * qsN)) - (Ms[1,0] + (Ms[1,1] * qsN)))
        b_s = (qs0 * (Ms[0,0] + (Ms[0,1] * qsN)) + (Ms[1,0] + (Ms[1,1] * qsN)))
        rs.append(a_s / b_s)
        Rs.append(np.abs(rs[t])**2) # Reflectance

        # P-polarization:
        qp0 = np.array(np.sqrt(e[0] - n[0]**2 * np.sin(t_rad)**2) / e[0])
        qpN = np.array(np.sqrt(e[-1] - n[0]**2 * np.sin(t_rad)**2) / e[-1])
        a_p = (qp0 * (Mp[0,0] + (Mp[0,1] * qpN)) - (Mp[1,0] + (Mp[1,1] * qpN)))
        b_p = (qp0 * (Mp[0,0] + (Mp[0,1] * qpN)) + (Mp[1,0] + (Mp[1,1] * qpN)))
        rp.append(a_p / b_p)
        Rp.append(np.abs(rp[t])**2) # Reflectance



        # RSN_0(ind_theta)=power(rs*conj(rs),Nref);
        # RPN_0(ind_theta)=power(rp*conj(rp),Nref);
        # dP0(ind_theta) = (n_core^2)*(sin(theta_rad))*(cos(theta_rad))/((1-((n_core^2)*(cos(theta_rad)^2)))^2);
        

#%%






#%%
# plt.close('all')
fig = plt.figure()
ax = plt.subplot(111)
ax.plot(t_gauss, Rp, 'b')

plt.xlabel("Incident Angle (\degree)")
plt.ylabel("Reflectance")
plt.title("SPR")
fig.legend((
        'm',
        ), loc='upper right')
plt.show()








