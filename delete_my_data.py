"""
TMM Model - Kretschmann configuration
"""

import numpy as np
import matplotlib.pyplot as plt
import DispDielectrics as dispD
import DispMetals as dispM
import DispOther as dispO

#%%
n_BK7, e_Ag, e_Pd, n_ta2o5, n_sio2M, n_sio2PY, rs, Rs, rp, Rp, l = ([] for i in range(11))

#l_array = np.linspace(450, 2000, 50)
l_array = [633]                                 # list or linspace??
t_array = np.arange(1, 90, 1)                      # list or linspace?

d = [0, 30, 100, 3, 0]            # Thickness of layers (in order from core / prism)
N = len(d)                       # Number of layers in the Multilayer

f_l = 0.01      # fibre length
f_d = 200e-6    # fibre diameter
NA = 0.22
num_int_points = 106  # number of sample points used to compute the Gaus-Legendre quadrature


# Include all wavelengths:
for i in range(len(l_array)):
    l.append(l_array[i])

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
        
        # Across All Layers in the Multilayer
        e, qs, qp, beta_j = ([] for i in range(4))
        for per in range(N):
            e.append((np.array(n[per]).real**2) - (np.array(n[per]).imag**2) + 1j * (2 * np.array(n[per]).real * np.array(n[per]).imag))
        
        for k in range(1, N-1):
#            e.append((np.array(n[k]).real**2) - (np.array(n[k]).imag**2) + 1j * (2 * np.array(n[k]).real * np.array(n[k]).imag))
            qs.append(np.sqrt(e[k] - n[0]**2 * np.sin(t_rad)**2))
            qp.append(np.sqrt(e[k] - (n[0]**2 * np.sin(t_rad)**2)) / e[k])
            beta_j.append((2 * np.pi / l[i]) * d[k]  * np.sqrt(e[k] - n[0]**2 * np.sin(t_rad)**2))
            
            m11s = np.cos(np.array(beta_j[k-1]))
            m12s = -1 * (1j / np.array(qs[k-1])) * np.sin(np.array(beta_j[k-1]))
            m21s = -1 * (1j * np.array(qs[k-1])) * np.sin(np.array(beta_j[k-1]))
            m22s = np.cos(np.array(beta_j[k-1]))
            Ms = Ms * np.matrix([[m11s, m12s], [m21s, m22s]])

            m11p = np.cos(np.array(beta_j[k-1]))
            m12p = -1 * (1j / np.array(qp[k-1])) * np.sin(np.array(beta_j[k-1]))
            m21p = -1 * (1j * np.array(qp[k-1])) * np.sin(np.array(beta_j[k-1]))
            m22p = np.cos(np.array(beta_j[k-1]))
            Mp = Mp * np.matrix([[m11p, m12p], [m21p, m22p]])
        
        qs0 = np.array(np.sqrt(e[0] - n[0]**2 * np.sin(t_rad)**2))
        qsN = np.array(np.sqrt(e[-1] - n[0]**2 * np.sin(t_rad)**2))
        a_s = (qs0 * (Ms[0,0] + (Ms[0,1] * qsN)) - (Ms[1,0] + (Ms[1,1] * qsN)))
        b_s = (qs0 * (Ms[0,0] + (Ms[0,1] * qsN)) + (Ms[1,0] + (Ms[1,1] * qsN)))
        rs.append(a_s / b_s)
        Rs.append(np.abs(rs[t])**2)

        qp0 = np.array(np.sqrt(e[0] - n[0]**2 * np.sin(t_rad)**2) / e[0])
        qpN = np.array(np.sqrt(e[-1] - n[0]**2 * np.sin(t_rad)**2) / e[-1])
        a_p = (qp0 * (Mp[0,0] + (Mp[0,1] * qpN)) - (Mp[1,0] + (Mp[1,1] * qpN)))
        b_p = (qp0 * (Mp[0,0] + (Mp[0,1] * qpN)) + (Mp[1,0] + (Mp[1,1] * qpN)))
        rp.append(a_p / b_p)
        Rp.append(np.abs(rp[t])**2)



        # RSN_0(ind_theta)=power(rs*conj(rs),Nref);
        # RPN_0(ind_theta)=power(rp*conj(rp),Nref);
        # dP0(ind_theta) = (n_core^2)*(sin(theta_rad))*(cos(theta_rad))/((1-((n_core^2)*(cos(theta_rad)^2)))^2);
        

#%%






#%%
# plt.close('all')
fig = plt.figure(figsize=[8,6])
ax = plt.subplot(111)
ax.plot(t_gauss, Rp, 'b')

plt.xlabel("Incident Angle (\degree)")
plt.ylabel("Reflectance")
plt.title("SPR")
fig.legend((
        'm',
        ), loc='upper right')
plt.show()








