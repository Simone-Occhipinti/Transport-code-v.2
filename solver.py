#solver

import numpy as np
import scipy as sp
import geometry
import material
import particle
import numpy.random as rnd

# reference energy
Eref = 2E6
# particle type
particle_type = 'neutron'
# number of stories
NN = 1E4

# nuclear data input
Cnat_total = np.loadtxt('cross_sections_Janis\\C\\Cnat_total.csv', delimiter=';',skiprows=3)
Cnat_scattering = np.loadtxt('cross_sections_Janis\\C\\Cnat_scattering.csv', delimiter=';', skiprows=3)

Onat_total = np.loadtxt('cross_sections_Janis\\O\\O_total.csv', delimiter=';',skiprows=3)
Onat_scattering = np.loadtxt('cross_sections_Janis\\O\\O_scattering.csv', delimiter=';', skiprows=3)

U238_total = np.loadtxt('cross_sections_Janis\\U238\\U238_total.csv', delimiter=';',skiprows=3)
U238_scattering = np.loadtxt('cross_sections_Janis\\U238\\U238_scattering.csv', delimiter=';', skiprows=3)
U238_fission = np.loadtxt('cross_sections_Janis\\U238\\U238_fission.csv', delimiter=';', skiprows=3)

U235_total = np.loadtxt('cross_sections_Janis\\U235\\U235_total.csv', delimiter=';',skiprows=3)
U235_scattering = np.loadtxt('cross_sections_Janis\\U235\\U235_scattering.csv', delimiter=';', skiprows=3)
U235_fission = np.loadtxt('cross_sections_Janis\\U235\\U235_fission.csv', delimiter=';', skiprows=3)

# isotope composition

carbon = material.isotope(6,12,1.1E23,Cnat_total[:,0],Cnat_total[:,1],Cnat_scattering[:,1])
oxigen = material.isotope(8,16,1.4E21,Onat_total[:,0],Onat_total[:,1],Onat_scattering[:,1])
uranium238 = material.isotope(92,238,1.76E22,U238_total[:,0],U238_total[:,1],U238_scattering[:,1],U238_fission[:,1],2.5)
uranium235 = material.isotope(92,235,4.4E21,U235_total[:,0],U235_total[:,1],U235_scattering[:,1],U235_fission[:,1],2.5)

# material composition
core_list = [uranium238, uranium235, oxigen]
core = material.material(core_list)
reflector_list = [carbon]
reflector = material.material(reflector_list)

# domain definition
n_interval = 100
LL = (0,50)
distribution = [('core',0,30),('reflector',30,50)]

mesh = geometry.mesh(n_interval,LL,distribution,True)
