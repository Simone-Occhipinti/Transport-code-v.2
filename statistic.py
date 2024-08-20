#statistic

import numpy as np
import scipy as sp
from math import pi
import numpy.random as rnd
import geometry as geo
import particle as prt

class tally:
    def __init__(self, ll=tuple, nn=int, log=bool, domain=geo.mesh, geotp=str):
        self.geometrytype = geotp
        self.space_range = domain.discretization
        if log is True:
            self.energy_range = np.logspace(ll[0],ll[1],nn)
        else:
            self.energy_range = np.linspace(ll[0],ll[1],nn)
        self.estimator = np.array([np.zeros(nn) for __ in range(len(domain.discretization))])
        self.iter = 0
    
    def newiter(self):
        self.iter+=1

    def newinteraction(self,pp=prt.particle,medium=geo.domain):
        if len(self.space_range)>1:
            spaceindex = np.where(pp.position.distance(self.geometrytype)<=self.space_range)[0][0]
        else:
            spaceindex = 0
        energyindex = np.where(pp.energy <= self.energy_range)[0][0]
        index = medium.findposition(pp.position)
        self.estimator[spaceindex][energyindex] += pp.weight/medium[index].macro_xs_scattering(pp.energy)
        
def scattering_step_function(energy_in=float, angle_in=geo.direction, alfa=float, nn=int):
    energy_range = np.logspace(np.log10(1E-5),np.log10(2E7),nn)
    mu_range = np.linspace(-1,1,nn)
    teta_range = np.arccos(mu_range)
    phi_range = np.linspace(0,2*pi,nn)
    out = np.zeros(nn)
    for ii in range(nn):
        if energy_range[ii]<=energy_in and energy_range[ii]>=alfa*energy_in:
            out[ii] += 1/(4*pi)/(1-alfa)/energy_in
    matrix_out = [[out for __ in range(nn)]for __ in range(nn)]
    return [phi_range, teta_range, energy_range, matrix_out]

def sample_newenergy_newangle(pp=prt.particle,medium=geo.domain,sfun=None):
    index = medium.findposition(pp.position)
    isotope = pp.sample_isotope(medium)
    alfa = medium.materials.composition[isotope].alpha
    precision = 100
    matrix = sfun(pp.energy,pp.direction,alfa,precision)
    if pp.type == 'neutron':
        out = 0
        while out == 0:
            rho_phi = 2*pi*rnd.rand()
            rho_mu = 2*rnd.rand()-1
            rho_teta = np.arccos(rho_mu)
            rho_en = pp.energy((1-alfa)+alfa*rnd.rand())
            rho = rnd.rand()
            index_teta = np.where(rho_teta<=matrix[1])[0][0]
            index_phi = np.where(rho_phi<=matrix[0])[0][0]
            index_en = np.where(rho_en<=matrix[2])[0][0]
            if rho<=matrix[3][index_phi][index_teta][index_en]:
                out+=1
        return (rho_phi,rho_teta,rho_en)
    else:
        out = 0
        possible_teta = np.linspace(0,pi,precision)
        possible_phi = np.linspace(0,2*pi,precision)
        possible_energy = np.logspace(np.log10(1E-5),np.log10(2E7),precision)
        dx = np.diff(possible_teta, append=possible_teta[-1])  # Append the last element to maintain array length
        dy = np.diff(possible_phi, append=possible_phi[-1])
        dz = np.diff(possible_energy, append=possible_energy[-1])
        # Create 3D grids for each interval array
        DX, DY, DZ = np.meshgrid(dx, dy, dz, indexing='ij')
        # Compute the 3D weight array
        weights = DX * DY * DZ
        LL = np.zeros(precision,precision,precision)
        for ii in range(len(possible_teta)):
            for jj in range(len(possible_phi)):
                for zz in range(len(possible_energy)):
                    II = scattering_step_function(possible_energy[zz],geo.direction(possible_phi[jj],possible_teta[ii]),alfa,precision)
                    LL += II[3]*weights[ii][jj][zz]                                        
        while out == 0:
            rho_phi = 2*pi*rnd.rand()
            rho_mu = 2*rnd.rand()-1
            rho_teta = np.arccos(rho_mu)
            rho_dir = geo.direction(rho_phi,rho_teta,False)
            rho_en = pp.energy*(1+(1/alfa-1)*rnd.rand())
            rho = rnd.rand()
            matrix1 = sfun(rho_en,rho_dir,alfa,precision)
            index_en = np.where(rho_en<=matrix1[2])[0][0]
            index_teta = np.where(rho_teta<=matrix1[1])[0][0]
            index_phi = np.where(rho_phi<=matrix1[0])[0][0]
            check = medium.materials[index].macro_xs_scattering(pp.energy)/LL[index_phi][index_teta][index_en]*matrix1[3][index_phi][index_teta][index_en]
            if rho <= check:
                out+=1
        return (rho_phi,rho_teta,rho_en)
            