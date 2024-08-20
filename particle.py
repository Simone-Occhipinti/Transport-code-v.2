import numpy as np
import numpy.random as rnd
import scipy as sp
import geometry as geo
import material
from math import pi

particle_type = ['neutron', 'adjunction']
neutron_mass = 1.675E-27

def rejection(fun=np.array, var=np.array, nn=int):
    if len(fun) != len(var):
        return 'length_error'
    if np.trapz(fun) != 1:
        return 'error_not_normalized'
    aa, bb = var[0], var[-1]
    max = np.max(fun)
    out = []
    for ii in range(nn):
        rho = 0
        while rho == 0:
            pp1, pp2 = aa + (bb-aa)*rnd.rand(), max*rnd.rand()
            index = np.where(pp1>=var)
            if index == 1:
                index+=1
            refernce = np.interp(pp1, [var[index-1], var[index]], [fun[index-1], fun[index]])
            if pp2 <= refernce:
                out.append(pp2)
    return out

class particle:
    def __init__(self, pos=geo.point, dir=geo.direction, ee=float, ww=float , tp=str, geotp=str):
        if tp is particle_type[0] or tp is particle_type[1]:
            self.type = tp
        self.position = pos
        self.direction = dir
        self.energy = ee
        @property
        def velocity(self):
            out = np.sqrt(2*self.energy/neutron_mass)
            return out 
        self.weight = ww
        self.geometry_type = geotp
    
    def sample_freeflight(self, medium=geo.domain):
        med_index = medium.findposition(self.position)
        delta_x = np.sin(self.direction.teta)*np.cos(self.direction.phi)
        delta_y = np.sin(self.direction.teta)*np.sin(self.direction.phi)
        delta_z = np.cos(self.direction.teta)
        delta = geo.point(delta_x,delta_y,delta_z)
        prova = geo.sumpos(self.position,delta)
        if self.position.distance(self.geometry_type) <= prova.distance(self.geometry_type):
            bb = [0]
            bb.append((medium.materialposition[med_index][1]-self.position.distance(self.geometry_type))*medium.materials[med_index].macro_xs_total(self.energy))
            for ii in range(med_index+1,len(medium.materials)):
                rr = medium.materialposition[ii][1]-medium.materialposition[ii][0]
                bb.append(medium.materials[ii].macro_xs_total(self.energy)*rr)
        else:
            bb = []
            for ii in range(med_index):
                rr = medium.materialposition[ii][1]-medium.materialposition[ii][0]
                bb.append(medium.materials[ii].macro_xs_total(self.energy)*rr)
            bb.append((self.position.distance(self.geometry_type)-medium.materialposition[med_index][0])*medium.materials[med_index].macro_xs_total(self.energy))
            bb.append(0)
            bb.reverse()
        BB = np.cumsum(bb)
        rho = rnd.rand()
        eta = -np.log(rho)
        if eta > BB[-1]:
            return geo.point(1E5,1E5,1E5)
        else:
            index = np.where(eta<=BB)[0][0]
            ll = (BB[index-1]-eta)/medium.materials[index-1].macro_xs_total(self.energy)
            delta_x = ll*np.sin(self.direction.teta)*np.cos(self.direction.phi)
            delta_y = ll*np.sin(self.direction.teta)*np.sin(self.direction.phi)
            delta_z = ll*np.cos(self.direction.teta)
            delta = geo.point(delta_x,delta_y,delta_z)
            new_pos = geo.sumpos(self.position,delta)
            return new_pos

    def sample_isotope(self,medium=geo.domain):
        index_med = medium.findposition(self.position)
        nn = len(medium.materials[index_med].composition)
        XS = np.zeros(nn)
        if self.energy <= 1E-5:
            index = 0
        else:
            index = np.where(self.energy<=medium.materials[index_med].energy)[0][0]
        for ii in range(nn):
            XS[ii] += medium.materials[index_med].composition[ii].micro_xs_scattering[index]*medium.materials[index_med].composition[ii].atomic_density/medium.materials[index_med].macro_scattering[index]
        rho = rnd.rand()
        cumulative = np.cumsum(XS)
        out = np.where(cumulative>rho)[0][0]
        return out
    
    def new_weight(self,medium=geo.domain):
        mat_index = medium.findposition(self.position)
        self.weight *= medium.materials[mat_index].macro_xs_scattering(self.energy)/medium.materials[mat_index].macro_xs_total(self.energy)
        if self.type=='adjunction':
            is_index = self.sample_isotope(medium)
            alfa = medium.materials[mat_index].composition[is_index].alpha
            low_i = np.where(self.energy>=medium.materials[mat_index].composition[is_index].energy)[0][0]
            up_i = np.where(self.energy/alfa<=medium.materials[mat_index].composition[is_index].energy)[0][0]
            SS = np.trapz(medium.materials[mat_index].macro_scattering[low_i:up_i]/medium.materials[mat_index].composition[is_index].energy[low_i:up_i],medium.materials[mat_index].composition[is_index].energy[low_i:up_i])
            self.weight *= 1/(1-alfa)/medium.materials[mat_index].macro_xs_scattering(self.energy)*SS

    def sample_newenergy_newangle(self,medium=geo.domain):
        index = medium.findposition(self.position)
        isotope = self.sample_isotope(medium)
        alfa = medium.materials[index].composition[isotope].alpha
        new_phi = 2*pi*rnd.rand()
        new_teta=np.arccos(2*rnd.rand()-1)
        new_dir = geo.direction(new_phi,new_teta,False)
        # energy
        if self.type=='neutron':
            new_energy = self.energy*(alfa+(1-alfa)*rnd.rand())
        elif self.type=='adjunction':
            low_i = np.where(self.energy>=medium.materials[index].composition[isotope].energy)[0][0]
            up_i = np.where(self.energy/alfa<=medium.materials[index].composition[isotope].energy)[0][0]
            SS = np.trapz(medium.materials[index].macro_scattering[low_i:up_i]/medium.materials[index].composition[isotope].energy[low_i:up_i],medium.materials[index].composition[isotope].energy[low_i:up_i])
            exponent = rnd.rand()*SS/medium.materials[index].macro_xs_scattering(self.energy)
            result_log = np.log(self.energy)+exponent
            new_energy = np.exp(result_log)
        return (new_dir,new_energy)

class source:
    def __init__(self, pos=geo.point, intens=float, type = str, range=None, distribution=None,map=None):
        self.position = pos
        self.intensity = intens
        self.type = type
        if type != 'Watt':
            self.distribution = distribution
            self.range = range

    def watt(self, aa, bb):
        kk = 1 + bb/(8*aa)
        LL = (kk + np.sqrt(kk**2 - 1))/aa
        MM = aa*LL-1
        out = 0
        while out<=0:
            xx = -np.log(rnd.rand())
            yy = -np.log(rnd.rand())
            if ((yy-MM*(xx+1))**2)<=bb*LL*xx:
                out += LL*xx
        return out*1E6
    
    def newparticle(self):
        if self.type != 'Watt':
            out = rejection(self.distribution, self.range,1)
        else:
            out = self.watt(0.988,2.249)
        return out

class tally:
    def __init__(self, ll=tuple, nn=int, log=bool, domain=geo.mesh, geotp=str):
        self.geometrytype = geotp
        self.space_range = domain.discretization
        if log is True:
            self.energy_range = np.logspace(ll[0],ll[1],nn)
        else:
            self.energy_range = np.linspace(ll[0],ll[1],nn)
        self.estimator = np.array([np.zeros(nn-1) for __ in range(len(domain.discretization)-1)])
        self.mom2 = np.array([np.zeros(nn-1) for __ in range(len(domain.discretization)-1)])
        self.iter = 0
    
    def newiter(self):
        self.iter+=1

    def newinteraction(self,pp=particle,medium=geo.domain):
        if len(self.space_range)>1:
            spaceindex = np.where(pp.position.distance(self.geometrytype)<=self.space_range)[0][0]-1
        else:
            spaceindex = 0
        energyindex = np.where(pp.energy <= self.energy_range)[0][0]-1
        index = medium.findposition(pp.position)
        self.estimator[spaceindex][energyindex] += pp.weight/medium.materials[index].macro_xs_scattering(pp.energy)
        self.mom2[spaceindex][energyindex] += (pp.weight/medium.materials[index].macro_xs_scattering(pp.energy))**2
    
    def sample_average(self,n_groups=tuple):
        out1 = [np.zeros(len(self.space_range)-1) for __ in range(len(n_groups))]
        out2 = [np.zeros(len(self.space_range)-1) for __ in range(len(n_groups))]
        low_i = int(0)
        up_i = int(0)
        for ii in range(len(n_groups)-1):
            up_i = np.where(n_groups[ii+1]<=self.energy_range)[0][0]
            for jj in range(len(self.space_range)-1):
                out1[ii][jj] += np.sum(self.estimator[jj][low_i:up_i])/self.iter
                out2[ii][jj] += np.sum(self.mom2[jj][low_i:up_i])/self.iter
            low_i = up_i
        out1 = np.array(out1)
        out2 = np.array(out2)
        diffE = np.diff(np.array(n_groups))
        diffE *= 1/2E6
        var = out2-out1**2
        for ii in range(len(diffE)):
            out1[ii] *= 1/diffE[ii]
        if len(self.space_range)>1:
            for ii in out1:
                for jj in range(len(ii)):
                    if self.geometrytype == 'sphere':
                        ii[jj] *= 1/(4/3*pi*(self.space_range[jj+1]**2-self.space_range[jj]**2))
                    else:
                        ii[jj] *= 1/(self.space_range[jj+1]-self.space_range[jj])
        return (out1,var)