#particle

import numpy as np
import numpy.random as rnd
import scipy as sp
import geometry as geo
import statistic
import material

particle_type = ['neutron', 'adjunction']
neutron_mass = 1.675E-27

class particle:
    def __init__(self, pos=geo.point, dir=geo.direction, ee=float, tp=str):
        if tp is particle_type[0] or tp is particle_type[1]:
            self.type = tp
        self.position = pos
        self.direction = dir
        self.energy = ee
        @property
        def velocity(self):
            out = np.sqrt(2*self.energy/neutron_mass)
            return out 
        self.weight = 1

class source:
    def __init__(self, range=np.array, prod=np.array, intens=float ,map=None):
        self.distribution = prod
        self.range = range
        if map is None:
            self.map = []
        else:
            self.map = map
        self.intensity = intens

    @property
    def new_particle(self):
        out = statistic.rejection(self.distribution, self.range,1)
        return out
    
    def watt(aa=float, bb=float):
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

class scattering_function:
    def __init__(self, range_en=np.array, range_ang=np.array, fun=np.matrix):
        self.distribution = fun
        self.energy_range = range_en
        self.angle_range = range_ang
    
def sample_free_flight(pp=particle, mat=material.material):
    rho = rnd.rand()