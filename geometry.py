#geometry

import numpy as np
import scipy as sp
from math import pi

geometry_type = ['plain', 'sphere']

class point:
    def __init__(self, ics=float, ips=float, zeta=float):
        self.x = ics
        self.y = ips
        self.z = zeta
    @property
    def distance(self, type=str):
        if type is geometry_type[1]:
            dd = np.sqrt(self.x^2 + self.y^2 + self.z^2)
        elif type is geometry_type[0]:
            dd = self.z
        else:
            dd = 'geometry_error'
        return dd

class direction:
    def __init__(self, PHI=float, TE=float, degr=bool):
        if degr is True:
            PHI *= pi/180
            TE *= pi/180
        if TE >= 0 and TE <= pi:
            self.teta = TE
            self.phi = PHI
        elif PHI >= 0 and PHI <= pi:
            self.teta = PHI
            self.phi = TE
        else:
            self.teta = 'errore'
            self.phi = 'errore'

class surface:
    def __init__(self, start=float, tp=str):
        self.position = start
        self.type = tp
    @property
    def distance_from(self, pp=point):
        if self.type is geometry_type[0]:
            return pp.z - self.position
        elif self.type is geometry_type[1]:
            return pp.distance - self.position
        else:
            return 'geometry_error'
        
class volume:
    def __init__(self, iniz=surface, final=surface, tp=str):
        self.surface1 = iniz
        self.surface2 = final
        self.width = final.position - iniz.position
    def isinside(self, pp=point):
        if pp.distance >= np.min([self.surface1.position,self.surface2.position]) and pp.distance <= np.max([self.surface1.position,self.surface2.position]):
            return True
        else:
            return False

class mesh:
    def __init__(self, nn=int, ll=tuple, form=list[tuple], homo=bool, vect=None):
        if homo is True:
            self.discretization = np.linspace(ll[0], ll[1], nn)
            self.delta = np.abs(ll[1]-ll[0])/nn
        else:
            self.discretization = vect
            self.delta = np.average([vect[ii+1]-vect[ii] for ii in range(len(vect)-1)])
        self.composition = []
        jj=0
        for ii in range(nn):
            self.composition.append(form[jj][0])
            if self.discretization[ii] >= form[jj][2]:
                jj+=1

        

        