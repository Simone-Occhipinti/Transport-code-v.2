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
    def distance(self, type=str):
        if type is geometry_type[1]:
            dd = np.sqrt(self.x**2 + self.y**2 + self.z**2)
        elif type is geometry_type[0]:
            dd = self.z
        else:
            dd = 'geometry_error'
        return dd
    
def sumpos(pp=point,qq=point):
    xx = pp.x + qq.x
    yy = pp.y + qq.y
    zz = pp.z + qq.z
    out = point(xx,yy,zz)
    return out

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

class mesh:
    def __init__(self, nDim=int, lDim=None, nInt=None, homo=None, vect=None):
        self.dimension = nDim
        if nDim > 0:
            if homo == True:
                self.discretization = np.linspace(lDim[0], lDim[1], nInt)
                self.delta = np.abs(lDim[1]-lDim[0])/nInt
            else:
                self.discretization = vect
                self.delta = np.average([vect[ii+1]-vect[ii] for ii in range(len(vect)-1)])
        else:
            self.discretization = 0

class domain:
    def __init__(self, disct=mesh, mat=list[tuple], simtype=str, en=tuple,log=bool):
        if log == True:
            self.energyrange = np.logspace(np.log10(en[0]),np.log10(en[1]),en[2])
        else:
            self.energyrange = np.linspace(en[0],en[1],en[2])
        self.geometrytype = simtype
        self.mesh = disct
        self.materials = []
        for ii in range(len(mat)):
            self.materials.append(mat[ii][0])
        self.materialposition = []
        for ii in range(len(mat)):
            self.materialposition.append((mat[ii][1],mat[ii][2]))
            
    def findposition(self, pp=point): # torna l'indice corrispettivo del materiale
        vett = []
        for ii in self.materialposition:
            vett.append(ii[1])
        index = np.where(pp.distance(self.geometrytype) <= vett)[0][0]
        return index
    
        


        