#statistic

import numpy as np
import scipy as sp
from math import pi
import numpy.random as rnd

def rejection(fun=np.array, var=np.array, nn=int):
    if len(fun) is not len(var):
        return 'length_error'
    if np.trapz(fun) is not 1:
        return 'error_not_normalized'
    aa, bb = var[0], var[-1]
    max = np.max(fun)
    out = []
    for ii in range(nn):
        rho = 0
        while rho == 0:
            pp1, pp2 = aa + (bb-aa)*rnd.rand(), max*rnd.rand()
            index = np.where(pp1>=var)
            if index is 1:
                index+=1
            refernce = np.interp(pp1, [var[index-1], var[index]], [fun[index-1], fun[index]])
            if pp2 <= refernce:
                out.append(pp2)
    return out

             
