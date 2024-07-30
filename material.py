#material

import numpy as np
import scipy as sp
from math import pi

class isotope:
    def __init__(self, atm_num=int, mas_num=int, dens=float, energy=np.array, tot=np.array, scatt=np.array, fiss=None, nu=None):
        self.atomic_number = atm_num
        self.mass_number = mas_num
        self.micro_xs_total = tot/1E24
        self.micro_xs_scattering = scatt/1E24
        self.micro_xs_absorption = self.micro_xs_total-self.micro_xs_scattering
        if fiss is not None:
            self.micro_xs_fission = fiss
            self.neutron_emitted = nu
        else:
            self.micro_xs_fission = np.zeros(len(energy))
            self.fission_spectrum = np.zeros(len(energy))
            self.neutron_emitted = 0
        self.micro_xs_capture = self.micro_xs_absorption-self.micro_xs_fission
        self.atomic_density = dens
        self.energy = energy

class material:
    def __init__(self, parts=list[isotope]):
        nn = len(parts)
        self.energy = parts[0].energy
        self.macro_total = np.zeros(len(parts[0].micro_xs_total))
        self.macro_scattering = np.zeros(len(parts[0].micro_xs_total))
        self.macro_absorption = np.zeros(len(parts[0].micro_xs_total))
        self.macro_fission = np.zeros(len(parts[0].micro_xs_total))
        self.macro_capture = np.zeros(len(parts[0].micro_xs_total))
        for ii in range(nn):
            self.macro_total += parts[ii].micro_xs_total*parts[ii].atomic_density
            self.macro_scattering += parts[ii].micro_xs_scattering*parts[ii].atomic_density
            self.macro_absorption += parts[ii].micro_xs_absorption*parts[ii].atomic_density
            self.macro_fission += parts[ii].micro_xs_fission*parts[ii].atomic_density
            self.macro_capture += parts[ii].micro_xs_capture*parts[ii].atomic_density
    def macro_xs_total(self, en=float):
        if en >= self.energy[-1]:
            out = self.macro_total[-1]
        elif en <= self.energy[0]:
            out = self.macro_total[0]
        else:
            index = np.where(self.energy>en)[0][0]
            out = np.interp(en, [self.energy[index-1],self.energy[index]], [self.macro_total[index-1],self.macro_total[index]])
        return out
    def macro_xs_scattering(self, en=float):
        if en >= self.energy[-1]:
            out = self.macro_scattering[-1]
        elif en <= self.energy[0]:
            out = self.macro_scattering[0]
        else:
            index = np.where(self.energy>en)[0][0]
            out = np.interp(en, [self.energy[index-1],self.energy[index]], [self.macro_scattering[index-1],self.macro_scattering[index]])
        return out
    def macro_xs_absorption(self, en=float):
        if en >= self.energy[-1]:
            out = self.macro_absorption[-1]
        elif en <= self.energy[0]:
            out = self.macro_absorption[0]
        else:
            index = np.where(self.energy>en)[0][0]
            out = np.interp(en, [self.energy[index-1],self.energy[index]], [self.macro_absorption[index-1],self.macro_absorption[index]])
        return out
    def macro_xs_fission(self, en=float):
        if en >= self.energy[-1]:
            out = self.macro_fission[-1]
        elif en <= self.energy[0]:
            out = self.macro_fission[0]
        else:
            index = np.where(self.energy>en)[0][0]
            out = np.interp(en, [self.energy[index-1],self.energy[index]], [self.macro_fission[index-1],self.macro_fission[index]])
        return out
    def macro_xs_capture(self, en=float):
        if en >= self.energy[-1]:
            out = self.macro_capture[-1]
        elif en <= self.energy[0]:
            out = self.macro_capture[0]
        else:
            index = np.where(self.energy>en)[0][0]
            out = np.interp(en, [self.energy[index-1],self.energy[index]], [self.macro_capture[index-1],self.macro_capture[index]])
        return out