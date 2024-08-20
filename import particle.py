import particle
import numpy as np
import geometry
import numpy.random as rnd
from solver import domain

phi0 = rnd.rand()*2*3.14
teta0 = rnd.rand()*3.14
pos_source = geometry.point(0,0,0)
p_source = particle.source(pos_source,1.0,'Watt')
e0 = p_source.newparticle()
w0 = 1
pos0 = geometry.point(0,0,0)
dir0 = geometry.direction(phi0,teta0,False)
nn = particle.particle(pos0,dir0,e0,w0,'neutron','sphere')

nn.sample_freeflight(domain)
new_data = nn.sample_newenergy_newangle(domain)
nn.direction, nn.energy = new_data[0],new_data[1]
nn.new_weight(domain)
print([nn.position.x,nn.position.y,nn.position.z,nn.energy])