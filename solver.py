import numpy as np
import scipy as sp
import geometry
import material
import particle
import numpy.random as rnd
import matplotlib.pyplot as plt

# geometry type
geometry_type = 'sphere'
# reference energy
Eref = 6E6
Erange = (1E-5,2E7,200)
# particle type
particle_type = 'neutron'
# number of stories
NN = 6E5
# particle source
pos_source = geometry.point(0,0,0)
p_source = particle.source(pos_source,1.0,'Watt')

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
#core_list = [uranium238, uranium235, oxigen]
core_list = [oxigen]
core = material.material(core_list,Erange)
reflector_list = [carbon]
reflector = material.material(reflector_list,Erange)

# domain definition
n_interval = 100
LL = (0,100)
n_dimensions = 1
distribution = [(core,0,50),(reflector,50,100)]
start = np.log10(1E-5)
stop = np.log10(2E7)
energy_range = (start,stop)

mesh = geometry.mesh(1,LL,n_interval,True)
domain = geometry.domain(mesh,distribution,geometry_type,(1E-5,2E7,10000))
counter = particle.tally((energy_range),100,True,mesh,geometry_type)

n_stories = NN
wd = 10
wgt_min = 1/wd
wgt_max = wd
particle_squeue = []
jj = 0

# MC Code
while counter.iter <= n_stories:
    # genero una nuova generazione
    if len(particle_squeue) == 0:
        counter.newiter()
        if jj >= 100:
            print(counter.iter)
            jj-=100
        else:
            jj+=1
        eof = 1
        # new particle's energy
        e0 = p_source.newparticle()
        w0 = 1
        # particle production is always isotropic
        phi0 = 2*np.pi*rnd.random()
        mu0 = 2*rnd.random()-1
        teta0 = np.arccos(mu0)
        pos0 = geometry.point(0,0,0)
        dir0 = geometry.direction(phi0,teta0,False)
        nn = particle.particle(pos0,dir0,e0,w0,particle_type,geometry_type)
    # generazione precedente
    else:
        nn = particle_squeue.pop(0)
        eof = 1
    # interactions simulation
    while eof > 0:
        # distanza percorsa prima dell'interazione
        nn.position = nn.sample_freeflight(domain)
        if nn.position.distance(nn.geometry_type) >= LL[1] or nn.position.distance(nn.geometry_type) <= LL[0]:
            eof = 0
        else:
            # take note of the interaction
            counter.newinteraction(nn,domain)
            # create new angle and energy
            new_data = nn.sample_newenergy_newangle(domain)
            nn.direction, nn.energy = new_data[0],new_data[1]
            # generate the new weight
            nn.new_weight(domain)
            # check for russian roulette
            if nn.weight<=wgt_min:
                rho = rnd.rand()
                if rho<=1/wd:
                    nn.weight *= wd
                else:
                    eof = 0
            # energy check
            if nn.energy<=Erange[0] or nn.energy>=Erange[1]:
                eof = 0
            # check for splitting
            if nn.weight>=wgt_max:
                N = nn.weight/wgt_max
                if N == int(N):
                    for ii in range(N-1):
                        particle_squeue.append(particle.particle(nn.position, nn.direction, nn.energy, nn.weight/N, nn.type, nn.geometry_type))
                    nn.weight *= 1/N
                else:
                    D = N - int(N)
                    if rnd.rand() <= 1-D:
                        for ii in range(int(N)-1):
                            particle_squeue.append(particle.particle(nn.position, nn.direction, nn.energy, nn.weight/int(N), nn.type, nn.geometry_type))
                        nn.weight *= 1/int(N)
                    else:
                        for ii in range(int(N)):
                            particle_squeue.append(particle.particle(nn.position, nn.direction, nn.energy, nn.weight/(int(N)+1), nn.type, nn.geometry_type))
                        nn.weight *= 1/(int(N)+1)
# computing sample average
avg, var = counter.sample_average((1E-5,0.625,2E7))
sigma = np.sqrt(var/counter.iter)

# plotting results
for ii in avg:
    plt.plot(counter.space_range[1:],ii)
plt.xlabel('position [cm]')
if nn.type == 'neutron':
    plt.title('Direct transport problem')
    plt.ylabel(r'$\Phi$ [n/cm3/s/eV]')
elif nn.type == 'adjunction':
    plt.title('Adjoint transport problem')
    plt.ylabel(r'$\Psi$ [-/cm3/s/eV]')
else:
    plt.ylabel('errore')
plt.yscale('log')
plt.show()
