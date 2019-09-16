from amuse.datamodel import Particles
from amuse.units import *
from amuse.io import read_set_from_file, write_set_to_file
import numpy as np

pos,mass = np.load('particles.npy')

print len(pos),mass

particles = Particles( len(pos) )
particles.x = pos[:,0] | units.km
particles.y = pos[:,1] | units.km
particles.z = pos[:,2] | units.km
particles.vx = 0.0 | units.km/units.s
particles.vy = 0.0 | units.km/units.s
particles.vz = 0.0 | units.km/units.s
particles.mass = [mass] * len(pos) | units.MSun

print particles
write_set_to_file(particles, 'InitialState.amuse', 'amuse')


