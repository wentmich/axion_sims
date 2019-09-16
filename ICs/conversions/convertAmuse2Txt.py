from amuse.lab import *

particle_data = read_set_from_file('InitialState.amuse', 'amuse')

x_pos = particle_data.x.value_in(units.cm)
y_pos = particle_data.y.value_in(units.cm)
z_pos = particle_data.z.value_in(units.cm)

x_vel = particle_data.vx.value_in(units.cm / units.s)
y_vel = particle_data.vy.value_in(units.cm / units.s)
z_vel = particle_data.vz.value_in(units.cm / units.s)

mass  = particle_data.mass.value_in(units.g)

N = len(mass)
print("Total of " + str(N) + " Particles")

f = open("initial_data_in_txt_form.txt", 'w')
for i in range(N):
    f.write(str(x_pos[i]) + '\t' + str(y_pos[i]) + '\t' + str(z_pos[i]) + '\t' + str(x_vel[i]) + '\t' + str(y_vel[i]) + '\t' + str(z_vel[i]) + '\t' + str(mass[i]) + '\n')

f.close()
