import numpy as np
import cell_class as cc
import matplotlib.pyplot as plt
import sys


def explicit_simulate_wind():


    Gamma = 2/3
    gamma = 5/3

    """will assume constant mass to start with"""
    T_corona = 10**6 # temperature in Kelvin at the closest cell

    G = 6.674*(10**-11) # N m**2/kg**2
    mass_sun = 1.989*(10**29) # kg

    radius_sun = 6.968*(10**8) # radius in meters

    corona_inner = radius_sun*0.85
    corona_outer = radius_sun


    """looking up plasma, we just missed the factor of n"""
    """but if we're treating it as an ideal gas, maybe go with 
        other form of specific heat"""

    # mass of a hydrogen atom
    mass_of_hydrogen_kg = 1.67*(10**-27) # kg
    mass_of_hydrogen_ev = 931.4*(10**6)  # eV

    avogadro = 6.022*(10**23)
    R = 8.314
    cp_hyd_stp = 14701 # J/K at STP, specific heat

    # boltzmann constant

    kb_J = 1.38*(10**-23)  # J/K
    kb_ev = 8.617*(10**-5) # eV/K

    #mass_of_cell= 10**8 # kg, conservative estimate

    """here 10**40 is the density of 10**15 particles/m**3 times the volume
        I've chosen to leave the mass constant for the moment"""
    mass_of_cell = (10**40)*mass_of_hydrogen_kg

    num_moles = (mass_of_cell/mass_of_hydrogen_kg)/avogadro

    """The number of cells and time step length are 
        mixed together by a ratio, 
        more cells => shorter length"""

    num_cells = 100
    time_step_length = 1
    num_time_steps = 1000


    cell_length = (corona_outer-corona_inner)/num_cells

    print(cell_length)

    """create cells for calculations"""

    cell_pressure = np.zeros(num_cells)
    cell_energy   = np.zeros(num_cells)
    cell_volume   = np.zeros(num_cells)
    cell_density  = np.zeros(num_cells)
    cell_position = np.zeros(num_cells)
    cell_particles = np.zeros(num_cells)
    cell_mass = np.zeros(num_cells)
    cell_speed_of_sound = np.zeros(num_cells)
    cell_massdensity = np.zeros(num_cells)

    node_area     = np.zeros(num_cells+1) # want a closing node
    node_velocity = np.zeros(num_cells+1) # these will stay zero until after the first time step
    node_position = np.zeros(num_cells+1) # should probably be called position
    node_gravity = np.zeros(num_cells+1)
    node_mass = np.zeros(num_cells+1)

    """matrices for debugging"""
    cell_pressure_matrix = np.zeros(shape=[num_time_steps+1,num_cells])
    cell_energy_matrix = np.zeros(shape=[num_time_steps+1,num_cells])
    cell_volume_matrix = np.zeros(shape=[num_time_steps+1,num_cells])
    cell_density_matrix = np.zeros(shape=[num_time_steps+1,num_cells])
    cell_position_matrix = np.zeros(shape=[num_time_steps+1,num_cells])
    cell_mass_matrix = np.zeros(shape=[num_time_steps+1,num_cells])
    cell_speed_of_sound_matrix = np.zeros(shape=[num_time_steps+1,num_cells])
    cell_massdensity_matrix = np.zeros(shape=[num_time_steps+1,num_cells])

    node_velocity_matrix = np.zeros(shape=[num_time_steps+1,num_cells+1])
    node_position_matrix = np.zeros(shape=[num_time_steps+1,num_cells+1])
    node_area_matrix = np.zeros(shape=[num_time_steps+1,num_cells+1])
    node_gravity_matrix = np.zeros(shape=[num_time_steps+1,num_cells+1])
    node_mass_matrix = np.zeros(shape=[num_time_steps+1,num_cells+1])


    node_position[0] = corona_inner
    node_area[0] = 4*np.pi*(node_position[0]**2)
    node_gravity[0] = -G*mass_sun/(node_position[0]**2)

    node_position_matrix[:,0] = corona_inner
    node_area_matrix[:,0] = node_area[0]


    """give initial values to the cells"""
    for i in range(num_cells):
        node_position[i+1] = node_position[i] + cell_length
        node_area[i+1] = 4*np.pi*(node_position[i+1]**2)
        node_gravity[i+1] = -G*mass_sun/(node_position[i+1]**2)

        node_velocity[i+1] = 10**5 # 100 km/s

        cell_volume[i] = (4/3)*np.pi*(node_position[i+1]**3
                                      - node_position[i]**3)

        cell_position[i] = (node_position[i] + node_position[i + 1]) / 2

        #cell_mass[i] = mass_of_cell*np.exp(-cell_position[i]/cell_position[0])
        cell_mass[i] = mass_of_cell * (cell_position[0] / cell_position[i])**2

        if i < num_cells-1 :
            node_mass[i+1] = (cell_mass[i]+cell_mass[i+1])/2
        #cell_density[i] = mass_of_cell/cell_volume[i]
        cell_density[i] = cell_mass[i]/cell_volume[i]
        cell_massdensity[i] = cell_mass[i]/cell_volume[i]

        cell_particles[i] = cell_mass[i] / mass_of_hydrogen_kg


        cell_energy[i] = (3 / 2) * cell_particles[i]*kb_J * T_corona




        #cell_density[i] = 1 / cell_volume[i]
        #cell_pressure[i] = (gamma-1)*cell_density[i]*cell_energy[i]

        cell_pressure[i] = (num_moles*R*T_corona)/cell_volume[i]

        cell_speed_of_sound[i] = np.sqrt(gamma*(cell_pressure[i]/cell_massdensity[i]))

        node_position_matrix[0,i+1] = node_position[i+1]
        node_area_matrix[0,i+1] = node_area[i+1]
        cell_volume_matrix[0,i] = cell_volume[i]
        cell_pressure_matrix[0,i] = cell_pressure[i]
        cell_energy_matrix[0,i] = cell_energy[i]
        cell_density_matrix[0,i] = cell_density[i]
        node_gravity_matrix[0,i+1] = node_gravity[i+1]
        node_velocity_matrix[0,i+1] = node_velocity[i+1]
        cell_speed_of_sound_matrix[0,i] = cell_speed_of_sound[i]
        cell_position_matrix[0,i] = cell_position[i]



    node_mass[0] = cell_mass[0]/2
    node_mass[num_cells] = cell_mass[num_cells-1]/2

    """initialization complete"""


    """make some initial plots to help verify how the system evolves"""
    #plt.plot(cell_position/radius_sun,cell_pressure,'go')
    #plt.title('initial pressure vs position')
    #plt.legend()
    #plt.show()

    #plt.plot(cell_position / radius_sun, cell_energy, 'go')
    #plt.title('initial energy vs position')
    #plt.legend()
    #plt.show()

    #plt.plot(node_position/radius_sun,node_velocity,'go')
    #plt.title('initial velocity vs position')
    #plt.legend()
    #plt.show()


    """now onto updating during timesteps"""

    for i in range(num_time_steps):
        """maybe put zeroth cell changes here"""
        for j in range(num_cells):
            cell_energy_matrix[i , j] = cell_energy[j]
            delta_energy = evolve_energy(velocity=node_velocity,pressure=cell_pressure,
                                         area=node_area,d_time=time_step_length,count=j,limit=num_cells,
                                         energy_matrix=cell_energy_matrix,time_step=i+1)

            #print(delta_energy)
            node_position[j + 1] = node_position[j + 1] + node_velocity_matrix[i, j + 1] * time_step_length

            node_velocity[j+1] = evolve_velocity(velocity=node_velocity,pressure=cell_pressure,
                                                 d_time=time_step_length,area=node_area,mass=node_mass,
                                                 gravity=node_gravity,j=j,limit=num_cells)

            node_area[j+1] = 4*np.pi*(node_position[j+1]**2)
            node_gravity[j+1] = -G*mass_sun/(node_position[j+1]**2)

            cell_energy[j] = cell_energy[j] + delta_energy
            cell_volume[j] = (4/3)*np.pi*(node_position[j+1]**3 - node_position[j]**3)
            """experimenting """

            cell_massdensity[j] = cell_mass[j] / cell_volume[j]

            cell_density[j] = 1 / cell_volume[j]
            cell_pressure[j] = (gamma-1)*cell_density[j]*cell_energy[j]

            cell_position[j] = (node_position[j] + node_position[j + 1]) / 2

            #cell_density[j] = mass_of_cell / cell_volume[j]
            #cell_pressure[j] = (num_moles * R * T_corona) / cell_volume[j]

            node_velocity_matrix[i+1,j+1] = node_velocity[j+1]
            node_position_matrix[i+1,j+1] = node_position[j+1]
            node_area_matrix[i+1,j+1] = node_area[j+1]
            cell_energy_matrix[i+1,j] = cell_energy[j]
            cell_volume_matrix[i+1,j] = cell_volume[j]
            cell_density_matrix[i+1,j] = cell_density[j]
            cell_pressure_matrix[i+1,j] = cell_pressure[j]
            cell_position_matrix[i+1,j] = cell_position[j]
            cell_speed_of_sound_matrix[i+1,j] = np.sqrt(gamma*(cell_pressure[j]/cell_massdensity[j]))

            #if(node_velocity[j+1] < 0):
                #print("negative velocity! likely error",j,i)
                #return
        plot_title = 'Time Evolution for timestep = %i' % i
        save_title = 'plots/vel_pos_%i.pdf' % i
        if i%(100) == 0 :
            plt.plot(node_position / radius_sun, node_velocity/1000,'bo')
            plt.title(plot_title)
            plt.xlabel('distance in solar radii')
            plt.ylabel('speed in km/s')
            plt.savefig(save_title)
            plt.show()

    #print(node_position[50],node_velocity[50])

    plt.plot(node_position_matrix[900,:] / radius_sun, node_velocity_matrix[900,:]/1000,'bo')
    plt.plot(cell_position_matrix[900,:]/radius_sun,cell_speed_of_sound_matrix[900,:]/1000,'go')
    plot_title = 'Velocity vs Distance'
    plt.title(plot_title)
    plt.xlabel('distance in solar radii')
    plt.ylabel('speed in km/s')
    plt.savefig('plots/final_vel_pos.pdf')
    plt.show()

    np.set_printoptions(threshold=np.inf)

    #print('\npressure \n',cell_pressure_matrix[0:1000:100,:])
    #print('\npressure \n', cell_pressure_matrix[0, :])
    #print('\nmass\n',cell_particles)
    print('\nvelocity \n',node_velocity_matrix[num_time_steps-15,:]/1000)
    #print('\nposition \n',node_position_matrix[:,:]/radius_sun)
    print('\nenergy\n',cell_energy_matrix[num_time_steps-15,:])
    #print('\nvolume\n',cell_volume_matrix[0,:])
    #print('\ndensity\n',cell_density_matrix[0,:])
    #print('\narea\n',node_area_matrix[0,:])
    #print('\ngravity\n',node_gravity_matrix[0,:])

    #plt.plot(node_position/radius_sun,node_velocity/1000,'bo')
    #plt.show()

    """total_energy_start = sum( node_velocity_matrix[0, :] ** 2) + sum(cell_energy_matrix[0, :])
    total_energy_end = sum(node_velocity_matrix[num_time_steps - 1, :] ** 2) + sum(
        cell_energy_matrix[num_time_steps - 1, :])
    print(total_energy_start)
    print(total_energy_end)
    print("End Total energy is ", total_energy_end / total_energy_start, " from start energy")"""

    total_energy_start = sum( 0.5 * node_mass * node_velocity_matrix[1,:] ** 2) + sum(cell_energy_matrix[1, :])
    total_energy_end = sum( 0.5 * node_mass * node_velocity_matrix[900, :] ** 2) + sum(
        cell_energy_matrix[900, :])
    print("End Total energy is ", total_energy_end / total_energy_start, " from start energy")

    return



def evolve_energy(velocity,pressure,area,d_time,count,limit,energy_matrix,time_step):

    cell_limit = limit-1
    node_limit = limit+1

    cell_pos = count
    node_pos = count


    if ( cell_pos < cell_limit and cell_pos > 0):
        """p_diff1 = pressure[cell_pos] - pressure[cell_pos-1]
        p_diff2 = pressure[cell_pos] - pressure[cell_pos+1]
        #p_area = p_diff1*area[node_pos] + p_diff2*area[node_pos+1]
        d_energy1 = p_diff1*area[node_pos]*velocity[node_pos]*d_time
        d_energy2 = p_diff2 * area[node_pos+1] * velocity[node_pos+1] * d_time"""

        """minus sign on energy2 because expanding"""
        """d_energy = d_energy1 - d_energy2"""

        d_energy1 = velocity[node_pos]*d_time*area[node_pos]*(pressure[cell_pos]-pressure[cell_pos-1])
        d_energy2 = velocity[node_pos+1]*d_time*area[node_pos+1]*(pressure[cell_pos]-pressure[cell_pos+1])

        d_energy = -d_energy1 - d_energy2

    # beginning boundary term
    if(cell_pos == 0):
        #p_diff = pressure[0] - pressure[1]

        """minus sign because expanding"""
        #d_energy = -velocity[1]*d_time*area[1]*p_diff
        d_energy = -velocity[1]*d_time*area[1]*pressure[0]

    # last cell boundary term
    if(cell_pos == cell_limit):
        """p_diff1 = pressure[cell_pos] - pressure[cell_pos-1]
        p_diff2 = pressure[cell_pos]
        p_area1 = p_diff1 * area[node_pos]
        p_area2 = p_diff2 * area[node_pos + 1]

        d_energy1 = velocity[count]*d_time*p_area1
        d_energy2 = velocity[cell_pos+1]*d_time*p_area2"""

        """minus sign on energy2 because it is expanding outward"""

        d_energy1 = velocity[node_pos]*d_time*area[node_pos]*pressure[cell_pos]
        d_energy2 = velocity[node_pos+1]*d_time*area[node_pos+1]*pressure[cell_pos]

        d_energy = d_energy1 - d_energy2

    #if(energy_matrix[time_step,cell_pos] + d_energy < 0):
        #print("negative energy at ",cell_pos,time_step,d_energy,energy_matrix[time_step,cell_pos])


    #if count == 1:
        #print('count1')
        #print(p_diff1)
        #print(p_diff2)
        #print(p_area)
        #print('\ndelta energy = ',d_energy,'\n')

    if d_energy == np.inf or d_energy == np.nan :
        print("overflow energy at ",time_step)
        #sys.exit()

    return d_energy

def evolve_velocity(velocity,pressure,d_time,area,mass,gravity,j,limit):

    cell_limit = limit-1
    node_limit = limit+1

    if(j<cell_limit):
        p_diff = pressure[j] - pressure[j+1]
        grav_t = gravity[j+1]*d_time
        p_vel = p_diff*area[j+1]*d_time/mass[j+1]

        new_vel = velocity[j+1] + p_vel + grav_t

        delta_v = p_vel + grav_t
        #print('change in velocity = ',delta_v)


    if(j==cell_limit):
        p_diff= pressure[j] - pressure[j]/2
        grav_t = gravity[j+1]*d_time
        p_vel = p_diff*area[j+1]*d_time/mass[j+1]

        new_vel = velocity[j+1] + p_vel + grav_t

        #if (pressure[j]) < 10**-5*pressure_matrix[0,j] :
        #    new_vel = velocity[j+1] + grav_t

    if new_vel == np.inf or new_vel == np.nan :
        print("overflow velocity at ",j)
        #sys.exit()


    return new_vel




explicit_simulate_wind()
