import numpy as np
import cell_class as cc
import matplotlib.pyplot as plt

def simulate_wind():


    Gamma = 2/3
    gamma = 5/3

    """will assume constant mass to start with"""
    T_corona = 10**6 # temperature in Kelvin at the closest cell

    G = 6.674*(10**-11) # N m**2/kg**2
    mass_sun = 1.989*(10**30) # kg

    radius_sun = 6.968*(10**8) # radius in meters

    corona_inner = radius_sun*0.85
    corona_outer = radius_sun


    """looking up plasma, we just missed the factor of n"""
    """but if we're treating it as an ideal gas, maybe go with 
        other form of specific heat"""

    # mass of a hydrogen atom
    mass_of_hydrogen_kg = 1.67*(10**-27) # kg
    mass_of_hydrogen_ev = 931.4*(10**6)  # eV

    cp_hyd_stp = 14701 # J/K at STP, specific heat

    # boltzmann constant

    kb_J = 1.38*(10**-23)  # J/K
    kb_ev = 8.617*(10**-5) # eV/K

    mass_of_cell= 10**8 # kg, conservative estimate
    """because all cells have the same mass, no need to 
        make an eqn for node mass"""

    num_cells = 1000
    num_time_steps = 10

    cell_length = (corona_outer-corona_inner)/num_cells

    print(cell_length)

    """create cells for calculations"""

    cell_velocity = np.zeros(num_cells) # not super important
    cell_pressure = np.zeros(num_cells)
    cell_energy   = np.zeros(num_cells)
    cell_volume   = np.zeros(num_cells)
    cell_density  = np.zeros(num_cells)

    node_area     = np.zeros(num_cells+1) # want a closing node
    node_velocity = np.zeros(num_cells+1)
    node_distance = np.zeros(num_cells+1)

    node_distance[0] = corona_inner
    node_area[0] = 4*np.pi*(node_distance[0]**2)

    """give initial values to the cells"""
    for i in range(num_cells):
        cell_energy[i] = (3/2)*kb_J*T_corona

        node_distance[i+1] = node_distance[i] + cell_length
        node_area[i+1] = 




    return


def initialize_arrays():
    """function to initialize my cells' values
        this function will do it for sets of arrays"""

    return

"""BOUNDARY CONDITION!!!
    Inner cell wall does not move!!!"""