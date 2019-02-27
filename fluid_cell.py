import numpy as np
import cell_class as cc
import matplotlib.pyplot as plt

def simulate_wind():

    num_cells = 1000
    num_time_steps = 1000
    cell_length = 1000 # meters

    """will assume constant mass to start with"""
    T_initial = 10**6 # temperature in Kelvin at the closest cell

    G = 6.674*(10**-11) # N m**2/kg**2
    mass_sun = 1.989*(10**30) # kg


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


    return


def kinetic_energy():
    


    return

def energy_initial():
    """initial energy of a cell except the boundary cell"""
    """length derivative of ideal gas"""

    return



"""BOUNDARY CONDITION!!!
    The pressure on the innermost wall does not have to be constant, but
    that is the important boundary condition for consistently adding energy 
    to the system
    An adiabatic gas' work would normally go into just expanding the volume,
    and the temperature would decrease, but work is being put into the system
    as it expands, and that must be being converted into kinetic energy"""