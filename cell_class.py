import numpy as np

class Cell:
    """base class for a 1D cell, basically a struct at the moment"""

    mass_hyd = 1.67*(10**-27) # kg

    kb_J = 1.38 * (10 ** -23)  # J/K


    def __init__(self, mass,x_1,x_2):
        # from research online, approx 10**36 particles taken away per second
        self.mass = mass # mass in kg
        self.x_1 = x_1 # inner coordinate position
        self.x_2 = x_2 # outer coordinate position
        self.length = self.x_2 - self.x_1 # length in meters
        self.num_particles = self.mass/self.mass_hyd

        # starting all cells with the same energy
        self.initial_internal_energy = (3/2)*self.num_particles*self.kb_J*(10**6)




class Sph_Cell(Cell):
    """a spherical based cell"""


    def __init__(self,mass,x_1,x_2):
        Cell.__init__(self,mass,x_1,x_2)
        self.volume = (4/3)*np.pi*(x_2**3 - x_1**3)
        self.inner_surface_area = 4*np.pi*(x_1**2)
        self.outer_surface_area = 4*np.pi*(x_2**2)
        self.particle_density = self.num_particles/self.volume
        self.mass_density = self.mass/self.volume


class Cyl_Cell(Cell):
    """a cylindrical cell"""


    def __init__(self,mass,x_1,x_2,height):
        Cell.__init__(self,mass,x_1,x_2)
        self.height = height
        self.volume = height*np.pi*(x_2**2 - x_1**2)
        self.inner_surface_area = height*2*np.pi*x_1
        self.outer_surface_area = height*2*np.pi*x_2


class Sqr_Cell(Cell):
    """a cubic cell"""


    def __init__(self,mass,x_1,x_2,area):
        Cell.__init__(self,mass,x_1,x_2)
        self.volume = area*self.length


"""Should I make a node class?"""
