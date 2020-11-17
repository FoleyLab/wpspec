"""
wpspec.py
A python package for simulating light and matter.

Handles the primary functions
"""
import numpy as np
from matplotlib import pyplot as plt
class Quantum:
    def __init__(self, args):
        if 'quantum_state' in args:
            self.n = args['quantum_state']
        else:
            self.n = 1 
        self.L = 1
        self.grid_points = 100
        self.x = np.linspace(0,self.L, self.grid_points)
        self.Psi = np.zeros(self.grid_points, dtype=complex)
        self.V = np.ones(self.grid_points)
        
    def split_operator(self):
        self.Psi = self.Psi * self.V
        return 1
            
            
class pib(Quantum):
    def __init__(self, args):
        Quantum.__init__(self,args)

    def eigenfunction(self):
        """ 
        function to compute the energy eigenfunction of the pib and
        store to attribute self.psi
        
        Parameters
        ----------
        self
        
        Returns
        -------
        self.psi : float array
            Energy eigenfunction of the PIB 
        """
        self.Psi = np.sqrt( 2. / self.L ) * np.sin(self.n * np.pi * self.x / self.L)
        return 1
    
    def plot_eigenfunction(self):
        plt.plot(self.x, self.Psi, 'red', label='Psi')
        plt.legend()
        plt.show()
    
    def delta_potential(self):
        self.V = np.zeros(self.grid_points)
        self.V[int(self.grid_points/2)] = 1
        plt.plot(self.x, self.V, 'blue', label='Potential')
        plt.legend()
        plt.show()
        return 1
        
        

'''
def canvas(with_attribution=True):
    """
    Placeholder function to show example docstring (NumPy format)

    Replace this function and doc string for your own project

    Parameters
    ----------
    with_attribution : bool, Optional, default: True
        Set whether or not to display who the quote is from

    Returns
    -------
    quote : str
        Compiled string including quote and optional attribution
    """
    quote = "When Zarathustra was thirty years old, he left his home and the lake of his home, and went into the Four Seasons Total Landscaping. There he enjoyed his spirit and solitude, and for ten years did not weary of it."
    if with_attribution:
        quote += "\n\t- Adapted from Anna Hughes' Twitter"
    return quote

""" parent class will be the quantum class """
class quantum:
    def __init__(self, args):
        if 'propagator' in args:
            self.method = args['propagator']
        else:
            self.method = 'split_operator'
        if 'grid_points' in args:
            self.grid_points = args['grid_points']
        else:
            self.grid_points = 100
        if 'box_length' in args:
            self.L = args['box_length']
        else:
            self.L = 1
        
        print(" Going to propagate using ",self.method)
        ### get the x-grid
        self.x = np.linspace(0, self.L, self.grid_points)
        ## make the potential a parabola centered at L/2
        self.V = 0.5 * (self.x - self.L/2)**2
        ## make initial wavefunction a PIB ground-state
        self.n = 1
        self.Psi = np.sqrt(2/self.L) * np.sin( self.n * np.pi * self.x / self.L)
        
    def split_operator(self):
        self.Psi = self.Psi * self.V
        return 1
        

        

class pib(quantum):
    def __init__(self, args):
        if 'quantum_state' in args:
            self.n = args['quantum_state']
        else:
            self.n = 1
        if 'box_length' in args:
            self.L = args['box_length']
        else:
            self.L = 1
        if 'particle_mass' in args:
            self.M = args['particle_mass']
        else:
            self.M =  1
        
        """ array of x-values """
        self.x = np.linspace(0, self.L, self.grid_points)
        
        """ reduced plancks constant in atomic units """
        self.hbar = 1

    def eigenfunction(self):
        """ 
        function to compute the energy eigenfunction of the pib and
        store to attribute self.psi
        
        Parameters
        ----------
        self
        
        Returns
        -------
        self.psi : float array
            Energy eigenfunction of the PIB 
        """
        self.Psi = np.sqrt( 2. / self.L ) * np.sin(self.n * np.pi * self.x / self.L)
        return 1
  
    def eigenvalue(self):
        """ 
        function to compute the energy eigenvalue of the pib and
        store to attribute self.E
        
        Parameters
        ----------
        self
        
        Returns
        -------
        self.E : float
            Energy eigenvalue of the PIB 
        """
        self.E = self.n ** 2 * np.pi **2 * self.hbar ** 2  / (2 * self.M * self.L **2) 
        

class pir:
    def __init__(self, args):
        if 'grid_points' in args:
            self.grid_points = args['grid_points']
        else:
            self.grid_points = 100
        if 'quantum_state' in args:
            self.m = args['quantum_state']
        else:
            self.m = 1

        if 'particle_mass' in args:
            self.M = args['particle_mass']
        else:
            self.M =  1
        if 'bond_length' in args:
            self.r = args['bond_length']
        else:
            self.r = 1.
        
        """ array of angles """    
        self.theta = np.linspace(0, np.pi, self.grid_points)
        
        """ moment of inertia """
        self.I = self.M * self.r **2
        
        """ reduced planck's constant in atomic units """
        self.hbar = 1

    def eigenfunction(self):
        """ 
        function to compute the energy eigenfunction of the pib and
        store to attribute self.psi
        
        Parameters
        ----------
        self
        
        Returns
        -------
        self.psi : float array
            Energy eigenfunction of the PIR 
        """
        ci = 0+1j
        self.Psi = np.sqrt( 1. / (2 * np.pi) ) * np.exp(-ci * self.m * self.theta)
        return 1
  
    def eigenvalue(self):
        """ 
        function to compute the energy eigenvalue of the pib and
        store to attribute self.E
        
        Parameters
        ----------
        self
        
        Returns
        -------
        self.E : float
            Energy eigenvalue of the PIB 
        """
        self.E = self.m ** 2 * self.hbar **2 /(2 * self.I)
        return 1
    
    
if __name__ == "__main__":
    # Do something if this file is invoked on its own
    print(canvas())

'''

