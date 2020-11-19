"""
wpspec.py
A python package for simulating light and matter.

Handles the primary functions
"""
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline

from matplotlib import pyplot as plt
class Quantum:
    def __init__(self, args):
        if 'quantum_state' in args:
            self.n = args['quantum_state']
        else:
            self.n = 1 
        if 'box_length' in args:
            self.L = args['box_length']
        else:
            self.L = 1
        if 'time_step' in args:
            self.dt = args['time_step']
        else:
            self.dt = 0.01
            
        self.grid_points = 100
        self.x = np.linspace(0,self.L, self.grid_points)
        self.Psi = np.zeros(self.grid_points, dtype=complex)
        self.Psi_p = np.zeros(self.grid_points, dtype=complex)
        self.Psi_pp = np.zeros(self.grid_points, dtype=complex)
        self.V = np.zeros(self.grid_points)
        self.hbar = 1
        self.m = 1
        
    #def split_operator(self):
    #    self.Psi = self.Psi * self.V
    #    return 1
    
    def derivatives(self):
        """ 
        function to compute the first and second derivatives of the
        wavefunction by interpolating with a cubic spline and using
        the built-in derivative methods of the splines.  The first
        derivative will be stored in the attribute self.Psi_p
        and the second derivative will be stored in the attribute self.Psi_pp
        """
        ci = 0+1j
        ### fit the spline to real and imaginary part of wf
        fr = InterpolatedUnivariateSpline(self.x, np.real(self.Psi))
        fi = InterpolatedUnivariateSpline(self.x, np.imag(self.Psi))
        ### get the derivative of the spline
        fr_p = fr.derivative()
        fi_p = fi.derivative()
        ### get the second derivative of the spline
        fr_pp = fr_p.derivative()
        fi_pp = fi_p.derivative()
        ### store the derivative of the spline to self.Psi_p
        self.Psi_p = fr_p(self.x) + ci * fi_p(self.x)
        ### sore the second derivative of the spline to self.Psi_pp
        self.Psi_pp = fr_pp(self.x) + ci * fi_pp(self.x)
        return 1
    
    def propagate(self):
        """ 
        function to propagate the wavefunction using finite-differences
        with an Euler update 
        """
        ci = 0+1j
        Psi_old = self.Psi 
        #self.derivatives()
        ### kinetic energy operator on wavefunction requires second derivative
        ### store T_hat on Psi to array T_Psi
        T_Psi = self.E * self.Psi
        ### store V_hat on Psi to array V_Psi
        #V_Psi = self.V * self.Psi
        ### Store -i/hbar *  H_hat on Psi to array Psi_dot
        k1 = -ci * T_Psi * self.dt
        self.Psi = Psi_old + k1/2
        
        #self.derivatives()
        T_Psi = self.E * self.Psi       
        k2 = -ci * T_Psi * self.dt
        
        self.Psi = Psi_old + k2/2
        
        #self.derivatives()
        T_Psi = self.E * self.Psi 
        k3 = -ci * T_Psi * self.dt
        
        self.Psi = Psi_old + k3
        
        #self.derivatives()
        #T_Psi = -0.5 * self.Psi_pp
        T_Psi = self.E *  self.Psi 
        k4 = -ci * T_Psi * self.dt
        
        self.Psi = Psi_old + 1/6 * (k1 + 2*k2 + 2*k3 + k4)
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
        self.E = self.n ** 2 * np.pi **2 * self.hbar ** 2  / (2 * self.m * self.L **2) 
    
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

