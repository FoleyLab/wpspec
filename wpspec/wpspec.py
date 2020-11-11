"""
wpspec.py
A python package for simulating light and matter.

Handles the primary functions
"""
import numpy as np

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


class pib:
    def __init__(self, args):
        if 'grid_points' in args:
            self.grid_points = args['grid_points']
        else:
            self.grid_points = 100
        if 'quantum_state' in args:
            self.n = args['quantum_state']
        else:
            self.n = 1
        if 'box_length' in args:
            self.L = args['box_length']
        else:
            self.L = 1
        if 'particle_mass' in args:
            self.m = args['particle_mass']
        else:
            self.m =  1
            
        self.x = np.linspace(0, self.L, self.grid_points)


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
        self.psi = np.sqrt( 2. / self.L ) * np.sin(self.n * np.pi * x / self.L)
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
        self.E = self.n ** 2 * np.pi **2 / (2 * self.m * self.L **2) 
        

if __name__ == "__main__":
    # Do something if this file is invoked on its own
    print(canvas())



