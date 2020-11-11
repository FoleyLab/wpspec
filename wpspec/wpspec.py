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


class wavefunction:
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


    def pib_eigenfunction(self):
      x = np.linspace(0, self.L, self.grid_points)
      fx = np.sqrt( 2. / self.L ) * np.sin(self.n * np.pi * x / self.L)
      return x, fx

if __name__ == "__main__":
    # Do something if this file is invoked on its own
    print(canvas())



