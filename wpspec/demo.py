# -*- coding: utf-8 -*-
"""
Created on Thu Nov 19 10:26:53 2020

@author: foleyj10
"""

from wpspec import pib
from matplotlib import pyplot as plt
import numpy as np
from matplotlib import animation


params = {'quantum_state': 5, 'box_length': 3}
#parent = quantum(params)
#print(parent.grid_points)
wf = pib(params)

wf.eigenfunction()
wf.eigenvalue()
print(wf.E)
#wf.plot_eigenfunction()

wf.derivatives()

T_Psi = -0.5 * wf.Psi_pp

#plt.plot(wf.x, wf.E * wf.Psi, 'blue', label='analytic')
#plt.plot(wf.x, T_Psi, 'r--', label='Derivative')
#plt.show()


#from matplotlib import animation, rc
#from IPython.display import HTML

# First set up the figure, the axis, and the plot element we want to animate
fig, ax = plt.subplots()
plt.close()


### parameters for plot
ax.set_xlim(( 0, wf.L))
ax.set_ylim((-1, 1))

line, = ax.plot([], [], lw=2)

# initialization function: plot the background of each frame
def init():
    line.set_data([], [])
    return (line,)

N_time = 100

# animation function. This is called sequentially  
def animate(i):
    #wf.propagate()
    line.set_data(wf.x, np.real(wf.Psi))
    return (line,)
  
anim = animation.FuncAnimation(fig, animate, init_func=init,
                             frames=N_time, interval=100, blit=True)

plt.show()