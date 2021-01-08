"""
wpspec.py
A python package for simulating light and matter.

Handles the primary functions
"""
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline
from matplotlib import pyplot as plt
import math
from numpy.polynomial.hermite import *


class Quantum:
    def __init__(self, args):

        if 'box_length' in args:
            self.L = args['box_length']
        else:
            self.L = 10.0
        if 'time_step' in args:
            self.dt = args['time_step']
        else:
            self.dt = 0.05
        if 'system' in args:
            self.system = args['system']
        else:
            self.system = 'harmonic'
            
        if 'grid_points' in args:
            self.grid_points =args['grid_points']
        else:
            self.grid_points = 256
        if 'v_offset' in args:
            self.voffset = args['v_offset']
        else:
            self.voffset = 0.
        if 'wfc_offset' in args:
            self.wfc_offset = args['wfc_offset']
        else:
            self.wfc_offset = 0.
        if 'basis_functions' in args:
            self.n_basis_functions = args['basis_functions']
        else:
            self.n_basis_functions = 100
        

        if self.system == 'harmonic':
            ''' create harmonic object '''
            self.dx = 2 * self.L / self.grid_points
            res = self.grid_points 
            self.x = np.arange(-self.L + self.L / self.grid_points, self.L, self.dx)
            self.dk = np.pi / self.L
            self.k = np.concatenate((np.arange(0, res / 2),
                                 np.arange(-res / 2, 0))) * self.dk
            self.V = 0.5 * (self.x - self.voffset) ** 2
            self.Psi = np.exp(-((self.x - self.wfc_offset) ** 2) / 2, dtype=complex)
            #self.Psi = complex(model.eigenfunction(2))
            self.K = np.exp(-0.5 * (self.k ** 2) * self.dt * 1j)
            self.R = np.exp(-0.5 * self.V * self.dt * 1j)
        elif self.system == 'pib':
            self.x = np.linspace(0, self.L, self.grid_points)
            self.dx = self.x[1]-self.x[0]
            res = self.grid_points 
            self.dk = 2 * np.pi / (res * self.dx)
            self.k = np.concatenate((np.arange(0, res / 2),
                                 np.arange(-res / 2, 0))) * self.dk
            
            self.V = np.zeros_like(self.x)
            #self.V = 0.5 * (self.x - self.voffset) ** 2
            self.Psi = np.sqrt(2/self.L) * np.sin((np.pi * self.x/self.L), dtype=complex)
            self.K = np.exp(-0.5 * (self.k ** 2) * self.dt * 1j)
            self.R = np.exp(-0.5 * self.V * self.dt * 1j)
            
            
        
        ''' some constants and parameters in atomic units '''
        # reduced planck's constant
        self.hbar = 1
        # electron mass
        self.m = 1
        # reduced mass for PIR system or a very light harmonic oscillator!
        self.mu = 1
        # very small spring constant!
        self.k = 1
        # position-space spread for a gaussian wavepacket
        self.sigma = self.L / 5
        # central position value for a gaussian wavepacket
        self.x0 = self.L / 2
        # central momentum for a gaussian wavepacket
        self.k0 = 0.4

        ''' some arrays for finite-difference methods '''
        self.T_matrix = np.zeros((self.grid_points, self.grid_points))
        self.V_matrix = np.zeros((self.grid_points, self.grid_points))
        self.H_matrix = np.zeros((self.grid_points, self.grid_points)) 
        
        ''' arrays for basis set expansion '''
        self.cn = np.zeros(self.n_basis_functions, dtype=complex)
        self.n = np.linspace(1,self.n_basis_functions,self.n_basis_functions)
        self.phi = np.zeros((len(self.x), self.n_basis_functions),dtype=complex)
        self.t_fac = np.zeros(self.n_basis_functions, dtype=complex)
        
        ''' some parameters for model potentials '''
        self.morse_D = 1
        self.morse_a = 1.5
        self.morse_re = 0
        self.periodic_A = 3.25 / 27.211
    
        
    ''' methods for model potentials that might be relevant for various systems '''
    def morse_potential(self):
        ''' method to generate the Morse potential and store to self.V '''
        self.V = self.morse_D * (1 - np.exp(-self.morse_a * (self.x - self.morse_a))) ** 2
        
    def periodic_potential(self):
        ''' Creates a periodic potential perturbing the particle-in-a-box
            system to mimic screened atomic attraction in cyanine dye
            systems as desribed in Figure 2 and Eq. 4 here:
            http://web.ist.utl.pt/berberan/PQF/JCE%202007%20Autschbach.pdf
        
            a = 2.49e-10
            # b constant in SI units related to the length of the terminal moities in 
            # the cyanine dyes
            b = 5.69e-10

           # create an array of L values, where L is the length of the box,
           # for each distinct dye defined by the k array
           L = a * (k + 1) + b
           m = 2*k + 5
           (L - b)/a - 1 = k
           
        '''
        # constant for the average -CH=CH- unit length in SI units
        a = 2.49e-10
        # convert to atomic units
        a /= 5.29e-11
        # constant for the terminal -N moity length in SI units
        b = 5.69e-10
        # convert to atomic units
        b /= 5.29e-11
        # number of -CH=CH- units
        k = (self.L - b)/a - 1
        # integer for number of atom sites
        m = 2 * np.floor(k) + 5
        # pre-factor for the potential in atomic units (3.25 eV in the paper
        # seems a reasonable value)
        self.V = self.A * np.cos( 2 * m * np.pi * self.x / self.L)
        
    ''' methods for numerical methods, including split operator and
        finite difference method '''
        
            
    def build_operator(self):
        self.R = np.exp(-0.5 * self.V * self.dt * 1j)
    
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
    
    def split_op(self):
        

        # Half-step in real space
        self.Psi *= self.R

        # FFT to momentum space
        self.Psi = np.fft.fft(self.Psi)

        # Full step in momentum space
        self.Psi *= self.K

        # iFFT back
        self.Psi = np.fft.ifft(self.Psi)

        # Final half-step in real space
        self.Psi *= self.R

        # Density for plotting and potential
        #density = np.abs(self.Psi) ** 2

        # Renormalizing for imaginary time
        '''
        if par.im_time:
            renorm_factor = sum(density) * par.dx
            opr.wfc /= sqrt(renorm_factor)

        # Outputting data to file. Plotting can also be done in a
        # similar way. This is set to output exactly 100 files, no
        # matter how many timesteps were specified.
        if i % (par.timesteps // 100) == 0:
            filename = "output{}.dat".format(str(i).rjust(5, str(0)))
            with open(filename, "w") as outfile:
                # Outputting for gnuplot. Any plotter will do.
                for j in range(len(density)):
                    template = "{}\t{}\t{}\n".format
                    line = template(par.x[j], density[j].real, opr.V[j].real)
                    outfile.write(line)
            print("Outputting step: ", i + 1)
        '''
    def finite_difference_T_matrix(self):
        ''' check out: http://sites.science.oregonstate.edu/~roundyd/COURSES/ph366/free-particle.html
        '''
        N = self.grid_points
        t0 = self.hbar ** 2 /(2 * self.m * self.dx ** 2)
        for i in range(0,N):
            self.T_matrix[i,i] = 2*t0
            if i<N-1 and i>0:
                self.T_matrix[i,i+1] = -t0
                self.T_matrix[i,i-1] = -t0
            elif i<N-1:
                self.T_matrix[i,i+1] = -t0
            elif i>0:
                self.T_matrix[i,i-1] = -t0
        return 1
    
    def finite_difference_V_matrix(self):
        N = self.grid_points
        for i in range(0,N):
            self.V_matrix[i,i] = self.V[i]
        return 1
    
    def finite_difference_H_matrix(self):
        self.H_matrix = self.T_matrix + self.V_matrix
        return 1
    
    ''' methods related to superpositions and Fourier analysis of
        wavefunctions '''
    
    ### use rectangle rule to integrate a function!
    def integrate(self, f_of_x):
        ''' rectangle rule integral method '''
        # get the width of each rectangle!
        w = self.x[1]-self.x[0]
        # initiate integrand
        integral = 0
        for i in range(1,len(self.x)):
            h = f_of_x[i]
            A = w * h
            integral = integral + A
        return integral
    
    def compute_coefficient(self, basis_state):
        ''' Determine the expansion coefficient for 
            pib eigenfunction n in wavefunction Psi 
        '''
        # product of PIB state _state_ and wavefunction
        integrand = np.conj(basis_state) * self.Psi
        ### get coefficient c_n from the integral of integrand
        return self.integrate(integrand)
    
    def expand_wavefunction_t(self):
        
        self.Psi = self.cn[0] * self.phi[:,0] * self.t_fac[0]
        for i in range(1, len(self.n)):
            self.Psi += self.cn[i] * self.phi[:,i] * self.t_fac[i]
        return 1
    
    ''' methods related to initial quantum states that may be 
        broadly applicable '''
    def gaussian_wavepacket(self):
        ''' method to compute a Gaussian wavepacket defined by 
            a central position value x0 with width sigma 
            and a central momentum value k0 and store it to 
            self.Psi
        '''
        ### need sqrt of -1!
        ci = 0+1j
        ### T1 will be the prefactor that is 1/(sigma * sqrt(2 * pi))
        T1 = 1/(self.sigma * np.sqrt(2 * np.pi))
        ### T2 will be the Gaussian function, exp(-0.5 * ((x-x0)/sigma)^2)
        T2 = np.exp(-0.5 * ((self.x-self.x0)/self.sigma)**2)
        ### T3 will be the complex exponential (aka the plane wave!)
        T3 = np.exp(ci * self.k0 * self.x)
        ### return the product of T1 * T2 * T3
        self.Psi = T1 * T2 * T3
        
    def position_eigenfunction(self, x_val):
        ''' method to approximate a position eigenfunction 
            at x=x_val with a very small width = self.L / 100.
            Note that the <p> = 0 for this model, but the 
            momentum uncertainty will be quite large particularly
            if the box is very small '''
        sig = self.L / 100
        T1 = 1/(sig * np.sqrt(2 * np.pi))
        T2 = np.exp(-0.5 * ((self.x-x_val)/sig)**2)
        self.Psi = T1 * T2 
        

    



class harmonic(Quantum):
    def __init__(self, args):
        Quantum.__init__(self,args)
        
    def eigenfunction(self, state):
            
        w = np.sqrt(self.k/self.mu)
        psi = np.zeros_like(self.x)
            
        herm_coeff = []
            
        for i in range(int(state)):
            herm_coeff.append(int(0))
        herm_coeff.append(int(1))
            
        for i in range(0,len(self.x)):
            psi[i] = math.exp(-self.mu*w*self.x[i]**2/(2*self.hbar)) * hermval((self.mu*w/self.hbar)**0.5 * self.x[i], herm_coeff)
        psi = np.multiply(psi, 1 / (math.pow(2, state) * math.factorial(state))**0.5 * (self.mu*w/(np.pi*self.hbar))**0.25)
        
        return psi
    
    def eigenvalue(self, state):
        return np.sqrt(self.k/self.mu) * (state+(1./2))
    
    def time_factor(self, time_step):
        t = time_step * self.dt
        ci = 0+1j
        En = self.eigenvalue(self.n)
        self.t_fac = np.exp(-ci*En*t) 
    
    def expand_harmonic(self):
        ''' Determine the expansion coefficient for 
            harmonic eigenfunction n in wavefunction Psi 
        '''
        for i in range(0,len(self.n)):
            self.phi[:,i] = self.eigenfunction(self.n[i])
            self.cn[i] = self.compute_coefficient(self.phi[:,i])
        
    
            
            
class pib(Quantum):
    def __init__(self, args):
        Quantum.__init__(self,args)
        

    def eigenfunction(self, n):
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
        psi = np.sqrt( 2. / self.L ) * np.sin(n * np.pi * self.x / self.L)
        return psi
    
    def plot_wavefunction(self):
        plt.plot(self.x, self.Psi, 'red', label='Psi')
        plt.legend()
        plt.show()
        
    def eigenvalue(self, n):
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
        E = n ** 2 * np.pi **2 * self.hbar ** 2  / (2 * self.m * self.L **2)
        return E
    
    def time_factor(self, time_step):
        t = time_step * self.dt
        ci = 0+1j
        En = self.eigenvalue(self.n)
        self.t_fac = np.exp(-ci*En*t) 
        
    
    def delta_potential(self):
        self.V = np.zeros(self.grid_points)
        self.V[int(self.grid_points/2)] = 1
        plt.plot(self.x, self.V, 'blue', label='Potential')
        plt.legend()
        plt.show()
        return 1
    
    def expand_pib(self):
        ''' Determine the expansion coefficient for 
            pib eigenfunction n in wavefunction Psi 
        '''
        for i in range(0,len(self.n)):
            self.phi[:,i] = self.eigenfunction(self.n[i])
            self.cn[i] = self.compute_coefficient(self.phi[:,i])
            
    def stern_gerlach_pib(self,time):
        ''' simulates the Stern-Gerlach experiment for a particle in a box
            system initialized as a Gaussian wavepacket that then undergoes
            position measurement and energy measurement '''
        # free evolution at first!
        if time<200:
            self.time_factor(time)
            self.expand_wavefunction_t()
        # measure position!
        elif time == 200:
            print("Measuring position!")
            # get probability density 
            P = np.real(np.conj(self.Psi) * self.Psi)
            # normalize the probability density
            norm = np.sum(P)
            P_norm = P / norm
            
            # get random position consistent with probability density
            p0 = np.random.choice(self.x, 1, p=P_norm)
            # collapose to that random position
            self.position_eigenfunction(p0)
            # determine coefficients for superposition 
            # of particle-in-a-box energy eigenfunctions 
            # for this position eigenfunction
            self.expand_pib()
            self.time_factor(time - 200)
            self.expand_wavefunction_t()
        elif time < 300:
            self.time_factor(time - 200)
            self.expand_wavefunction_t()
        elif time == 300:
            pn = np.real( np.conj(self.cn) * self.cn )
            norm = np.sum(pn)
            pn_norm = pn / norm 
            nval = np.random.choice(self.n, 1, p=pn_norm)
            self.cn = np.zeros(len(self.cn),dtype=complex)
            self.cn[int(nval[0])-1] = 1+0j
            self.time_factor(time - 300)
            self.expand_wavefunction_t()
        else:
            self.time_factor(time-300)
            self.expand_wavefunction_t()
            



    
class rigid_rotor(Quantum):
    def __init__(self, args):
        Quantum.__init__(self,args)
        self.R = 1
        self.I = self.m * self.R **2
        
        ### quantum numbers
        self.n = np.linspace(-50,50,101)
        self.cn = np.zeros(len(self.n),dtype=complex)
        self.t_fac = np.zeros(len(self.n),dtype=complex) 
        self.x = np.linspace(0, np.pi * 2, self.grid_points)
    

    def eigenfunction(self, n):
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
        ci = 0+1j
        psi = 1/np.sqrt(2*np.pi) * np.exp(ci * n * self.x)
        return psi
    
    def plot_wavefunction(self):
        plt.plot(self.x, self.Psi, 'red', label='Psi')
        plt.legend()
        plt.show()
        
    def eigenvalue(self, n):
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
        E = self.hbar **2 / (2 * self.I) * n ** 2
        return E
    
    def time_factor(self, time_step):
        t = time_step * self.dt
        ci = 0+1j
        En = self.eigenvalue(self.n)
        self.t_fac = np.exp(-ci*En*t) 
        
    
    def delta_potential(self):
        self.V = np.zeros(self.grid_points)
        self.V[int(self.grid_points/2)] = 1
        plt.plot(self.x, self.V, 'blue', label='Potential')
        plt.legend()
        plt.show()
        return 1
    
    def expand_pir(self):
        ''' Determine the expansion coefficient for 
            pib eigenfunction n in wavefunction Psi 
        '''
        for i in range(0,len(self.n)):
            self.phi[:,i] = self.eigenfunction(self.n[i])
            self.cn[i] = self.compute_coefficient(self.phi[:,i])
        

        
        


