#!/bin/env python
#
# use consistent units

import scipy
from scipy import optimize
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import animation

x_range = [-3*np.pi, 3*np.pi]
m=1.
hbar=1.
dt=np.pi/150

#V0 = 1e10
#a = 1e-10
#V0 = -1
#V0 = 0.54
#a = 1.
#E = 1

class Potential(object):
    def __init__(self, V0, a):
        self.V0 = V0
        self.a = a

class WaveFunction(object):
    def __init__(self, E, potential):
        self.E = E
        self.potential = potential
        k1,k2,A1,B1,A2,B2,A3 = self._calc_coefficients(E, potential)
        self.k1 = k1      
        self.k2 = k2
        self.A1 = A1
        self.B1 = B1
        self.A2 = A2
        self.B2 = B2
        self.A3 = A3
        self.T = self._calc_transmission_coefficient()

    # define stationary wave functions pieces
    def psi1_R(self, x): return self.A1 * np.exp(1j*self.k1*x) 
    def psi1_L(self, x): return self.B1 * np.exp(-1j*self.k1*x) 
    def psi1(self, x): return self.psi1_R(x) + self.psi1_L(x)
    def psi2_L(self, x): return self.A2 * np.exp(1j*self.k2*x) # k2 can be imag
    def psi2_R(self, x): return self.B2 * np.exp(-1j*self.k2*x) # k2 can be imag
    def psi2(self, x): return self.psi2_R(x) + self.psi2_L(x)
    def psi3(self, x): return self.A3 * np.exp(1j*self.k1*x) 

    # define actual wave function
    def Psi(self, x, t):
        # may not work for scalars because of this bug:
        # https://github.com/numpy/numpy/issues/5729
        a = self.potential.a
        x = x+0j
        phase = np.exp(-1j*self.E*t/hbar)
        return phase * np.piecewise(x, [x.real<=0, (x.real>0) & (x.real<a), x.real>=a], 
                                    [self.psi1, self.psi2, self.psi3])

    def _calc_coefficients(self, E, potential):
        # coefficients that satisfy continuity of psi and psi'
        # k2 can be real or imag depending on E-V0
        # real case corresponds to plane waves because psi2: e^ix, e^-ix
        V0 = self.potential.V0
        a = self.potential.a
        k1=np.sqrt(2.*m*E)/hbar
        if V0-E<0:
            k2 = np.sqrt(2.*m*(E-V0))/hbar
        if V0-E>=0:
            k2 = 1j*np.sqrt(2.*m*(V0-E))/hbar
        print "k1,k2:          ",k1,k2
        A1 = -np.exp(1j*k1*a)/4/k1/k2 * ((k2-k1)**2*np.exp(1j*k2*a)-(k1+k2)**2*np.exp(-1j*k2*a))
        B1 = np.exp(1j*k1*a)/4/k1/k2 * (k2-k1) * (k1+k2) * (np.exp(1j*k2*a)-np.exp(-1j*k2*a))
        A2 = (k1+k2)/2/k2*np.exp(1j*(k1-k2)*a)
        B2 = (k2-k1)/2/k2*np.exp(1j*(k1+k2)*a)
        A3 = 1.
        
        # normalize so norm of incident wave is 1
        const = A1
        A1 = A1/const
        B1 = B1/const
        A2 = A2/const
        B2 = B2/const
        A3 = A3/const

        print "A1,B1,A2,B2,A3:  ",A1,B1,A2,B2,A3
        return k1,k2,A1,B1,A2,B2,A3

    def _calc_transmission_coefficient(self):
        V0 = self.potential.V0
        a = self.potential.a
        E = self.E
        print "T (analytical, from scherrer)", 1/(1+(V0*V0/4/E/(V0-E))*np.sinh(1j*self.k2*a)**2.)
        return np.absolute(self.A3)**2 / np.absolute(self.A1)**2
       
    def find_a_pseudopotential(self, guess=0):
        # find alternative V0 such that transmission coefficient is unchanged
        V0 = self.potential.V0
        a = self.potential.a
        E = self.E
        def F(V0):
            # taken from Scherrer's analytical equation for transmission coefficient
            # zeros of this function represent V0 that will give transmission const = T
            T = self.T
            #return np.real(T*V0**2*np.sinh(a*np.sqrt(complex(2*m*(V0-E)))/hbar)**2+(T-1)*4*E*(V0-E))
            return np.real(T*V0**2/(4*E*(V0-E))*np.sinh(a*np.sqrt(complex(2*m*(V0-E)))/hbar)**2+(T-1))
        try:
            V0_pseudo = scipy.optimize.broyden1(F, guess, f_tol=1e-14)
            return float(V0_pseudo) # may not work if multiple solutions returned?
        except:
            return None

    def find_some_pseudopotentials(self, Vmin=-5, Vmax=5, steps=50):
        # find all or most of the pseudopotentials in the range Vmin to Vmax
        # steps is number of guesses to check
        V_list = [self.find_a_pseudopotential(guess) for guess in np.linspace(Vmin, Vmax,steps)]
        V_list_valid = [V for V in V_list if V is not None] # remove Nones
        V_list_unique = np.unique(np.array(V_list_valid).round(decimals=8))
        return V_list_unique


def main():
    pot = Potential(V0, a)
    wf = WaveFunction(E, pot)

    print "Transmission coefficient", wf.T
    print wf.k2
    plot_stationary(wf)
    animate_scattering(wf)

def plot_stationary(wave_function):
    a = wave_function.potential.a
    V0 = wave_function.potential.V0
    x = np.linspace(x_range[0], x_range[1], 500)

    plt.figure()

    # plot x,y axes
    plt.axvline(linewidth=0.5, color = 'k')
    plt.axhline(linewidth=0.5, color = 'k')
    
    # manually plot potential
    plt.plot([x_range[0], 0], [0, 0], 'k-', lw=3)
    plt.plot([0, 0], [0, V0], 'k-', lw=3)
    plt.plot([0, a], [V0, V0], 'k-', lw=3)
    plt.plot([a, a], [V0, 0], 'k-', lw=3)
    plt.plot([a, x_range[1]], [0,0], 'k-', lw=3)
    
    # plot wave function
    plt.plot(x, np.real(wave_function.Psi(x, 0.)))
    plt.plot(x, np.imag(wave_function.Psi(x, 0.)))

    plt.show(block=False)
    #raw_input('...')

def animate_scattering(wave_function):
    # all wave functions should have same potential for this to be valid
    a = wave_function.potential.a
    V0 = wave_function.potential.V0
    x = np.linspace(x_range[0], x_range[1], 5000)

    # First set up the figure, the axis, and the plot element we want to animate
    fig = plt.figure()
    ax = plt.axes(xlim=(x_range[0], x_range[1]), ylim=(-2, 2))
    #x = np.linspace(0, 2., 1000)
    lines = []
    for i in range(2):
        lobj = ax.plot([],[],lw=2)[0]
        lines.append(lobj)

    def init():
        # manually plot potential
        ax.plot([x_range[0], 0], [0, 0], 'k-', lw=3)
        ax.plot([0, 0], [0, V0], 'k-', lw=3)
        ax.plot([0, a], [V0, V0], 'k-', lw=3)
        ax.plot([a, a], [V0, 0], 'k-', lw=3)
        ax.plot([a, x_range[1]], [0,0], 'k-', lw=3)

        lines[0].set_data([], [])
        lines[1].set_data([], [])
        return lines[0], lines[1]


    # animation function.  This is called sequentially
    def animate(t_step):
        t = t_step*dt
        #x = np.linspace(0, 2, 1000)
        #Psi_sum = np.sum(1
        Psi_real = np.real(wave_function.Psi(x, t))
        Psi_imag = np.imag(wave_function.Psi(x, t))
        lines[0].set_data(x, Psi_real)
        lines[1].set_data(x, Psi_imag)
        return lines[0], lines[1]
    
    # call the animator.  blit=True means only re-draw the parts that have changed.
    anim = animation.FuncAnimation(fig, animate, init_func=init,
                                   frames=10000, interval=20., blit=True)
    
    plt.title('E: '+str(wave_function.E)+', a: '+str(a)+', V0: '+str(V0))
    plt.show()

if __name__ == "__main__":
    main()
