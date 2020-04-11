""" Implementation of numerical integration schemes.

Author: Daniel Rothenberg <darothen@mit.edu>
Version: February 13, 2015

"""

import abc
import numpy as np

class Integrator:
    """ This is an abstract class providing a template for 
    integration routine logic.

    All integration routines need to implement an `integrate`
    method, which march a system of ODEs forward in time over 
    discrete timesteps from some initial conditions, given a
    function which implements the RHS of that ODE system.

    """

    __metaclass__ = abc.ABCMeta
    
    @staticmethod
    @abc.abstractmethod
    def integrate(self, func, y0, t):
        """Performs the integration"""
        
    def __call__(self, func, y0, t):
        return self.integrate(func, y0, t)        

class EulerIntegrator(Integrator):
    
    def integrate(self, func, y0, t):
        out_y = [np.asarray(y0), ]
        for i, ti in enumerate(t[:-1]):
            y = out_y[-1]
            delta_t = t[i+1] - ti
            new_y = y + delta_t*func(y, ti)
            out_y.append(new_y)
        return np.array(out_y)
    
class HeunIntegrator(Integrator):
    
    def integrate(self, func, y0, t):
        out_y = [np.asarray(y0), ]
        for i, ti in enumerate(t[:-1]):
            y = out_y[-1]
            delta_t = t[i+1] - ti
            
            y_tilde = y + delta_t*func(y, ti)
            new_y = y + \
                (delta_t/2.)*(  func(y, ti)                \
                              + func(y_tilde, ti+delta_t) )
            out_y.append(new_y)
        return np.array(out_y)
    
class RK4Integrator(Integrator):
    
    def integrate(self, func, y0, t):
        out_y = [np.asarray(y0), ]
        for i, ti in enumerate(t[:-1]):
            y = out_y[-1]
            delta_t = t[i+1] - ti
            
            k1 = func(y, ti)
            k2 = func(y + 0.5*k1*delta_t, ti + 0.5*delta_t)
            k3 = func(y + 0.5*k2*delta_t, ti + 0.5*delta_t)
            k4 = func(y + k3*delta_t, ti + delta_t)
            
            new_y = y + (delta_t/6.)*(k1 + 2.*k2 + 2.*k3 + k4)
            out_y.append(new_y)
        return np.array(out_y)
