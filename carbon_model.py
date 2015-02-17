""" Implementation of Compartmental of Carbon Cycle model.

Author: Daniel Rothenberg <darothen@mit.edu>
Version: February 13, 2015

"""

import numpy as np
import pandas as pd

class CarbonModel(object):
    """ Container class implementing Compartmental Carbon Cycle
    model.

    This class encapsulates the logic for instantiating, integrating,
    and analyzing a multi-component model of the carbon cycle. It is
    designed to simplify the task of repeateadly running the model 
    with different parameters - including emissions functions - and
    for saving and storing the results. Running the model is very simple:

    """

    ## Default Parameters

    beta_2 = 9.4
    beta_3 = 10.2
    gamma  = 62.0
    Gamma  = 198.
    a_d    = 0.230
    a_r    = 1.0

    k_12   = 0.0931
    k_13   = 0.0311
    k_15   = 147.0
    k_21   = 58.*(720.**(-beta_2))
    k_23   = 0.0781
    k_24   = 0.0164
    k_31   = 18.*(140.**(-beta_3))
    k_34   = 0.714
    k_42   = 0.00189
    k_43   = 0.00114
    k_51   = 0.0862
    k_56   = 0.0862
    k_61   = 0.0333

    def __init__(self, M1, M2, M3, M4, M5, M6, M7, G, emis_func,
                 **kwargs):

        ## Initial Conditions
        self.M1 = M1
        self.M2 = M2
        self.M3 = M3
        self.M4 = M4
        self.M5 = M5
        self.M6 = M6
        self.M7 = M7
        self.G  = G

        self.emis_func = emis_func

        self.y0 = [ M1, M2, M3, M4, M5, M6, M7, G ]

        for key, value in kwargs.iteritems():
            self.__dict__[key] = value

    def __call__(self, y, t):
        """ Alias to call the ``rhs'' method.

        """
        return self.rhs(y, t)

    def rhs(self, y, t):
        """ Evaluate the RHS of the model governing equations

        """
    
        M1, M2, M3, M4, M5, M6, M7, G = y[:]
        F_r, F_d, F_f = self.emis_func(t)
        
        dM1_dt = -1.*(self.k_12 + self.k_13)*M1 \
               - self.k_15*G*(M1 - self.gamma)/(M1 + self.Gamma) \
               + self.k_21*(M2**self.beta_2) + self.k_31*(M3**self.beta_3) \
               + self.k_51*M5 \
               + self.k_61*M6 + F_f + F_d - F_r
        dM2_dt = self.k_12*M1 \
               - (self.k_23 + self.k_24)*M2 - self.k_21*(M2**self.beta_2) \
               + self.k_42*M4
        dM3_dt = self.k_13*M1 + self.k_23*M2 - self.k_34*M3 \
               - self.k_31*(M3**self.beta_3) \
               + self.k_43*M4
        dM4_dt = self.k_24*M2 + self.k_34*M3 - (self.k_42 + self.k_43)*M4
        dM5_dt = self.k_15*G*(M1 - self.gamma)/(M1 + self.Gamma) \
               - (self.k_51 + self.k_56)*M5 \
               - F_d + F_r
        dM6_dt = self.k_56*M5 - self.k_61*M6
        dM7_dt = -1.*F_f
        dG_dt  = -1.*(self.a_d*F_d - self.a_r*F_r)/self.M5 # note that self.M5 is
                                                           # frozen to the initial
                                                           # M5
        
        dy_dt = [ dM1_dt, dM2_dt, dM3_dt, dM4_dt, dM5_dt, dM6_dt,
                  dM7_dt, dG_dt ]
        
        return np.array(dy_dt)

    def _output(self, t, u, t_offset):
        """ Process the model output into a descriptive DataFrame 
        """
        M_cols = ["M%1d" % i for i in range(1, 8)]
        tt = pd.Series(t, name="time") + t_offset
        yy = pd.DataFrame(u, 
                          columns=M_cols + ["G", ],
                          index=tt)

        ppm_fac = 1./2.13 #: ppmv / PgC

        yy['atm_ppm'] = yy.M1*ppm_fac

        yy_emis = pd.DataFrame([self.emis_func(ti) for ti in tt],
                               columns=["F_r", "F_d", "F_f"],
                               index=tt)

        yy = pd.concat([yy, yy_emis], axis=1)

        return yy

    def integrate(self, integrator, t_end, dt=1./365.,
                  t_offset=1850., ):
        """ Interface for integrating the model using custom
        integration schemes.
        """

        t = np.arange(0., t_end+dt, dt)
        y = integrator(self, self.y0, t)

        return self._output(t, y, t_offset)

    def integrate_odespy(self, integrator, t_end, dt=1./365.,
                         t_offset=1850.):
        """ Interface for running the model using integration
        schemes from odespy package
        """

        t = np.arange(0., t_end+dt, dt)
        solver = integrator(self)
        solver.set_initial_condition(self.y0)
        y, t = solver.solve(t)

        return self._output(t, y, t_offset)


## Boundary conditions
def emissions_22p6(t, switch=100.):
    """ Return 3 values, corresponding to fluxes F_r, F_d, and F_f 
    evaluated at a given time 't' 
    
    Fluxes are in Pg(C)/yr; time 't' is in years since
    a particular baseline (1850 in this case)
    """
    F_r = 0.
    F_d = 0.3 + 0.01*t
    F_f = 0.014*t if t <= switch else 1.4 + (4.6/40.)*(t - 100.)
    
    return F_r, F_d, F_f

def emissions_22p7(t):
    """ Similar to above, except for the emissions scenario in problem
    22.7, which has a baseline reference year of 2010
    """

    F_r = 0.
    F_d = 1.8 + 0.01*t
    F_f = 6. + (4.6/40)*t

    return F_r, F_d, F_f
    