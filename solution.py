#
# Class to solve the differential equations
#

import numpy as np
import scipy.integrate

class Solution:
    def __init__(self, model, t_0, t_end, step_size):
        self.t_eval = list(range(t_0, t_end, step_size))
        self.model = model
        if self.model.PD == 'suvr':
            self.name = self.model.SUVr_PD
            self.y0 = [0, 0, self.model.baseline_SUVr]
        elif self.model.PD == 'ab':
            self.name = self.model.AB_PD
            self.y0 = [0, 0, self.model.baseline_Abeta]
        elif self.model.PD == 'tau':
            self.name = self.model.tau_PD
            self.y0 = [0, 0, self.model.baseline_tau]
        elif self.model.PD == 'all':
            self.name = self.model.combined_PD
            self.y0 = [0, 0, self.model.baseline_SUVr, self.model.baseline_Abeta, self.model.baseline_tau]

    def solve(self):
        solution = scipy.integrate.solve_ivp(fun=lambda t, y: self.name(t, y),
                                             t_span=[self.t_eval[0], self.t_eval[-1]],
                                             y0=self.y0,
                                             t_eval=self.t_eval,
                                             max_step = 1,
                                             first_step = 1,
                                             rtol = 1e90,
                                             atol = 1e90)
        
        return solution