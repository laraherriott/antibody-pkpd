#
# Class to solve the differential equations
#

import numpy as np
import scipy.integrate

class Solution:
    def __init__(self, model, t_0, t_end, step_size):
        self.t_start = t_0
        self.t_end = t_end 
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
        elif self.model.PD == 'constant':
            self.name = self.model.constant
            self.y0 = [0, 0]
        elif self.model.PD == 'sine':
            self.name = self.model.sine
            self.y0 = [0, 0]
        elif self.model.PD == 'narrow_sine':
            self.name = self.model.narrow_sine
            self.y0 = [0, 0]
        elif self.model.PD == 'narrower_sine':
            self.name = self.model.narrower_sine
            self.y0 = [0, 0]
        elif self.model.PD == 'square':
            self.name = self.model.smooth_square
            self.y0 = [0, 0]

    def solve(self):
        solution = scipy.integrate.solve_ivp(fun=lambda t, y: self.name(t, y),
                                             t_span=[self.t_eval[0], self.t_eval[-1]],
                                             y0=self.y0,
                                             t_eval=self.t_eval,
                                             max_step = 0.1)

        #solution = self.euler([self.t_start, self.t_end], self.y0, self.name)
        
        return solution

    # forward euler method
    def euler(self, t_span, y0, equations):
        # other parameters
        h = 0.01 # step size
        t0 = t_span[0] # initial time
        tstop = t_span[1] # final time
        y_1 = [0]*int(((tstop-t0)/h)+1)
        y_2 = [0]*int(((tstop-t0)/h)+1)
        y_3 = [0]*int(((tstop-t0)/h)+1) # pre-allocate memory for output

        y_1[0] = y0[0]
        y_2[0] = y0[1]
        y_3[0] = y0[2]

        
        t_list = []
        t_list.append(t0)

        # for i in range(3):
        t = t0 # initialise time variable
        n = 1

        while t < tstop: # condition such that stops when tn > tstop
            
            t_list.append(t)
            y_1[n] = y_1[n-1] + h*(equations(t, [y_1[n-1], y_2[n-1], y_3[n-1]])[0])
            y_2[n] = y_2[n-1] + h*(equations(t, [y_1[n-1], y_2[n-1], y_3[n-1]])[1])
            y_3[n] = y_3[n-1] + h*(equations(t, [y_1[n-1], y_2[n-1], y_3[n-1]])[2]) # forward euler formula
            t += h
            n += 1
            #t = round(t, 1)
        return [y_1, y_2, y_3], t_list