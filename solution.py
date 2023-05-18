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
            self.y0 = [0, 0, self.model.baseline_SUVr, 0]
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
        elif self.model.PD == 'suvr_ab':
            self.name = self.model.model
            self.y0 = [39, 0, 0, 0, 1] #, 1.38]
        elif self.model.PD == 'adu_path':
            self.name = self.model.pathway
            self.y0 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        elif self.model.PD == 'brain':
            self.name = self.model.brain
            self.y0 = [0.2, 370, 5500, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        elif self.model.PD == 'brain2':
            self.name = self.model.brain
            self.y0 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        elif self.model.PD == 'half':
            self.name = self.model.half_life
            self.y0 = [0]

    def solve(self):
        solution = scipy.integrate.solve_ivp(fun=lambda t, y: self.name(t, y),
                                             t_span=[self.t_eval[0], self.t_eval[-1]],
                                             y0=self.y0,
                                             t_eval=self.t_eval)

        #solution = self.euler([self.t_start, self.t_end], self.y0, self.name)
        
        return solution

    # forward euler method
    def euler(self, t_span, y0, equations): # need 12
        # other parameters
        h = 0.01 # step size
        t0 = t_span[0] # initial time
        tstop = t_span[1] # final time
        y_1 = [0]*int(((tstop-t0)/h)+1)
        y_2 = [0]*int(((tstop-t0)/h)+1)
        y_3 = [0]*int(((tstop-t0)/h)+1)
        y_4 = [0]*int(((tstop-t0)/h)+1)
        y_5 = [0]*int(((tstop-t0)/h)+1) # pre-allocate memory for output
        y_6 = [0]*int(((tstop-t0)/h)+1)
        y_7 = [0]*int(((tstop-t0)/h)+1)
        y_8 = [0]*int(((tstop-t0)/h)+1)
        y_9 = [0]*int(((tstop-t0)/h)+1)
        y_10 = [0]*int(((tstop-t0)/h)+1)
        y_11 = [0]*int(((tstop-t0)/h)+1)
        #y_12 = [0]*int(((tstop-t0)/h)+1)

        y_1[0] = y0[0]
        y_2[0] = y0[1]
        y_3[0] = y0[2]
        y_4[0] = y0[3]
        y_5[0] = y0[4]
        y_6[0] = y0[5]
        y_7[0] = y0[6]
        y_8[0] = y0[7]
        y_9[0] = y0[8]
        y_10[0] = y0[9]
        y_11[0] = y0[10]
        #y_12[0] = y0[11]

        
        t_list = []
        t_list.append(t0)

        # for i in range(3):
        t = t0 # initialise time variable
        n = 1

        while t < tstop: # condition such that stops when tn > tstop
            
            t_list.append(t)
            y_1[n] = y_1[n-1] + h*(equations(t, [y_1[n-1], y_2[n-1], y_3[n-1], y_4[n-1], y_5[n-1], y_6[n-1], y_7[n-1], y_8[n-1], y_9[n-1], y_10[n-1], y_11[n-1]])[0])
            y_2[n] = y_2[n-1] + h*(equations(t, [y_1[n-1], y_2[n-1], y_3[n-1], y_4[n-1], y_5[n-1], y_6[n-1], y_7[n-1], y_8[n-1], y_9[n-1], y_10[n-1], y_11[n-1]])[1])
            y_3[n] = y_3[n-1] + h*(equations(t, [y_1[n-1], y_2[n-1], y_3[n-1], y_4[n-1], y_5[n-1], y_6[n-1], y_7[n-1], y_8[n-1], y_9[n-1], y_10[n-1], y_11[n-1]])[2]) 
            y_4[n] = y_4[n-1] + h*(equations(t, [y_1[n-1], y_2[n-1], y_3[n-1], y_4[n-1], y_5[n-1], y_6[n-1], y_7[n-1], y_8[n-1], y_9[n-1], y_10[n-1], y_11[n-1]])[3]) 
            y_5[n] = y_5[n-1] + h*(equations(t, [y_1[n-1], y_2[n-1], y_3[n-1], y_4[n-1], y_5[n-1], y_6[n-1], y_7[n-1], y_8[n-1], y_9[n-1], y_10[n-1], y_11[n-1]])[4]) # forward euler formula
            y_6[n] = y_6[n-1] + h*(equations(t, [y_1[n-1], y_2[n-1], y_3[n-1], y_4[n-1], y_5[n-1], y_6[n-1], y_7[n-1], y_8[n-1], y_9[n-1], y_10[n-1], y_11[n-1]])[5])
            y_7[n] = y_7[n-1] + h*(equations(t, [y_1[n-1], y_2[n-1], y_3[n-1], y_4[n-1], y_5[n-1], y_6[n-1], y_7[n-1], y_8[n-1], y_9[n-1], y_10[n-1], y_11[n-1]])[6])
            y_8[n] = y_8[n-1] + h*(equations(t, [y_1[n-1], y_2[n-1], y_3[n-1], y_4[n-1], y_5[n-1], y_6[n-1], y_7[n-1], y_8[n-1], y_9[n-1], y_10[n-1], y_11[n-1]])[7])
            y_9[n] = y_9[n-1] + h*(equations(t, [y_1[n-1], y_2[n-1], y_3[n-1], y_4[n-1], y_5[n-1], y_6[n-1], y_7[n-1], y_8[n-1], y_9[n-1], y_10[n-1], y_11[n-1]])[8])
            y_10[n] = y_10[n-1] + h*(equations(t, [y_1[n-1], y_2[n-1], y_3[n-1], y_4[n-1], y_5[n-1], y_6[n-1], y_7[n-1], y_8[n-1], y_9[n-1], y_10[n-1], y_11[n-1]])[9])
            y_11[n] = y_11[n-1] + h*(equations(t, [y_1[n-1], y_2[n-1], y_3[n-1], y_4[n-1], y_5[n-1], y_6[n-1], y_7[n-1], y_8[n-1], y_9[n-1], y_10[n-1], y_11[n-1]])[10])
            #y_12[n] = y_12[n-1] + h*(equations(t, [y_1[n-1], y_2[n-1], y_3[n-1], y_4[n-1], y_5[n-1], y_6[n-1], y_7[n-1], y_8[n-1], y_9[n-1], y_10[n-1], y_11[n-1], y_12[n-1]])[11])
            t += h
            n += 1
            #t = round(t, 1)
        return [y_1, y_2, y_3, y_4, y_5, y_6, y_7, y_8, y_9, y_10, y_11], t_list