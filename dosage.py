#
# Class to outline the doasage regimen
#

import numpy as np


class GaussConvFn:
    """
    This class represents the convolution of a delta function with Gaussian function
    By Gaussian function we mean the pdf of Gaussian distribution.
    The purpose of this function is build up a smooth function which is close to a delta function,
    and enable the numerical solution of ODE using Runge-Kutta method.
    Formally, Let us denote \rho(t)=A * \delta(t-t_center), where A is the magnitude
    \omega(t) = \int \frac{1}{\sqrt{2\pi}\sigma}e^{-\frac{t^2}{2\sigma^2}},
    then this class just represent \rho * \omega (t), where * is the convolution.
    See the definition in https://en.wikipedia.org/wiki/Convolution for more information.
    """  
    def __init__(self, center: float, magnitude: float, sigma=0.02):
        """
        params:
        center: t_center
        magnitude: magnitude in the delta function
        sigma: sigma in the pdf of Gaussian distribution.
        Note that sigma could not be too small, otherwise ODE solving might fail using Runge-Kutta method.
        """
        self.center = center
        self.magnitude = magnitude
        self.sigma = sigma

    def eval_at(self, x: float) -> float:
        '''
        Return the value of this function at x.
        '''
        return self.magnitude / (2 * np.pi) ** 0.5 / self.sigma * np.exp(-(x - self.center)**2 / (2 * self.sigma ** 2))

class DeltaFn:
    def __init__(self, center: float, magnitude: float):
        self.center = center
        self.magnitude = magnitude

    def eval_at(self, t):
        if t == self.center:
            return self.magnitude
        else:
            return 0
        
class StepFn:
    def __init__(self, center: float, magnitude: float, time_list, step=1):
        self.center = center
        self.magnitude = magnitude/step
        self.step = step
        self.time_list = time_list

    def eval_at(self, t):
        t_range = [self.center, self.center + self.step]
        if t >= t_range[0] and t <= t_range[1]:
            return self.magnitude
            
        else:
            return 0

class DoseFn:
    """
    This class represents the dose function in the ODE.
    The dose function DOSE(t) should be a linear combination of several pseudo delta function(see GaussConvFn),
    plus a constant value.
    It represents consist of instantaneous doses of X ng of the drug at one or more time points,
    or a steady application of X ng per hour over a given time period, or some combination.
    """
    
    def __init__(self, dose, interval, final_dose, patient):
        '''
        params:
        dose: size of instantaneous doses, in mg/kg
        interval: time between doses, in time units of choice
        final_dose: time of final dose
        '''

        # note to self it seems still to be doing this at every time step

        self.dose = dose * patient.weight_SUVr
        self.dose_times = []
        t = 1
        self.dose_times.append(t)
        while t <= (final_dose):
            next_time = t + (interval)
            self.dose_times.append(next_time)
            t = next_time

        self.deltainput = []
        for i in range(len(self.dose_times)):
            # write a line here to define and call the correct dosage function
            #self.deltainput.append(GaussConvFn(self.dose_times[i], self.dose))
            self.deltainput.append(DeltaFn(self.dose_times[i], self.dose))
            #self.deltainput.append(StepFn(self.dose_times[i], self.dose))

            # note that in the way the solver is currently set up (fixed time step) the delta and gausconv functions are equivalent.
            # problem with the step function arose when were able to sample from same dosage peak multiple times when allowing the solver to choose its own timestep.

    def get_profile(self, n):
        # for plotting the dosage function
        result = []
        for t in range(n):
            result_t = 0
            if t in self.dose_times:
                for i in range(len(self.deltainput)):
                    if t == self.dose_times[i]:
                        result_t += self.deltainput[i].eval_at(t)
            result.append(result_t)
        return result

    def eval_at(self, t):
        '''
        Return the dose function value at time t
        '''
        
        result = 0
        if t in self.dose_times:
            for i in range(len(self.deltainput)):
                if t == self.dose_times[i]:
                    result += self.deltainput[i].eval_at(t)

        return result