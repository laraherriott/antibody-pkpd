import math
import numpy as np
import random

from parameters.parameters import Parameters
from parameters.day_parameters import DayParameters
from parameters.noIIV import NoIIVParameters
from patient import Patient, FixedPatient
from dosage import DoseFn

class LecanemabModel:
    def __init__(self, PD_model: str, dose, dose_interval, final_dose, dose_change = None, median_patient = False, param_type = 'day'):
        if param_type == 'day':
            self.PK_params = DayParameters()
        elif param_type == 'hour':
            self.PK_params = Parameters()
        elif param_type == 'noiiv':
            self.PK_params = NoIIVParameters()

        if median_patient:
            self.patient = FixedPatient()
        else:
            self.patient = Patient()
        self.dose = dose

        if dose_change is not None:
            self.next_dose = dose_change[0]
            next_interval = dose_change[1]
            self.switch_time = dose_change[2]
        else:
            self.next_dose = self.dose
            next_interval = dose_interval
            self.switch_time = final_dose + 1

        self.dose_list = []
        i = 0
        while i <= final_dose:
            self.dose_list.append(i)
            if i < self.switch_time:
                i += dose_interval
            else:
                i += next_interval
        self.PD = PD_model

    def PK_model(self, albumin, sex, weight, ada, race, process):
        # the random variance will be added per model run which means per individual
        self.CL = ((self.PK_params.clearance * ((weight/73.7)**self.PK_params.weight_clearance) *
                   ((albumin/42.9)**self.PK_params.albumin_clearance) * (self.PK_params.sex_clearance**sex) *
                   (self.PK_params.ADA_clearance**ada)) *
                   math.exp(random.normalvariate(0, math.sqrt(self.PK_params.random_clearance)))) 

        self.V1 = ((self.PK_params.volume1 * ((weight/73.7)**self.PK_params.weight_volume1) *
                   (self.PK_params.sex_volume1**sex)) *
                   math.exp(random.normalvariate(0, math.sqrt(self.PK_params.random_volume1))))
        
        self.V2 = ((self.PK_params.volume2 * (self.PK_params.japanese_volume2**race)) *
                   math.exp(random.normalvariate(0, math.sqrt(self.PK_params.random_volume1))))
        
        self.Q = self.PK_params.inter_compartment

        if process == 0:
            self.F = 1
        else:
            self.F = ((self.PK_params.processB) *
                      math.exp(random.normalvariate(0, math.sqrt(self.PK_params.random_F))))

    def SUVr_fixed(self, apoe, age):
        self.baseline_SUVr = (self.PK_params.SUVr_baseline * (self.PK_params.APOE_SUVr_baseline ** apoe) *
                              math.exp(random.normalvariate(0, self.PK_params.random_SUVr_baseline)))
        # note baseline acts as y0 for SUVr
        self.Emax = (self.PK_params.SUVr_Emax * ((age/72) ** self.PK_params.age_Emax) *
                     math.exp(random.normalvariate(0, self.PK_params.random_Emax)))

        #self.SUVr_Kin = self.PK_params.SUVr_Kin / 365 
        self.SUVr_Kin = self.PK_params.SUVr_Kin / 8760 # this is number for per hour model
        self.SUVr_EC50 = self.PK_params.SUVr_EC50
        self.SUVr_Kout = self.SUVr_Kin/self.baseline_SUVr

    def Abeta_fixed(self):
        self.baseline_Abeta = (self.PK_params.AB42_40_baseline * 
                               math.exp(random.normalvariate(0, self.PK_params.random_AB42_40_baseline)))
        self.Abeta_Kout = self.PK_params.AB42_40_Kout / 730
        self.Abeta_slope = (self.PK_params.AB42_40_slope * 
                            math.exp(random.normalvariate(0, self.PK_params.random_AB42_40_slope)))
        self.Abeta_Kin = self.Abeta_Kout * self.baseline_Abeta

    def tau_fixed(self, weight):
        self.baseline_tau = (self.PK_params.tau_baseline * ((weight/72.2)**self.PK_params.weight_tau_baseline) *
                             math.exp(random.normalvariate(0, self.PK_params.random_tau_baseline)))
        self.tau_Kout = self.PK_params.tau_Kout / 730
        self.tau_slope = (self.PK_params.tau_slope * 
                          math.exp(random.normalvariate(0, self.PK_params.random_tau_slope)))
        self.tau_Kin = self.tau_Kout * self.baseline_tau
    
    def SUVr_PD(self, t, y): # y = [Lcent, Lper, SUVr]
        C1L = y[0]/self.V1
        dLcent_dt = (self.Q/self.V2)*y[1] - ((self.Q/self.V1) + (self.CL/self.V1))*y[0] + self.dosefn(self.dose, self.dose_list, self.next_dose, self.switch_time, self.patient, t)
        dLper_dt = (self.Q/self.V1)*y[0] - (self.Q/self.V2)*y[1]
        dSUVr_dt = self.SUVr_Kin - (y[2] * self.SUVr_Kout * (1 + ((self.Emax * C1L)/(self.SUVr_EC50 + C1L))))

        dYdt = [dLcent_dt, dLper_dt, dSUVr_dt]
    
        return dYdt
    
    def AB_PD(self, t, y):
        C1L = y[0]/self.V1
        dLcent_dt = (self.Q/self.V2)*y[1] - ((self.Q/self.V1) + (self.CL/self.V1))*y[0] + self.dosefn(self.dose, self.dose_list, self.next_dose, self.switch_time, self.patient, t)
        dLper_dt = (self.Q/self.V1)*y[0] - (self.Q/self.V2)*y[1]
        dAB_dt = self.Abeta_Kin*(1+self.Abeta_slope*C1L) - self.Abeta_Kout*y[2]

        dYdt = [dLcent_dt, dLper_dt, dAB_dt]

        return dYdt
    
    def tau_PD(self, t, y):
        C1L = y[0]/self.V1
        dLcent_dt = (self.Q/self.V2)*y[1] - ((self.Q/self.V1) + (self.CL/self.V1))*y[0] + self.dosefn(self.dose, self.dose_list, self.next_dose, self.switch_time, self.patient, t)
        dLper_dt = (self.Q/self.V1)*y[0] - (self.Q/self.V2)*y[1]
        dtau_dt = self.tau_Kin*(1-self.tau_slope*C1L) - self.tau_Kout*y[2]

        dYdt = [dLcent_dt, dLper_dt, dtau_dt]

        return dYdt
    
    def combined_PD(self, t, y):
        dose = DoseFn(self.dose, self.interval, self.final_dose, self.patient)
        C1L = y[0]/self.V1
        dLcent_dt = (self.Q/self.V2)*y[1] - ((self.Q/self.V1) + (self.CL/self.V1))*y[0] + dose.eval_at(t)
        dLper_dt = (self.Q/self.V1)*y[0] - (self.Q/self.V2)*y[1]
        dSUVr_dt = self.SUVr_Kin - (y[2] * self.SUVr_Kout * (1 + ((self.Emax * C1L)/(self.SUVr_EC50 + C1L))))
        dAB_dt = self.Abeta_Kin*(1+self.Abeta_slope+C1L) - self.Abeta_Kout*y[3]
        dtau_dt = self.tau_Kin*(1+self.tau_slope+C1L) - self.tau_Kout*y[4]

        dYdt = [dLcent_dt, dLper_dt, dSUVr_dt, dAB_dt, dtau_dt]

        return dYdt
    
    def dosefn(self, dose, dose_list, new_dose, switch_time, patient, t):
        infusion_1 = dose * patient.weight_SUVr
        infusion_2 = new_dose * patient.weight_SUVr
        f = 0.5
        delta = 0.1

        sol = 0
        for n in dose_list:
            if t >=(n-0.5) and t<=(n+1.5):
                if n < switch_time:
                    sol =  ((infusion_1/2)/math.atan(1/delta))*(math.atan(math.sin(2*math.pi*(t-n-0.5)*f)/delta)) + infusion_1/2
                else:
                    sol =  ((infusion_2/2)/math.atan(1/delta))*(math.atan(math.sin(2*math.pi*(t-n-0.5)*f)/delta)) + infusion_2/2
        return sol

    def __call__(self):
        self.PK_model(self.patient.albumin, self.patient.sex, self.patient.weight_SUVr,
                      self.patient.ADA, self.patient.race, self.patient.process)
        if self.PD == 'suvr':
            self.SUVr_fixed(self.patient.APOE, self.patient.age_SUVr)
        elif self.PD == 'ab':
            self.Abeta_fixed()
        elif self.PD == 'tau':
            self.tau_fixed(self.patient.weight_SUVr) # should I change this so the individuals being simulated for each
            # biomarker are representative of the population that was used to generate the initial model? i.e. from the clinical trials?
        elif self.PD == 'all':
            self.SUVr_fixed(self.patient.APOE, self.patient.age_SUVr)
            self.Abeta_fixed()
            self.tau_fixed(self.patient.weight_SUVr)

class LecanemabModel2:
    def __init__(self, PD_model: str, dose, dose_interval, final_dose, dose_change = None, median_patient = False, param_type = 'day'):
        if param_type == 'day':
            self.PK_params = DayParameters()
        elif param_type == 'hour':
            self.PK_params = Parameters()
        elif param_type == 'noiiv':
            self.PK_params = NoIIVParameters()

        if median_patient:
            self.patient = FixedPatient()
        else:
            self.patient = Patient()
        self.dose = dose

        if dose_change is not None:
            self.next_dose = dose_change[0]
            next_interval = dose_change[1]
            self.switch_time = dose_change[2]
        else:
            self.next_dose = self.dose
            next_interval = dose_interval
            self.switch_time = final_dose + 1

        self.dose_list = []
        i = 0
        while i <= final_dose:
            self.dose_list.append(i)
            if i < self.switch_time:
                i += dose_interval
            else:
                i += next_interval
        self.PD = PD_model

    def PK_model(self, albumin, sex, weight, ada, race, process):
        # the random variance will be added per model run which means per individual
        self.CL = ((self.PK_params.clearance * ((weight/73.7)**self.PK_params.weight_clearance) *
                   ((albumin/42.9)**self.PK_params.albumin_clearance) * (self.PK_params.sex_clearance**sex) *
                   (self.PK_params.ADA_clearance**ada)) *
                   math.exp(random.normalvariate(0, math.sqrt(self.PK_params.random_clearance)))) 

        self.V1 = ((self.PK_params.volume1 * ((weight/73.7)**self.PK_params.weight_volume1) *
                   (self.PK_params.sex_volume1**sex)) *
                   math.exp(random.normalvariate(0, math.sqrt(self.PK_params.random_volume1))))
        
        self.V2 = ((self.PK_params.volume2 * (self.PK_params.japanese_volume2**race)) *
                   math.exp(random.normalvariate(0, math.sqrt(self.PK_params.random_volume1))))
        
        self.Q = self.PK_params.inter_compartment

        if process == 0:
            self.F = 1
        else:
            self.F = ((self.PK_params.processB) *
                      math.exp(random.normalvariate(0, math.sqrt(self.PK_params.random_F))))

    def SUVr_fixed(self, apoe, age):
        self.baseline_SUVr = (self.PK_params.SUVr_baseline * (self.PK_params.APOE_SUVr_baseline ** apoe) *
                              math.exp(random.normalvariate(0, self.PK_params.random_SUVr_baseline)))
        # note baseline acts as y0 for SUVr
        self.Emax = (self.PK_params.SUVr_Emax * ((age/72) ** self.PK_params.age_Emax) *
                     math.exp(random.normalvariate(0, self.PK_params.random_Emax)))

        #self.SUVr_Kin = self.PK_params.SUVr_Kin / 365 
        self.SUVr_Kin = self.PK_params.SUVr_Kin / 8760 # this is number for per hour model
        self.SUVr_EC50 = self.PK_params.SUVr_EC50
        self.SUVr_Kout = self.SUVr_Kin/self.baseline_SUVr

    def Abeta_fixed(self):
        self.baseline_Abeta = (self.PK_params.AB42_40_baseline * 
                               math.exp(random.normalvariate(0, self.PK_params.random_AB42_40_baseline)))
        self.Abeta_Kout = self.PK_params.AB42_40_Kout / 8760
        self.Abeta_slope = (self.PK_params.AB42_40_slope * 
                            math.exp(random.normalvariate(0, self.PK_params.random_AB42_40_slope)))
        self.Abeta_Kin = self.Abeta_Kout * self.baseline_Abeta

    def tau_fixed(self, weight):
        self.baseline_tau = (self.PK_params.tau_baseline * ((weight/72.2)**self.PK_params.weight_tau_baseline) *
                             math.exp(random.normalvariate(0, self.PK_params.random_tau_baseline)))
        self.tau_Kout = self.PK_params.tau_Kout / 8760
        self.tau_slope = (self.PK_params.tau_slope * 
                          math.exp(random.normalvariate(0, self.PK_params.random_tau_slope)))
        self.tau_Kin = self.tau_Kout * self.baseline_tau
    
    def SUVr_PD(self, t, y): # y = [Lcent, Lper, SUVr]
        if t < self.switch_time:
            C1L = y[0]/self.V1
        else:
            C1L = ((self.next_dose/1e9)*147181.62*3.22*1000)/self.V1
        dLcent_dt = (self.Q/self.V2)*y[1] - ((self.Q/self.V1) + (self.CL/self.V1))*y[0] + self.dosefn(self.dose, self.dose_list, self.next_dose, self.switch_time, self.patient, t)
        dLper_dt = (self.Q/self.V1)*y[0] - (self.Q/self.V2)*y[1]
        dSUVr_dt = self.SUVr_Kin - (y[2] * self.SUVr_Kout * (1 + ((self.Emax * C1L)/(self.SUVr_EC50 + C1L))))

        dYdt = [dLcent_dt, dLper_dt, dSUVr_dt]
    
        return dYdt
    
    def AB_PD(self, t, y):
        if t < self.switch_time:
            C1L = y[0]/self.V1
        else:
            C1L = ((self.next_dose/1e9)*147181.62*3.22*1000)/self.V1
        dLcent_dt = (self.Q/self.V2)*y[1] - ((self.Q/self.V1) + (self.CL/self.V1))*y[0] + self.dosefn(self.dose, self.dose_list, self.next_dose, self.switch_time, self.patient, t)
        dLper_dt = (self.Q/self.V1)*y[0] - (self.Q/self.V2)*y[1]
        dAB_dt = self.Abeta_Kin*(1+self.Abeta_slope*C1L) - self.Abeta_Kout*y[2]

        dYdt = [dLcent_dt, dLper_dt, dAB_dt]

        return dYdt
    
    def tau_PD(self, t, y):
        if t < self.switch_time:
            C1L = y[0]/self.V1
        else:
            C1L = ((self.next_dose/1e9)*147181.62*3.22*1000)/self.V1
        dLcent_dt = (self.Q/self.V2)*y[1] - ((self.Q/self.V1) + (self.CL/self.V1))*y[0] + self.dosefn(self.dose, self.dose_list, self.next_dose, self.switch_time, self.patient, t)
        dLper_dt = (self.Q/self.V1)*y[0] - (self.Q/self.V2)*y[1]
        dtau_dt = self.tau_Kin*(1-self.tau_slope*C1L) - self.tau_Kout*y[2]

        dYdt = [dLcent_dt, dLper_dt, dtau_dt]

        return dYdt
    
    def combined_PD(self, t, y):
        dose = DoseFn(self.dose, self.interval, self.final_dose, self.patient)
        C1L = y[0]/self.V1
        dLcent_dt = (self.Q/self.V2)*y[1] - ((self.Q/self.V1) + (self.CL/self.V1))*y[0] + dose.eval_at(t)
        dLper_dt = (self.Q/self.V1)*y[0] - (self.Q/self.V2)*y[1]
        dSUVr_dt = self.SUVr_Kin - (y[2] * self.SUVr_Kout * (1 + ((self.Emax * C1L)/(self.SUVr_EC50 + C1L))))
        dAB_dt = self.Abeta_Kin*(1+self.Abeta_slope+C1L) - self.Abeta_Kout*y[3]
        dtau_dt = self.tau_Kin*(1+self.tau_slope+C1L) - self.tau_Kout*y[4]

        dYdt = [dLcent_dt, dLper_dt, dSUVr_dt, dAB_dt, dtau_dt]

        return dYdt
    
    def dosefn(self, dose, dose_list, new_dose, switch_time, patient, t):
        infusion_1 = dose * patient.weight_SUVr
        infusion_2 = new_dose * patient.weight_SUVr
        f = 0.5
        delta = 0.1

        sol = 0
        if t < switch_time:
            for n in dose_list:
                if t >=(n-0.5) and t<=(n+1.5):
                    sol =  ((infusion_1/2)/math.atan(1/delta))*(math.atan(math.sin(2*math.pi*(t-n-0.5)*f)/delta)) + infusion_1/2
        # else:
        #     sol= (new_dose/1e9)*147181.62*3.22*1000 # which is actually just a continuous antibody concentration from this point on.
        return sol

    def __call__(self):
        self.PK_model(self.patient.albumin, self.patient.sex, self.patient.weight_SUVr,
                      self.patient.ADA, self.patient.race, self.patient.process)
        if self.PD == 'suvr':
            self.SUVr_fixed(self.patient.APOE, self.patient.age_SUVr)
        elif self.PD == 'ab':
            self.Abeta_fixed()
        elif self.PD == 'tau':
            self.tau_fixed(self.patient.weight_SUVr) # should I change this so the individuals being simulated for each
            # biomarker are representative of the population that was used to generate the initial model? i.e. from the clinical trials?
        elif self.PD == 'all':
            self.SUVr_fixed(self.patient.APOE, self.patient.age_SUVr)
            self.Abeta_fixed()
            self.tau_fixed(self.patient.weight_SUVr)