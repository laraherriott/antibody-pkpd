# for running model with antibody profiles

import math
import numpy as np
import pandas as pd
import random
from  scipy.interpolate import CubicSpline

from parameters.parameters import Parameters
from parameters.day_parameters import DayParameters
from parameters.noIIV import NoIIVParameters
from patient import Patient, FixedPatient
from dosage import DoseFn

class AbModel:
    def __init__(self):
        self.PK_params = Parameters()
        self.patient = FixedPatient()

        self.PD = 'suvr_ab'

        data = pd.read_csv('data/vaccinity_profile_UB311Q3M.csv')
        data['conc'] = data['titre'].apply(lambda x: x*(120/2.6))
        data['hours'] = data['time'].apply(lambda x: x*7*24)
        self.spline = CubicSpline(data['time'], data['conc'])

    def SUVr_fixed(self, apoe, age):
        self.baseline_SUVr = (self.PK_params.SUVr_baseline * (self.PK_params.APOE_SUVr_baseline ** apoe) *
                              math.exp(random.normalvariate(0, self.PK_params.random_SUVr_baseline)))
        # note baseline acts as y0 for SUVr
        self.Emax = (self.PK_params.SUVr_Emax * ((age/72) ** self.PK_params.age_Emax) *
                     math.exp(random.normalvariate(0, self.PK_params.random_Emax)))

        #self.SUVr_Kin = self.PK_params.SUVr_Kin / 365 
        self.SUVr_Kin = self.PK_params.SUVr_Kin / 8760 # this is number for per hour model
        self.SUVr_Kout = self.SUVr_Kin/self.baseline_SUVr

        self.SUVr_EC50_1 = self.PK_params.SUVr_EC50 * 1
        self.SUVr_EC50_2 = self.PK_params.SUVr_EC50 * 0.5

        Kd_P = 1.79 # nM
        Kd_Fc = 0.12 # nM

        self.K_onP = 0.001*60*60 # nM/h
        self.K_onFc = 0.001*60*60 # nM/h
        self.K_offP = Kd_P * self.K_onP
        self.K_offFc = Kd_Fc * self.K_onFc
        self.clear = 0.000000008*60*60
        self.SCALE = 200
        self.C1 = 1
        self.C3 = 1.3
        self.C4 = 3.5

    def model(self, t, y): # y = [plaque, Ab, plaque-Ab, plaque-Ab-Fc, Fc]
        dPlaque = self.SUVr_Kin - self.clear*y[0] - self.K_onP*y[0]*y[1] + self.K_offP*y[2]
        dAntibody = self.antibody_fn(t) - self.K_onP*y[0]*y[1] + self.K_offP*y[2] 
        dAntibodyPlaque = self.K_onP*y[0]*y[1] - self.K_offFc*y[2] - self.K_onFc*y[2]*y[4] + self.K_offFc*y[3]
        dAntibodyPlaqueFc = self.K_onFc*y[2]*y[4] - self.K_offFc*y[3] - self.SCALE*self.SUVr_Kout*y[3]
        dFc = - self.K_onFc*y[2]*y[4] + self.K_offFc*y[3] + self.SCALE*self.SUVr_Kout*y[3]

        #dSUVr = 1 + self.C1*(y[0]**self.C3)/((y[0]**self.C3)+(self.C4**self.C3))

        # SCALE = how much more the presence of the antibody increases the rate of clearance
        # Kout = 5.78 per year, 0.00000018 per s, vs k_clear = 0.000037 per s (so Ab ~ 200x greater clearance)

        # y0 values: Plaque = 39 nM OR baseline SUVr?; Ab = 0; AbP = 0; AbPFc = 0, Fc = 1

        dYdt = [dPlaque, dAntibody, dAntibodyPlaque, dAntibodyPlaqueFc, dFc]#, dSUVr]
        return dYdt
    
    def antibody_fn(self, t):
        V1 = ((self.PK_params.volume1 * ((self.patient.weight_SUVr/73.7)**self.PK_params.weight_volume1) *
               (self.PK_params.sex_volume1**self.patient.sex)) *
               math.exp(random.normalvariate(0, math.sqrt(self.PK_params.random_volume1))))
        return 0.83 #(((120*1000)/1000000000)/145181)*1000000000#((self.spline(t) * (V1*1000))/1000)


    def __call__(self):
        self.SUVr_fixed(self.patient.APOE, self.patient.age_SUVr)