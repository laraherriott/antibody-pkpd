#
# Class to define the model parameters - use as template to generate other parameter files
# note that these parameters are in time units of hour
#

class Parameters:
    def __init__(self):
        self.clearance = 0.0181
        self.volume1 = 3.22
        self.volume2 = 2.19
        self.inter_compartment = 0.0349
        self.albumin_clearance = -0.243
        self.sex_clearance = 0.792
        self.weight_clearance = 0.403
        self.ADA_clearance = 1.09
        self.weight_volume1 = 0.606
        self.sex_volume1 = 0.893
        self.japanese_volume2 = 0.455
        self.processB = 0.998

        self.random_clearance = 0.2 
        self.random_volume1 = 0.5
        self.random_volume2 = 1
        self.random_F = 1
    
        self.SUVr_baseline = 1.34
        self.SUVr_Kin = 0.232 # note in model /8760 to scale per year (8760 hours per y therefore a timestep is an hour)
        self.SUVr_Emax = 1.54
        self.SUVr_EC50 = 75
        self.APOE_SUVr_baseline = 1.04
        self.age_Emax = 1.58

        self.random_SUVr_baseline = 0.015
        self.random_Emax = 0.05

        self.AB42_40_baseline = 0.0842
        self.AB42_40_Kout = 0.367 # /8760
        self.AB42_40_slope = 0.00155

        self.random_AB42_40_baseline = 0.1
        self.random_AB42_40_slope = 0.1

        self.tau_baseline = 4.06
        self.tau_Kout = 0.468 # /8760
        self.tau_slope = 0.00313
        self.weight_tau_baseline = -0.3 

        self.random_tau_baseline = 0.1
        self.random_tau_slope = 0.1

