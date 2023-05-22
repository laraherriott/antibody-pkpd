#
# Class to define the model parameters
#

class NoIIVParameters:
    def __init__(self):
        self.clearance = 0.0181 / 360
        self.volume1 = 3.22
        self.volume2 = 2.19
        self.inter_compartment = 0.0349 /360
        self.albumin_clearance = -0.243
        self.sex_clearance = 0.792
        self.weight_clearance = 0.403
        self.ADA_clearance = 1.09
        self.weight_volume1 = 0.606
        self.sex_volume1 = 0.893
        self.japanese_volume2 = 0.455
        self.processB = 0.998

        # variance values for adding inter-individual variability on the parameters
        self.random_clearance = 0
        self.random_volume1 = 0
        self.random_volume2 = 0
        self.random_F = 0
    
        self.SUVr_baseline = 1.34
        self.SUVr_Kin = 0.232 # note in model /8760 to scale per year (8760 hours per y therefore a timestep is an hour)
        self.SUVr_Emax = 1.54
        self.SUVr_EC50 = 75
        self.APOE_SUVr_baseline = 1.04
        self.age_Emax = 1.58

        # variance values for adding inter-individual variability on the parameters
        self.random_SUVr_baseline = 0
        self.random_Emax = 0

        self.AB42_40_baseline = 0.0842
        self.AB42_40_Kout = 0.367 # /8760
        self.AB42_40_slope = 0.00155

        # variance values for adding inter-individual variability on the parameters
        self.random_AB42_40_baseline = 0
        self.random_AB42_40_slope = 0

        self.tau_baseline = 4.06
        self.tau_Kout = 0.468 # /8760
        self.tau_slope = 0.00313
        self.weight_tau_baseline = -0.3 

        # variance values for adding inter-individual variability on the parameters
        self.random_tau_baseline = 0
        self.random_tau_slope = 0

        # variance values for proportional residual variability

