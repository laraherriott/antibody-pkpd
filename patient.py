import random

class Patient:

    def __init__(self):
        chance_male = 0.521
        chance_ADA_neg = 0.57
        chance_apoe4 = 0.698
        chance_japanese = 0.069
        chance_processA = 0.946
        mean_age_SUVr = 71.5
        sd_age_SUVr = 8.1
        mean_weight_SUVr = 74.3
        sd_weight_SUVr = 15.1
        mean_albumin = 42.9
        sd_albumin = 2.9
        # ada_titer = {1: 0.187, 2: 0.012, 5: 0.26, 16: 0.016, 25: 0.224, 80: 0.008,
        #              125: 0.191, 400: 0.02, 625: 0.073, 3125: 0.004, 15625: 0.004}
        
        self.age_SUVr = random.normalvariate(mean_age_SUVr, sd_age_SUVr)
        self.weight_SUVr = random.normalvariate(mean_weight_SUVr, sd_weight_SUVr)
        self.albumin = random.normalvariate(mean_albumin, sd_albumin)

        r = random.random()
        if r < chance_male:
            self.sex = 0
        else:
            self.sex = 1
        if r < chance_ADA_neg:
            self.ADA = 0
        else:
            self.ADA = 1
        if r < chance_apoe4:
            self.APOE = 0
        else:
            self.APOE = 1
        if r < chance_japanese:
            self.race = 1
        else:
            self.race = 0
        if r < chance_processA:
            self.process = 0
        else:
            self.process = 1

class FixedPatient:

    def __init__(self):

        self.age_SUVr = 72
        self.weight_SUVr = 74.3
        self.albumin = 42.9
        self.sex = 0
        self.ADA = 0
        self.APOE = 0
        self.race = 0
        self.process = 0