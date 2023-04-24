class AducanumabModel:
    def __init__(self, dose, dose_interval, final_dose):
        self.dose = dose
        self.interval = dose_interval
        self.final_dose = final_dose
        self.params = AducannumabParams()


    def plasma_eq(self, t):
        dose = DoseFn(self.dose, self.interval, self.final_dose, self.patient)
        CAPP = y[0]/self.params.plasma_vol
        CBACE = 
        CAPPBACE = 
        dYdt = []
        dAPP = self.params.k_synth - self.params.k_clear_app*CAPP - self.params.k_onPP*CAPP*CBACE + self.params.k_offBACE*CAPPBACE
        dAPP_BACE
        dBACE
        dBACEs
        dCTFB
        dCTFB_Gamma
        dGamma
        dABeta
        dABeta_mAb
        dOligomer
        dOligomer_mAb

        return dYdt
    
    def csf_eq(self):
        dYdt = []
        return dYdt
    
    def isf_eq(self):
        dYdt = []
        return dYdt
    
    def all_eq(self):
        plasma = self.plasma_eq()
        csf = self.csf_eq()
        isf = self.isf_eq()

        dYdt = [plasma, csf, isf]
        return dYdt