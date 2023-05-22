import math
import numpy as np
import pandas as pd
import random
from  scipy.interpolate import CubicSpline

from parameters.parameters import MITParams
from patient import FixedPatient

class MITModel:
    def __init__(self, model_type):
        self.params = MITParams()

        # self.params.bisf_vol = 0.261

        # # pathway parameters
        # self.params.k_synth = 0.002
        # self.params.k_synth_bace = 0.0001
        # self.params.k_synth_gamma = 0.000052
        
        # self.params.k_clear_app = 0.0000693
        # self.params.k_clear_bace = 0.000012
        # self.params.k_clear_gamma = 0.000008
        # self.params.k_clear_ctf = 0.000055
        
        # self.params.k_onPP = 0.001
        # self.params.k_offBACE = 120
        # self.params.k_off_gamma = 0.4
        # self.params.k_cleave = 0.0000002
        # self.params.k_cat_bace = 0.0072
        # self.params.k_cat_gamma = 0.0012
        
        self.PD = model_type # of 'adu_path' or 'brain'

        self.half_life_dose = [1]

        self.dose_list = []
        i = 0
        while i <= (24*365*360*1.5):
            self.dose_list.append(int(i))
            i += (14*24*360) # lec = 14, adu = 28

    def pathway(self, t, y): # y = [APP, BACE, BACEs]
        App = y[0]/self.params.bisf_vol
        Bace = y[1]/self.params.bisf_vol
        Baces = y[2]/self.params.bisf_vol
        Gamma = y[3]/self.params.bisf_vol
        Ctfb = y[4]/self.params.bisf_vol
        GtfbGamma = y[5]/self.params.bisf_vol
        AppBace = y[6]/self.params.bisf_vol
        #Abeta = y[7]/self.params.bisf_vol


        dAPP = self.params.k_synth - self.params.k_clear_app*App - self.params.k_onPP*App*Bace + self.params.k_offBACE*AppBace
        dBACE = self.params.k_synth_bace - self.params.k_clear_bace*Bace - self.params.k_onPP*App*Bace + self.params.k_offBACE*AppBace - self.params.k_cleave*Bace + self.params.k_cat_bace*AppBace
        dBACEs = self.params.k_cleave*Bace - self.params.k_clear_baces*Baces
        dGamma = self.params.k_synth_gamma - self.params.k_clear_gamma*Gamma + self.params.k_cat_gamma*GtfbGamma
        dCTFB = self.params.k_cat_bace*AppBace - self.params.k_onPP*Ctfb*Gamma + self.params.k_off_gamma*GtfbGamma - self.params.k_clear_ctf*Ctfb
        
        dCTFB_Gamma = self.params.k_onPP*Ctfb*Gamma - self.params.k_off_gamma*GtfbGamma - self.params.k_cat_gamma*GtfbGamma
        dAPP_BACE = self.params.k_onPP*App*Bace - self.params.k_offBACE*AppBace - self.params.k_cat_bace*AppBace

        #dABeta = self.params.k_cat_gamma*СTFBGAMMA - self.params.k_olig_inc*Abeta + self.params.k_olig_sep*Oligomer - self.params.k_onPP*Abeta*mAb + k_off_ma0*ABeta_mAb


        dYdt = [dAPP, dBACE, dBACEs, dGamma, dCTFB, dCTFB_Gamma, dAPP_BACE]
        return dYdt

    def brain(self, t, y):
        Abeta = y[0]#/self.params.bisf_vol
        Oligomer = y[1]#/self.params.bisf_vol
        Plaque = y[2]#/self.params.bisf_vol
        FcR = y[3]#/self.params.bisf_vol
        mAb = y[4]#/self.params.bisf_vol
        ABeta_mAb = y[5]#/self.params.bisf_vol
        Oligomer_mAb = y[6]#/self.params.bisf_vol
        Plaque_mAb = y[7]#/self.params.bisf_vol
        Oligomer_mAb_FcR = y[8]#/self.params.bisf_vol
        Plaque_mAb_FcR = y[9]#/self.params.bisf_vol
        mAb_plasma = y[10]
        #Lec_P = y[10]
        #C1L = (((y[10]/1000/self.params.volume1)/147181.62))*1e9 # nM
        
        dABeta = self.params.k_in - self.params.k_olig_inc*Abeta + self.params.k_olig_sep*Oligomer + self.params.k_off_ma0*ABeta_mAb - self.params.k_clear_Abeta*Abeta - self.params.k_onPP*Abeta*mAb
        dOligomer = self.params.k_olig_inc*Abeta - self.params.k_olig_sep*Oligomer - self.params.k_plaque_inc*Oligomer + self.params.k_plaque_sep*Plaque +self.params.k_off_ma1*Oligomer_mAb - self.params.k_clear_olig*Oligomer - self.params.k_onPP*Oligomer*mAb
        dPlaque =   + self.params.k_off_ma2*Plaque_mAb - self.params.k_clear_P*Plaque + self.params.k_plaque_inc*Oligomer - self.params.k_plaque_sep*Plaque - self.params.k_onPD*Plaque*mAb

        dFcR = self.params.k_synth_FcR - self.params.k_clear_FcR*FcR - self.params.k_onPF*Oligomer_mAb*FcR + self.params.k_offPF*Oligomer_mAb_FcR - self.params.k_onPF*Plaque_mAb*FcR + self.params.k_offPF*Plaque_mAb_FcR + self.params.k_ADCP*Plaque_mAb_FcR + self.params.k_ADCP*Oligomer_mAb_FcR
        dmAb = self.params.mAb_transport*mAb_plasma - self.params.mAb_transport_back*mAb - self.params.k_mAbcomplex_clear*mAb  + self.params.k_off_ma0*ABeta_mAb  + self.params.k_off_ma1*Oligomer_mAb  + self.params.k_off_ma2*Plaque_mAb - self.params.k_onPP*Abeta*mAb - self.params.k_onPP*Oligomer*mAb - self.params.k_onPD*Plaque*mAb

        dABeta_mAb = - self.params.k_off_ma0*ABeta_mAb - self.params.k_mAbcomplex_clear*ABeta_mAb + self.params.k_onPP*Abeta*mAb
        dOligomer_mAb = - self.params.k_off_ma1*Oligomer_mAb - self.params.k_mAbcomplex_clear*Oligomer_mAb - self.params.k_onPF*Oligomer_mAb*FcR + self.params.k_offPF*Oligomer_mAb_FcR + self.params.k_onPP*Oligomer*mAb
        dPlaque_mAb = - self.params.k_off_ma2*Plaque_mAb - self.params.k_onPF*Plaque_mAb*FcR + self.params.k_offPF*Plaque_mAb_FcR + self.params.k_onPD*Plaque*mAb

        dOligomer_mAb_FcR = self.params.k_onPF*Oligomer_mAb*FcR - self.params.k_offPF*Oligomer_mAb_FcR - self.params.k_ADCP*Oligomer_mAb_FcR
        dPlaque_mAb_FcR = self.params.k_onPF*Plaque_mAb*FcR - self.params.k_offPF*Plaque_mAb_FcR - self.params.k_ADCP*Plaque_mAb_FcR

        
        dmAb_plasma = self.dosefn(self.dose_list, t) - self.params.clearance*mAb_plasma + self.params.mAb_transport_back*mAb - self.params.mAb_transport*mAb_plasma
        #dLper_dt = (self.params.inter_compartment/self.params.volume1)*y[4] - (self.params.inter_compartment/self.params.volume2)*y[10]

        dPhagocytosed = self.params.k_ADCP*Plaque_mAb_FcR
        dNewPlaqueAgg = self.params.k_plaque_inc*Oligomer - self.params.k_plaque_sep*Plaque - self.params.k_clear_P*Plaque

        dYdt = [dABeta, dOligomer, dPlaque, dFcR, dmAb, dABeta_mAb, dOligomer_mAb, dPlaque_mAb, dOligomer_mAb_FcR, dPlaque_mAb_FcR, dmAb_plasma, dPhagocytosed, dNewPlaqueAgg]
        return dYdt
    

    def half_life(self, t, y):
        dmAb_plasma = self.dosefn(self.half_life_dose, t) - self.params.clearance*y[0]

        return [dmAb_plasma]
    
    # def csf_eq(self):
    #     dYdt = []
    #     return dYdt
    
    # def isf_eq(self):
    #     dYdt = []
    #     return dYdt
    
    # def all_eq(self):
    #     plasma = self.plasma_eq()
    #     csf = self.csf_eq()
    #     isf = self.isf_eq()

    #     dYdt = [plasma, csf, isf]
    #     return dYdt
    
    def dosefn(self, dose_list, t):
        infusion = (((((10 * 70)/1000/3.22)/147181.62))*1e9)/360 # nM
        f = 0.5
        delta = 0.1

        sol = 0
        for n in dose_list:
            if t >=(n-(0.5*360)) and t<=(n+(1.5*360)):
                sol =  ((infusion/2)/math.atan(1/delta))*(math.atan(math.sin(2*math.pi*(t-n-0.5)*f)/delta)) + infusion/2

        return sol