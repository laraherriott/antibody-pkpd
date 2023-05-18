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

class MITParams:
    def __init__(self):
        self.bisf_vol = 0.261

        # pathway parameters
        self.k_synth = 0.002
        self.k_synth_bace = 0.0001
        self.k_synth_gamma = 0.000052
        
        self.k_clear_app = 0.0000693
        self.k_clear_bace = 0.000012
        self.k_clear_baces = 0.000012 # not given, assuming same as above
        self.k_clear_gamma = 0.000008
        self.k_clear_ctf = 0.000055
        
        self.k_onPP = 0.001#*360
        self.k_offBACE = 120
        self.k_off_gamma = 0.4
        self.k_cleave = 0.0000002
        self.k_cat_bace = 0.0072
        self.k_cat_gamma = 0.0012

        # brain parameters

        
        self.k_synth_FcR = 0.0000503/0.261#*360 # = 0.009000000000000001 vs 0.01086183411851932 fitted
        self.k_clear_FcR = 0.000193#*360 # = 0.03456 vs 0.010861834118586933 fitted
        self.k_onPD = 0.001#*360
        self.k_onPF = 0.001#*360
        self.k_olig_inc = 0.000015#*360
        self.k_olig_sep = 0.000000014#*360
        self.k_plaque_inc = 0.00000007#*360
        self.k_plaque_sep = 0.00000000007#*360
        
        self.k_Abetaplasma = 0.000000148#*360
        self.k_AbetaBrain = 0.0000000148#*360
        self.k_oligplasma = 0.000000148#*360
        self.k_oligBrain = 0.0000000148#*360
        self.k_AbetaCSF = 0.000075#*360
        self.k_oligCSF = 0.00000225#*360
        self.k_ADCP = 0.0036#*360
        self.k_offPF = self.k_onPF * 0.12

        # antibody-specific parameters
        # lecanemab
        self.k_off_ma0 = 2.29#*360
        self.k_off_ma1 = 0.0673#*360
        self.k_off_ma2 = 0.00179#*360

        # aducanumab
        # self.k_off_ma0 = 8#*360
        # self.k_off_ma1 = 0.004#*360
        # self.k_off_ma2 = 0.004#*360

        # #"donanemab"/plaque only
        # self.k_off_ma0 = 0
        # self.k_off_ma1 = 0
        # self.k_off_ma2 = 0.0015

        # # oligomer only
        # self.k_off_ma0 = 0
        # self.k_off_ma1 = 0.0015
        # self.k_off_ma2 = 0

        # monomer only
        # self.k_off_ma0 = 0.0015
        # self.k_off_ma1 = 0
        # self.k_off_ma2 = 0

        # un-fitted
        self.k_in = 0.000055#*360 # = 0.0198 vs 0.004523823527794346 fitted
        self.k_clear_Abeta = 0.000055#*360 # 0.0198 = vs 0.00066133270361457 fitted
        self.k_clear_olig = 0.000000022#*360 # 0.00000791 = vs 0.000890437909069397 fitted
        self.k_clear_P = 0.00000000441#*360 # 0.00000288= vs 0.000002817876923133511 fitted

        # new - fit with restrictions to keep order of magnitude
        # self.k_in = 0.0636#0.000055919 # = 0.0198 vs 0.004523823527794346 fitted
        # self.k_clear_Abeta = 0.0847#0.000059641 # 0.0198 = vs 0.00066133270361457 fitted
        # self.k_clear_olig = 0.0000024734 # 0.00000791 = vs 0.000890437909069397 fitted
        # self.k_clear_P = 0.0000000078274 # 0.00000288= vs 0.000002817876923133511 fitted
        # self.k_synth_FcR = 0.0000983#0.00006131 # = 0.009000000000000001 vs 0.01086183411851932 fitted
        # self.k_clear_FcR = 0.0000983#0.000006131 # = 0.03456 vs 0.010861834118586933 fitted

        # fitted but not for FcR for 3 yr:
        # self.k_in = 0.4469159655135057
        # self.k_clear_Abeta = 0.5905175220791339
        # self.k_clear_olig = 0.0008904379091214665
        # self.k_clear_P = 0.000002817876923078
        # moved away from this since FcR was not staying constant in this regimen when it should certainly be...

        # fitted inc for FcR
        # self.k_in = 0.004523823527794346
        # self.k_clear_Abeta = 0.00066133270361457
        # self.k_clear_olig = 0.000890437909069397
        # self.k_clear_P = 0.000002817876923133511
        # self.k_synth_FcR = 0.01086183411851932
        # self.k_clear_FcR = 0.010861834118586933
        # self.k_mAbcomplex_clear = 0.00000015*360
        # self.mAb_transport = 0.00000000172 * 360
        # self.mAb_transport_back = 0.00000810 *360


        # fitted for higher FcR (10nm, so not limiting)
        # self.k_in = 0.196
        # self.k_clear_Abeta = 0.256
        # self.k_clear_olig = 0.000890437909069397
        # self.k_clear_P = 0.000002817876923133511
        # self.k_synth_FcR = 1
        # self.k_clear_FcR = 0.1
        # self.k_mAbcomplex_clear = 0.00000015*360
        # self.mAb_transport = 0.00000000172 * 360
        # self.mAb_transport_back = 0.00000810 *360

        # fitted with equations not including mAb etc
        # self.k_in = 0.00999775
        # self.k_clear_Abeta = 0.00796
        # self.k_clear_olig = 0.0008881
        # self.k_clear_P = 0.000000008*360
        # self.k_synth_FcR = 0.09999
        # self.k_clear_FcR = 0.09999
        self.k_mAbcomplex_clear = 0.00000015#*360
        self.mAb_transport = 0.0000016 #* 360 NB this value now slightly increased.
        self.mAb_transport_back = 0.0032 #*360

        # fitted to Abeta monomer only?

        # fitted to PK value, otherwise clearance =  0.0181 OR 0.00000015*360?
        self.clearance = 0.00000146
        self.volume1 = 3.22
        self.volume2 = 2.19
        self.inter_compartment = 0.0349

