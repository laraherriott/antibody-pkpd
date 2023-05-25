import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy
from lmfit import minimize, Parameters, report_fit
from scipy.integrate import odeint
from scipy import stats
from SALib.sample import latin
from SALib.analyze import sobol
from datetime import datetime

startTime = datetime.now()

k_in = 0.000055 # take initial value to balance baseline clearance (expect to get higher to account for binding and aggregation)
k_clear_Abeta = 0.000055 # Cirrito et al 2003
k_clear_olig = 0.000000022  #BioMath fitted value (ref 1,3 SILK from adu paper)
k_clear_P = 0.00000000441 # ref 28 from Adu paper
k_synth_FcR = 0.0000503/0.261 # nmol/s from Adu paper (so FcR not limiting) and convert to nM by dividing by vol ISF (Shah & Betts 2012)
k_clear_FcR = 0.000193 # ref 29 typical receptor turnover
k_onPP = 0.001 # this and below typical protein-protein association
k_onPD = 0.001
k_onPF = 0.001
k_olig_inc = 0.000015 # Adu ref 21,22 (in vitro oligomerization take up to 20 hrs)
k_olig_sep = 0.000000014 # fitted to SILK (ref 1,3 from adu)
k_plaque_inc = 0.00000007 # fitted to SILK (ref 1,3)
k_plaque_sep = 0.00000000007 # assumed 1000x slower than formation
k_ADCP = 0.0036 # fitted to SUVr
k_mAbcomplex_clear = 0.00000015
k_mAb_transport = 0.0000016 # fitted to aducanumab CSF PK
k_mAb_transport_back = 0.0032 # fitted to aducanumab CSF PK
clearance = 0.00003 #0.00000146 # fitted to aducanumab PK


# lecanemab
k_off_ma0 = 2290 * k_onPP
k_off_ma1 = 67.3 * k_onPP
k_off_ma2 = 1.79 * k_onPD
k_offPF = 0.12 * k_onPF

def ode_model(t, y, p, dose_list):
    k_in = p[0] # take initial value to balance baseline clearance (expect to get higher to account for binding and aggregation)
    k_clear_Abeta = p[1] # Cirrito et al 2003
    k_clear_olig = p[4]  #BioMath fitted value (ref 1,3 SILK from adu paper)
    k_clear_P = p[3] # ref 28 from Adu paper
    k_synth_FcR = p[5] # nmol/s from Adu paper (so FcR not limiting) and convert to nM by dividing by vol ISF (Shah & Betts 2012)
    k_clear_FcR = p[2] # ref 29 typical receptor turnover
    k_onPP = p[6] # this and below typical protein-protein association
    k_onPD = p[7]
    k_onPF = p[8]
    k_olig_inc = p[9] # Adu ref 21,22 (in vitro oligomerization take up to 20 hrs)
    k_olig_sep = p[10] # fitted to SILK (ref 1,3 from adu)
    k_plaque_inc = p[11] # fitted to SILK (ref 1,3)
    k_plaque_sep = p[12] # assumed 1000x slower than formation
    k_ADCP = p[13] # fitted to SUVr
    k_mAbcomplex_clear = p[14]
    k_mAb_transport = p[15] # fitted to aducanumab CSF PK
    k_mAb_transport_back = p[16] # fitted to aducanumab CSF PK
    clearance = p[17] #0.00000146 # fitted to aducanumab PK

    # lecanemab
    k_off_ma0 = p[18]
    k_off_ma1 = p[19]
    k_off_ma2 = p[20]
    k_offPF = p[21]
    
    Abeta = y[0]#/bisf_vol
    Oligomer = y[1]#/bisf_vol
    Plaque = y[2]#/bisf_vol
    FcR = y[3]#/bisf_vol
    mAb = y[4]#/bisf_vol
    ABeta_mAb = y[5]#/bisf_vol
    Oligomer_mAb = y[6]#/bisf_vol
    Plaque_mAb = y[7]#/bisf_vol
    Oligomer_mAb_FcR = y[8]#/bisf_vol
    Plaque_mAb_FcR = y[9]#/bisf_vol
    mAb_plasma = y[10]
    
    dABeta = k_in - k_olig_inc*Abeta + k_olig_sep*Oligomer - k_clear_Abeta*Abeta - k_onPP*Abeta*mAb + k_off_ma0*ABeta_mAb 
    dOligomer = k_olig_inc*Abeta - k_olig_sep*Oligomer - k_plaque_inc*Oligomer + k_plaque_sep*Plaque - k_clear_olig*Oligomer - k_onPP*Oligomer*mAb +k_off_ma1*Oligomer_mAb 
    dPlaque = k_plaque_inc*Oligomer - k_plaque_sep*Plaque - k_clear_P*Plaque - k_onPD*Plaque*mAb + k_off_ma2*Plaque_mAb 

    dFcR = k_synth_FcR - k_clear_FcR*FcR - k_onPF*Oligomer_mAb*FcR + k_offPF*Oligomer_mAb_FcR - k_onPF*Plaque_mAb*FcR + k_offPF*Plaque_mAb_FcR + k_ADCP*Oligomer_mAb_FcR + k_ADCP*Plaque_mAb_FcR
    
    dmAb_plasma = dosefn(dose_list, t) - clearance*mAb_plasma + k_mAb_transport_back*mAb - k_mAb_transport*mAb_plasma
    dmAb = - k_mAb_transport_back*mAb + k_mAb_transport*mAb_plasma - k_onPP*Abeta*mAb + k_off_ma0*ABeta_mAb - k_onPP*Oligomer*mAb + k_off_ma1*Oligomer_mAb - k_onPD*Plaque*mAb + k_off_ma2*Plaque_mAb - k_mAbcomplex_clear*mAb

    dABeta_mAb = k_onPP*Abeta*mAb - k_off_ma0*ABeta_mAb - k_mAbcomplex_clear*ABeta_mAb
    dOligomer_mAb = k_onPP*Oligomer*mAb - k_off_ma1*Oligomer_mAb - k_mAbcomplex_clear*Oligomer_mAb + k_offPF*Oligomer_mAb_FcR - k_onPF*Oligomer_mAb*FcR
    dPlaque_mAb = k_onPD*Plaque*mAb - k_off_ma2*Plaque_mAb - k_onPF*Plaque_mAb*FcR + k_offPF*Plaque_mAb_FcR

    dOligomer_mAb_FcR = k_onPF*Oligomer_mAb*FcR - k_offPF*Oligomer_mAb_FcR - k_ADCP*Oligomer_mAb_FcR
    dPlaque_mAb_FcR = k_onPF*Plaque_mAb*FcR - k_offPF*Plaque_mAb_FcR - k_ADCP*Plaque_mAb_FcR

    dYdt = [dABeta, dOligomer, dPlaque, dFcR, dmAb, dABeta_mAb, dOligomer_mAb, dPlaque_mAb, dOligomer_mAb_FcR, dPlaque_mAb_FcR, dmAb_plasma]
    return dYdt

def ode_solver(p):
    t_span = np.arange(0, (24*56*360), (360*24)) # one data point per day

    initAbeta = 0.2
    initOlig = 370
    initPlaque = 5500
    initFcR = 1
    initmAb = 0
    initAbetamAb = 0
    initOligmAb = 0
    initPlaquemAb = 0
    initOligmAbFcR = 0
    initPlaquemAbFcR = 0
    initPlasmamAb = 0

    initial_conditions = [initAbeta, initOlig, initPlaque, initFcR, initmAb, initAbetamAb, initOligmAb, initPlaquemAb, initOligmAbFcR, initPlaquemAbFcR, initPlasmamAb]

    dose_list = []
    i = 0
    while i <= (24*56*360):
        dose_list.append(int(i))
        i += (14*24*360)

    res = scipy.integrate.solve_ivp(fun=lambda t, y: ode_model(t, y, p, dose_list),
                                             t_span=[t_span[0], t_span[-1]],
                                             y0=initial_conditions,
                                             t_eval=t_span)
    return res

def dosefn(dose_list, t):
    infusion = (((((10 * 70)/1000/3.22)/147181.62))*1e9)/360 # nM
    f = 0.5
    delta = 0.1

    sol = 0
    for n in dose_list:
        if t >=(n-(0.5*360)) and t<=(n+(1.5*360)):
            sol =  ((infusion/2)/math.atan(1/delta))*(math.atan(math.sin(2*math.pi*(t-n-0.5)*f)/delta)) + infusion/2

    return sol

initial_guess = [k_in, k_clear_Abeta, k_clear_FcR, k_clear_P, 
                 k_clear_olig, k_synth_FcR, k_onPP, k_onPD, k_onPF, 
                 k_olig_inc, k_olig_sep, k_plaque_inc, k_plaque_sep, 
                 k_ADCP, k_mAbcomplex_clear, k_mAb_transport, 
                 k_mAb_transport_back, clearance, 
                 k_off_ma0, k_off_ma1, k_off_ma2, k_offPF]

bounds = []
for i in initial_guess:
    lower = initial_guess[i]*0.9
    upper = initial_guess[i]*1.1
    sample_within = [lower, upper]
    bounds.append(sample_within)

problem = {
  'num_vars': 22,
  'names': ['k_in', 'k_clear_Abeta', 'k_clear_FcR', 'k_clear_P', 
                 'k_clear_olig', 'k_synth_FcR', 'k_onPP', 'k_onPD', 'k_onPF', 
                 'k_olig_inc', 'k_olig_sep', 'k_plaque_inc', 'k_plaque_sep', 
                 'k_ADCP', 'k_mAbcomplex_clear', 'k_mAb_transport', 
                 'k_mAb_transport_back', 'clearance', 
                 'k_off_ma0', 'k_off_ma1', 'k_off_ma2', 'k_offPF'],
  'bounds': bounds
}

# Generate samples
param_values = latin.sample(problem, 1000)
N = len(param_values) # number of parameter samples
Y = np.zeros(N)

# Run model for each parameter set, save the output in array Y
for i in range(N):
  if i % 1000 == 0:
    print(i)

  Y[i] = ode_solver(param_values[i])

# Perform analysis
Si = sobol.analyze(problem, Y, print_to_console=True)