import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy
import statistics
from lmfit import minimize, Parameters, report_fit
from scipy.integrate import odeint
from scipy import stats
from SALib.sample import fast_sampler
from SALib.analyze import fast
from datetime import datetime

startTime = datetime.now()

k_in = 0.000055*360 # take initial value to balance baseline clearance (expect to get higher to account for binding and aggregation)
k_clear_Abeta = 0.000055*360 # Cirrito et al 2003
k_clear_olig = 0.000000022*360  #BioMath fitted value (ref 1,3 SILK from adu paper)
k_clear_P = 0.00000000441*360 # ref 28 from Adu paper
k_synth_FcR = (0.0000503/0.261)*360 # nmol/s from Adu paper (so FcR not limiting) and convert to nM by dividing by vol ISF (Shah & Betts 2012)
k_clear_FcR = 0.000193*360 # ref 29 typical receptor turnover
k_onPP = 0.001*360 # this and below typical protein-protein association
k_onPD = 0.001*360
k_onPF = 0.001*360
k_olig_inc = 0.000015*360 # Adu ref 21,22 (in vitro oligomerization take up to 20 hrs)
k_olig_sep = 0.000000014*360 # fitted to SILK (ref 1,3 from adu)
k_plaque_inc = 0.00000007*360 # fitted to SILK (ref 1,3)
k_plaque_sep = 0.00000000007*360 # assumed 1000x slower than formation
k_ADCP = 0.0036*360 # fitted to SUVr
k_mAbcomplex_clear = 0.00000015*360
k_mAb_transport = 0.0000016*360 # fitted to aducanumab CSF PK
k_mAb_transport_back = 0.0032*360 # fitted to aducanumab CSF PK
clearance = 0.00003*360 #0.00000146 # fitted to aducanumab PK


# lecanemab
k_off_ma0 = 2290 * k_onPP * 360
k_off_ma1 = 67.3 * k_onPP *360
k_off_ma2 = 1.79 * k_onPD *360
k_offPF = 0.12 * k_onPF*360

def ode_model(t, y, p, dose_list):
    k_in = p[0] # take initial value to balance baseline clearance (expect to get higher to account for binding and aggregation)
    #k_clear_Abeta = p[1] # Cirrito et al 2003
    k_clear_olig = p[2]#2.2e-8*360  #BioMath fitted value (ref 1,3 SILK from adu paper)
    k_clear_P = p[1]#4.41e-9*360 # ref 28 from Adu paper
    #k_synth_FcR = p[5] # nmol/s from Adu paper (so FcR not limiting) and convert to nM by dividing by vol ISF (Shah & Betts 2012)
    #k_clear_FcR = p[2] # ref 29 typical receptor turnover
    #k_onPP = p[6] # this and below typical protein-protein association
    #k_onPD = p[7]
    #k_onPF = p[8]
    #k_olig_inc = p[9] # Adu ref 21,22 (in vitro oligomerization take up to 20 hrs)
    #k_olig_sep = p[10] # fitted to SILK (ref 1,3 from adu)
    k_plaque_inc = p[3]#7e-8*360 # fitted to SILK (ref 1,3)
    #k_plaque_sep = p[12] # assumed 1000x slower than formation
    k_ADCP = p[4]#0.0036*360 # fitted to SUVr
    #k_mAbcomplex_clear = p[14]
    #k_mAb_transport = p[15] # fitted to aducanumab CSF PK
    #k_mAb_transport_back = p[16] # fitted to aducanumab CSF PK
    clearance = p[5]#1.1e-5*360 #0.00000146 # fitted to aducanumab PK

    # lecanemab
    k_off_ma0 = p[6]#2290*0.001*360
    k_off_ma1 = p[7]#67.3*0.001*360
    k_off_ma2 = p[8]#1.79*0.001*360
    k_onPF = 0.001*360
    k_offPF = 0.12 * k_onPF * 360
    k_mAbcomplex_clear = 0.00000015*360
    k_mAb_transport = 0.0000016*360 # fitted to aducanumab CSF PK
    k_mAb_transport_back = 0.0032*360
    k_clear_Abeta = 0.000055*360 # Cirrito et al 2003
    k_synth_FcR = (0.0000503/0.261)*360 # nmol/s from Adu paper (so FcR not limiting) and convert to nM by dividing by vol ISF (Shah & Betts 2012)
    k_clear_FcR = 0.000193*360 # ref 29 typical receptor turnover
    k_onPP = 0.001*360 # this and below typical protein-protein association
    k_onPD = 0.001*360
    
    k_olig_inc = 0.000015*360 # Adu ref 21,22 (in vitro oligomerization take up to 20 hrs)
    k_olig_sep = 0.000000014*360 # fitted to SILK (ref 1,3 from adu)
    k_plaque_sep = 0.00000000007*360

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
    t_span = np.arange(0, (24*56), (24)) # one data point per day

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
    while i <= (24*56):
        dose_list.append(int(i))
        i += (14*24)

    res = scipy.integrate.solve_ivp(fun=lambda t, y: ode_model(t, y, p, dose_list),
                                             t_span=[t_span[0], t_span[-1]],
                                             y0=initial_conditions,
                                             t_eval=t_span,
                                             max_step = 0.1, 
                                             method = 'LSODA')
    return res

def dosefn(dose_list, t):
    infusion = (((((10 * 70)/1000/3.22)/147181.62))*1e9) # nM
    f = 0.5
    delta = 0.1

    sol = 0
    for n in dose_list:
        if t >=(n-(0.5)) and t<=(n+(1.5)):
            sol =  ((infusion/2)/math.atan(1/delta))*(math.atan(math.sin(2*math.pi*(t-n-0.5)*f)/delta)) + infusion/2

    return sol

initial_guess = [k_in, k_clear_P, 
                 k_clear_olig, k_plaque_inc,  
                 k_ADCP, clearance, 
                 k_off_ma0, k_off_ma1, k_off_ma2]

bounds = []
for i in range(len(initial_guess)):
    lower = initial_guess[i]*1e-2
    upper = initial_guess[i]*1e2
    sample_within = [lower, upper]
    bounds.append(sample_within)

problem = {
  'num_vars': 9,
  'names': ['k_in', 'k_clear_P', 'k_clear_olig',  'k_plaque_inc', 'k_ADCP', 'clearance', 'k_off_ma0', 'k_off_ma1', 'k_off_ma2'],
  'bounds': bounds
}

bounds[6] = [0, initial_guess[6]*1e2]
bounds[7] = [0, initial_guess[6]*1e2]
bounds[8] = [0, initial_guess[6]*1e2]

# Generate samples
param_values = fast_sampler.sample(problem, 100)
N = len(param_values) # number of parameter samples
Y = np.zeros([N,1])
#print(N)
# Run model for each parameter set, save the output in array Y
for i in range(N):
  #print(i)

  sol = ode_solver(param_values[i])

  #plaque = [(sol.y[2][i] + sol.y[7][i] + sol.y[9][i]) for i in range(len(sol.y[2]))]
  antibody = statistics.mean(sol.y[10])

  Y[i,0] = antibody

# Perform analysis
Si = fast.analyze(problem, Y, print_to_console=True)#, calc_second_order=False)

df = Si.to_df()

df.to_csv('sensitivity_out_other_plasma.csv')

Si.plot()
fig = plt.gcf()
plt.savefig('sensitivity_out_other_plasma.png')