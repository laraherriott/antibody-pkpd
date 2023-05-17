import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy
from lmfit import minimize, Parameters, report_fit
from scipy.integrate import odeint

# all units in nM/s (I think - but what does it mean when the units are quoted as /s? Am I not supposed to multiply by conc?)

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
clearance = 0.00000146 # fitted to aducanumab PK

# lecanemab
k_off_ma0 = 2290 * k_onPP
k_off_ma1 = 67.3 * k_onPP
k_off_ma2 = 1.79 * k_onPD
k_offPF = 0.12 * k_onPF

# when fitting, many of these values can now roughly be trusted, so should be reasonably tightly constrained
# the two parameters which really do need to be fitted and are the real 'unknowns' are k_in and k_ADCP
# note though that in a more complicated model which takes the full synthesis pathway into account k_in could be fitted for separately, as were the aggregation params

params = Parameters()
params.add('k_in', value=k_in, min=0, max=0.1)
params.add('k_synth_FcR', value=k_synth_FcR, min=0, max=0.1)
params.add('k_clear_FcR', value=k_clear_FcR, min=0, max=0.1)
params.add('k_clear_Abeta', value=k_clear_Abeta, min=0, max=0.1)
params.add('k_clear_olig', value=k_clear_olig, min=0, max=0.1)
params.add('k_clear_P', value=k_clear_P, min=0, max=0.1)
# need to add all params to here and as inputs to ode_model function

def ode_model(y, t, k_in, k_clear_Abeta, k_clear_olig, k_synth_FcR, k_clear_FcR, k_clear_P):
    
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
    
    dABeta = k_in - k_olig_inc*Abeta + k_olig_sep*Oligomer - k_clear_Abeta*Abeta - k_onPP*Abeta*mAb + k_off_ma0*ABeta_mAb 
    dOligomer = k_olig_inc*Abeta - k_olig_sep*Oligomer - k_plaque_inc*Oligomer + k_plaque_sep*Plaque - k_clear_olig*Oligomer - k_onPP*Oligomer*mAb +k_off_ma1*Oligomer_mAb 
    dPlaque = k_plaque_inc*Oligomer - k_plaque_sep*Plaque - k_clear_P*Plaque - k_onPD*Plaque*mAb + k_off_ma2*Plaque_mAb 

    dFcR = k_synth_FcR - k_clear_FcR*FcR - k_onPF*Oligomer_mAb*FcR + k_offPF*Oligomer_mAb_FcR - k_onPF*Plaque_mAb*FcR + k_offPF*Plaque_mAb_FcR + k_ADCP*Oligomer_mAb_FcR + k_ADCP*Plaque_mAb_FcR
    
    dmAb_plasma = dosefn(dose_list, t) - clearance*mAb_plasma + mAb_transport_back*mAb - mAb_transport*mAb_plasma
    dmAb = -mAb_transport_back*mAb + mAb_transport*mAb_plasma - k_onPP*Abeta*mAb + k_off_ma0*ABeta_mAb - k_onPP*Oligomer*mAb + k_off_ma1*Oligomer_mAb - k_onPD*Plaque*mAb + k_off_ma2*Plaque_mAb - k_mAbcomplex_clear*mAb

    dABeta_mAb = k_onPP*Abeta*mAb - k_off_ma0*ABeta_mAb - k_mAbcomplex_clear*ABeta_mAb
    dOligomer_mAb = k_onPP*Oligomer*mAb - k_off_ma1*Oligomer_mAb - k_mAbcomplex_clear*Oligomer_mAb + k_offPF*Oligomer_mAb_FcR - k_onPF*Oligomer_mAb*FcR
    dPlaque_mAb = k_onPD*Plaque*mAb - k_off_ma2*Plaque_mAb - k_onPF*Plaque_mAb*FcR + k_offPF*Plaque_mAb_FcR

    dOligomer_mAb_FcR = k_onPF*Oligomer_mAb*FcR - k_offPF*Oligomer_mAb_FcR - k_ADCP*Oligomer_mAb_FcR
    dPlaque_mAb_FcR = k_onPF*Plaque_mAb*FcR - k_offPF*Plaque_mAb_FcR - k_ADCP*Plaque_mAb_FcR

    dYdt = [dABeta, dOligomer, dPlaque, dFcR, dmAb, dABeta_mAb, dOligomer_mAb, dPlaque_mAb, dOligomer_mAb_FcR, dPlaque_mAb_FcR]
    return dYdt

def ode_solver(t, initial_conditions, params):
    initAbeta, initOlig, initPlaque, initFcR = initial_conditions # need to add to initial conditions and to parameters being read in
    k_in, k_clear_Abeta, k_clear_olig, k_synth_FcR, k_clear_FcR, k_clear_P = params['k_in'].value, params['k_clear_Abeta'].value, params['k_clear_olig'].value, params['k_synth_FcR'].value, params['k_clear_FcR'].value, params['k_clear_P'].value
    res = odeint(ode_model, [initAbeta, initOlig, initPlaque, initFcR], t, args=(k_in, k_clear_Abeta, k_clear_olig, k_synth_FcR, k_clear_FcR, k_clear_P))
    return res
#need to develop this function to directly compare plasma mAb to PK profile and to fit for the % reduction in plaque as can be measured from SUVr
def error(params, initial_conditions, tspan, data):
    sol = ode_solver(tspan, initial_conditions, params)
    return (sol[:, 0:4] - data).ravel()

months = 12
tspan = np.arange(0, (12*28*24*360), 1)
Abeta = [0.75]*(12*28*24*360)
Olig = [4.4]*(12*28*24*360)
Plaque = [39]*(12*28*24*360)
Receptor = [1]*(12*28*24*360)
df = pd.DataFrame({'monomer': Abeta,
                     'oligomer': Olig,
                     'plaque': Plaque,
                     'receptor': Receptor
                    })
data = df.values

df2 = pd.read_csv('PK_fit_data.csv')
data2 = df2['PK'].values
times2 = df2['time'].values

# need to decide what input concentrations to use for each species - could try with these and then also with the values used in the other model (much higher)
initAbeta = 0.75
initOlig = 4.4
initPlaque = 39
initFcR = 1
initmAb = 0
initAbetamAb = 0
initOligmAb = 0
initPlaquemAb = 0
initOligmAbFcR = 0
initPlaquemAbFcR = 0

initial_conditions = [initAbeta, initOlig, initPlaque, initFcR, initmAb, initAbetamAb, initOligmAb, initPlaquemAb, initOligmAbFcR, initPlaquemAbFcR]

result = minimize(error, params, args=(initial_conditions, tspan, data), method='leastsq')


print(result.params)
print(report_fit(result))
