import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy
from lmfit import minimize, Parameters, report_fit
from scipy.integrate import odeint

k_in = 0.000055#*360
k_clear_Abeta = 0.000055#*360
k_clear_olig = 0.000000022#*360
k_clear_P = 0.000000008#*360
k_synth_FcR = 0.000025#*360
k_clear_FcR = 0.000096#*360

params = Parameters()
params.add('k_in', value=k_in, min=0, max=0.1)
params.add('k_synth_FcR', value=k_synth_FcR, min=0, max=0.1)
params.add('k_clear_FcR', value=k_clear_FcR, min=0, max=0.1)
params.add('k_clear_Abeta', value=k_clear_Abeta, min=0, max=0.1)
params.add('k_clear_olig', value=k_clear_olig, min=0, max=0.1)
params.add('k_clear_P', value=k_clear_P, min=0, max=0.1)

def ode_model(y, t, k_in, k_clear_Abeta, k_clear_olig, k_synth_FcR, k_clear_FcR, k_clear_P):
    k_onPP = 0.001
    k_onPD = 0.001#*360
    k_onPF = 0.001#*360
    k_olig_inc = 0.000015#*360
    k_olig_sep = 0.000000014#*360
    k_plaque_inc = 0.00000007#*360
    k_plaque_sep = 0.00000000007#*360
    k_offPF = k_onPF * 15.6#*360
    k_off_ma0 = 2.29#*360
    k_off_ma1 = 0.0673#*360
    k_off_ma2 = 0.00179#*360
    k_ADCP = 0.000037#*360
    # k_clear_Abeta = 0.000055*360
    # k_clear_olig = 0.000000022*360
    #k_clear_P = 0.000000008#*360
    # k_synth_FcR = 0.000025*360
    # k_clear_FcR = 0.000096*360
    
    Abeta = y[0]#/bisf_vol
    Oligomer = y[1]#/bisf_vol
    Plaque = y[2]#/bisf_vol
    FcR = y[3]#/bisf_vol
    # mAb = y[4]#/bisf_vol
    # ABeta_mAb = y[5]#/bisf_vol
    # Oligomer_mAb = y[6]#/bisf_vol
    # Plaque_mAb = y[7]#/bisf_vol
    # Oligomer_mAb_FcR = y[8]#/bisf_vol
    # Plaque_mAb_FcR = y[9]#/bisf_vol
    
    dABeta = k_in - k_olig_inc*Abeta + k_olig_sep*Oligomer - k_clear_Abeta*Abeta #- k_onPP*Abeta*mAb + k_off_ma0*ABeta_mAb 
    dOligomer = k_olig_inc*Abeta - k_olig_sep*Oligomer - k_plaque_inc*Oligomer + k_plaque_sep*Plaque - k_clear_olig*Oligomer #- k_onPP*Oligomer*mAb +k_off_ma1*Oligomer_mAb 
    dPlaque = k_plaque_inc*Oligomer - k_plaque_sep*Plaque - k_clear_P*Plaque #- k_onPD*Plaque*mAb + k_off_ma2*Plaque_mAb 

    dFcR = k_synth_FcR - k_clear_FcR*FcR #- k_onPF*Oligomer_mAb*FcR + k_offPF*Oligomer_mAb_FcR - k_onPF*Plaque_mAb*FcR + k_offPF*Plaque_mAb_FcR
    # dmAb = k_onPP*Abeta*mAb + k_off_ma0*ABeta_mAb - k_onPP*Oligomer*mAb + k_off_ma1*Oligomer_mAb - k_onPD*Plaque*mAb + k_off_ma2*Plaque_mAb

    # dABeta_mAb = k_onPP*Abeta*mAb - k_off_ma0*ABeta_mAb
    # dOligomer_mAb = k_onPP*Oligomer*mAb - k_off_ma1*Oligomer_mAb
    # dPlaque_mAb = k_onPD*Plaque*mAb - k_off_ma2*Plaque_mAb

    # dOligomer_mAb_FcR = k_onPF*Oligomer_mAb*FcR - k_offPF*Oligomer_mAb_FcR - k_ADCP*Oligomer_mAb_FcR
    # dPlaque_mAb_FcR = k_onPF*Plaque_mAb*FcR - k_offPF*Plaque_mAb_FcR - k_ADCP*Plaque_mAb_FcR

    dYdt = [dABeta, dOligomer, dPlaque, dFcR]#, dmAb, dABeta_mAb, dOligomer_mAb, dPlaque_mAb, dOligomer_mAb_FcR, dPlaque_mAb_FcR]
    return dYdt

def ode_solver(t, initial_conditions, params):
    initAbeta, initOlig, initPlaque, initFcR = initial_conditions
    k_in, k_clear_Abeta, k_clear_olig, k_synth_FcR, k_clear_FcR, k_clear_P = params['k_in'].value, params['k_clear_Abeta'].value, params['k_clear_olig'].value, params['k_synth_FcR'].value, params['k_clear_FcR'].value, params['k_clear_P'].value
    res = odeint(ode_model, [initAbeta, initOlig, initPlaque, initFcR], t, args=(k_in, k_clear_Abeta, k_clear_olig, k_synth_FcR, k_clear_FcR, k_clear_P))
    return res

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

initial_conditions = [initAbeta, initOlig, initPlaque, initFcR]#, initmAb, initAbetamAb, initOligmAb, initPlaquemAb, initOligmAbFcR, initPlaquemAbFcR]

result = minimize(error, params, args=(initial_conditions, tspan, data), method='leastsq')


print(result.params)
print(report_fit(result))
