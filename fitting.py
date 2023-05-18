import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy
from lmfit import minimize, Parameters, report_fit
from scipy.integrate import odeint
from scipy import stats

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

dose_list = []
i = 0
while i <= (24*100*360):
    dose_list.append(int(i))
    i += (14*24*360)

# when fitting, many of these values can now roughly be trusted, so should be reasonably tightly constrained
# the two parameters which really do need to be fitted and are the real 'unknowns' are k_in and k_ADCP
# note though that in a more complicated model which takes the full synthesis pathway into account k_in could be fitted for separately, as were the aggregation params

params = Parameters()
params.add('k_in', value=k_in, min=0, max=0.1)
params.add('k_clear_Abeta', value=k_clear_Abeta, min=0, max=0.1)
params.add('k_clear_olig', value=k_clear_olig, min=0, max=0.1)
params.add('k_clear_P', value=k_clear_P, min=0, max=0.1)
params.add('k_synth_FcR', value=k_synth_FcR, min=0, max=0.1)
params.add('k_clear_FcR', value=k_clear_FcR, min=0, max=0.1)
params.add('k_onPP', value=k_onPP, vary=False)#, min=0, max=0.1)
params.add('k_onPD', value=k_onPD, vary=False)#, min=0, max=0.1)
params.add('k_onPF', value=k_onPF, vary=False)#, min=0, max=0.1)
params.add('k_olig_inc', value=k_olig_inc, min=0, max=0.1)
params.add('k_olig_sep', value=k_olig_sep, min=0, max=0.1)
params.add('k_plaque_inc', value=k_plaque_inc, min=0, max=0.1)
params.add('k_plaque_sep', value=k_plaque_sep, min=0, max=0.1)
params.add('k_ADCP', value=k_ADCP, min=0, max=1)
params.add('k_mAbcomplex_clear', value=k_mAbcomplex_clear, min=0, max=0.1)
params.add('k_mAb_transport', value=k_mAb_transport, min=0, max=0.1)
params.add('k_mAb_transport_back', value=k_mAb_transport_back, min=0, max=0.1)
params.add('clearance', value=clearance, min=0, max=1)
params.add('k_off_ma0', value=k_off_ma0, vary=False)
params.add('k_off_ma1', value=k_off_ma1, vary=False)
params.add('k_off_ma2', value=k_off_ma2, vary=False)
params.add('k_offPF', value=k_offPF, vary=False)

def ode_model(t, y, k_in, k_clear_Abeta, k_clear_FcR, k_clear_P, k_clear_olig, k_synth_FcR, 
              k_onPP, k_onPD, k_onPF, k_olig_inc, k_olig_sep, k_plaque_inc, k_plaque_sep, k_ADCP,
              k_mAbcomplex_clear, k_mAb_transport, k_mAb_transport_back, clearance, k_off_ma0,
              k_off_ma1, k_off_ma2, k_offPF, dose_list):
    
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

def ode_solver(t, initial_conditions, params, dose_list):
    #initAbeta, initOlig, initPlaque, initFcR, initmAb, initAbetamAb, initOligmAb, initPlaquemAb, initOligmAbFcR, initPlaquemAbFcR, initPlasmamAb = initial_conditions # need to add to initial conditions and to parameters being read in
    k_in, k_clear_Abeta, k_clear_olig, k_synth_FcR, k_clear_FcR, k_clear_P = params['k_in'].value, params['k_clear_Abeta'].value, params['k_clear_olig'].value, params['k_synth_FcR'].value, params['k_clear_FcR'].value, params['k_clear_P'].value
    k_onPP, k_onPD, k_onPF, k_olig_inc, k_olig_sep  = params['k_onPP'].value, params['k_onPD'].value, params['k_onPF'].value, params['k_olig_inc'].value, params['k_olig_sep'].value
    k_plaque_inc, k_plaque_sep, k_ADCP, k_mAbcomplex_clear, k_mAb_transport = params['k_plaque_inc'].value, params['k_plaque_sep'].value, params['k_ADCP'].value, params['k_mAbcomplex_clear'].value, params['k_mAb_transport'].value
    k_mAb_transport_back, clearance, k_off_ma0, k_off_ma1, k_off_ma2, k_offPF = params['k_mAb_transport_back'].value, params['clearance'].value, params['k_off_ma0'].value, params['k_off_ma1'].value, params['k_off_ma2'].value, params['k_offPF'].value 
    # res = odeint(ode_model, initial_conditions, t,
    #                          args=(k_in, k_clear_Abeta, k_clear_olig, k_synth_FcR, 
    #                                k_clear_FcR, k_clear_P,
    #                                k_onPP, k_onPD, k_onPF, k_olig_inc, k_olig_sep, 
    #                                k_plaque_inc, k_plaque_sep, k_ADCP, k_mAbcomplex_clear, 
    #                                k_mAb_transport, k_mAb_transport_back, clearance, 
    #                                k_off_ma0, k_off_ma1, k_off_ma2, k_offPF, dose_list),
    #                                hmax=0.1)
    res = scipy.integrate.solve_ivp(fun=lambda t, y: ode_model(t, y, k_in, k_clear_Abeta, k_clear_olig, k_synth_FcR, 
                                   k_clear_FcR, k_clear_P,
                                   k_onPP, k_onPD, k_onPF, k_olig_inc, k_olig_sep, 
                                   k_plaque_inc, k_plaque_sep, k_ADCP, k_mAbcomplex_clear, 
                                   k_mAb_transport, k_mAb_transport_back, clearance, 
                                   k_off_ma0, k_off_ma1, k_off_ma2, k_offPF, dose_list),
                                             t_span=[t[0], t[-1]],
                                             y0=initial_conditions,
                                             t_eval=t)
    return res
#need to develop this function to directly compare plasma mAb to PK profile and to fit for the % reduction in plaque as can be measured from SUVr
def error(params, initial_conditions, tspan, data, dose_list):
    sol = ode_solver(tspan, initial_conditions, params, dose_list)
    plaque = sol.y[2] 
    percentage = [(((plaque[0]-plaque[i])/plaque[0])*100) for i in range(len(plaque))]
    plasma_mAb = sol.y[10]
    prediction = np.column_stack((percentage, plasma_mAb))
    return (plasma_mAb - data).ravel()

def dosefn(dose_list, t):
    infusion = (((((10 * 70)/1000/3.22)/147181.62))*1e9)/360 # nM
    f = 0.5
    delta = 0.1

    sol = 0
    for n in dose_list:
        if t >=(n-(0.5*360)) and t<=(n+(1.5*360)):
            sol =  ((infusion/2)/math.atan(1/delta))*(math.atan(math.sin(2*math.pi*(t-n-0.5)*f)/delta)) + infusion/2

    return sol

def per_iteration(params, iteration, resid, *args, **kws):
    print('ITER ', iteration, 'RESID ', sum(resid**2))

tspan = np.arange(0, (24*100*360), 360)
time = [i*(7*24*360) for i in [0, 79]]
fall = [0, 77.5]
df1 = pd.DataFrame({'Observed': fall,
                     'Time': time
                    })

slope, intercept, r_value, p_value, std_err = stats.linregress(df1['Time'], df1['Observed'])

line = slope*tspan+intercept
data1 = line
df2 = pd.read_csv('PK_fit_data.csv')
data2 = df2['PK'].values
times2 = df2['time'].values
data = np.column_stack((data1, data2))

# plt.plot(times2, data2)
# plt.show()

# plt.plot(times2, line)
# plt.show()

# need to decide what input concentrations to use for each species - could try with these and then also with the values used in the other model (much higher)
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
count = 0
result = minimize(error, params, args=(initial_conditions, tspan, data2, dose_list), method='leastsq', iter_cb=per_iteration)#, nan_policy='omit')


print(result.params)
print(report_fit(result))

plt.figure(2)
plt.plot(times2, data2, color='r', linestyle='dashed')
fit = ode_solver(tspan, initial_conditions, result.params, dose_list)
plt.plot(tspan, fit[:,10], color='b')
plt.show()
