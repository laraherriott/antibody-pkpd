import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy
from lmfit import minimize, Parameters, report_fit
from scipy.integrate import odeint

clearance = 0.0181

dose_list = []
i = 0
while i <= (24*100):
    dose_list.append(int(i))
    i += (14*24)

params = Parameters()
params.add('k_in', value=k_in, min=0, max=1)
params.add('k_synth_FcR', value=k_synth_FcR, min=0, max=1)
params.add('k_clear_FcR', value=k_clear_FcR, min=0, max=1)
params.add('k_clear_Abeta', value=k_clear_Abeta, min=0, max=1)
params.add('k_clear_olig', value=k_clear_olig, min=0, max=1)
params.add('k_clear_P', value=k_clear_P, min=0, max=1)
params.add('clearance', value=clearance, min=0, max=1)

def ode_model(y, t, clearance, dose_list):
    mAb_plasma = y[0]
    dmAb_plasma = dosefn(dose_list, t) - clearance*mAb_plasma

    dYdt = [dmAb_plasma]
    return dYdt

def ode_solver(t, initial_conditions, params, dose_list):
    #initmAb = initial_conditions
    clearance = params['clearance'].value
    res = odeint(ode_model, initial_conditions, t, args=(clearance, dose_list), hmax=0.1)
    return res

def error(params, initial_conditions, tspan, data, dose_list):
    sol = ode_solver(tspan, initial_conditions, params, dose_list)
    return (sol - data).ravel()

def dosefn(dose_list, t):
    infusion = (((((10 * 70)/1000/3.22)/147181.62))*1e9)/360 # nM
    f = 0.5
    delta = 0.1

    sol = 0
    for n in dose_list:
        if t >=(n-(360*0.5)) and t<=(n+(360*1.5)):
            sol =  ((infusion/2)/math.atan(1/delta))*(math.atan(math.sin(2*math.pi*(t-n-0.5)*f)/delta)) + infusion/2

    return sol

tspan = np.arange(0, (24*100*360), 360)

df2 = pd.read_csv('PK_fit_data.csv')
data2 = df2['PK'].values
times2 = df2['time'].values

initmAb = 0
initPlasmamAb = 0

result = minimize(error, params, args=([initmAb, initPlasmamAb], tspan, data2, dose_list), method='leastsq')


print(result.params)
print(report_fit(result))

plt.figure(1)
plt.plot(times2, data2)
plt.show()

fit = ode_solver(tspan, initmAb, result.params, dose_list)
plt.figure(3)
plt.plot(tspan, fit)
plt.show()