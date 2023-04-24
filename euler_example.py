#
# example flow - here sampling for 1000 different patients with covariates, inter-individual variability (and, to be added, residual variability)
#

import pandas as pd
import numpy as np
import random
import matplotlib.pyplot as plt

from lecanemab_model import LecanemabModel
from solution import Solution

#random.seed(1)
n = 1
model_type = 'suvr'

# bi-weekly
results = pd.DataFrame()
dict_of_cols = {}

# monthly
results_m = pd.DataFrame()
dict_of_cols_m = {}

model_bw = LecanemabModel(model_type, 10, (14*24), (12*1080), median_patient=True, param_type='noiiv')
model_bw()

model_m = LecanemabModel(model_type, 10, (24*28), (12*1080), median_patient=True, param_type='noiiv')
model_m()

# order: baseline SUVr, Kin, Emax, EC50; baseline AB, Kout, slope; baseline ptau, Kout, slope
means = [1.34, 0.232, 1.54, 75.0, 0.0842, 0.367, 0.00155, 4.06, 0.468, 0.00313]
rse = [0.873, 11.1, 11.8, 19.6, 4.28, 1.97, 9.32, 1.61, 20.7, 15.6]
se = [estimate*(error/100) for estimate,error in zip(means,rse)]
var = [value**2 for value in se]
cov = np.zeros((len(var), len(var)))
np.fill_diagonal(cov, var)

parameter_sample = np.random.multivariate_normal(means, cov, n)

for i in range(n):
    sample = parameter_sample[i]
    models = [model_m, model_bw]
    for model in models:
        # Overwrite parameters for this simulation
        # model.CL = 0.0181 * 24
        # model.V1 = 3.22
        # model.V2 = 2.19
        model.baseline_SUVr = 1.38 #sample[0]
        model.SUVr_Kin = sample[1] / (365*24)
        model.SUVr_Kout = model.SUVr_Kin/model.baseline_SUVr
        model.Emax = sample[2]
        model.SUVr_EC50 = sample[3]
        model.baseline_Abeta = 0.0842 #sample[4]
        model.Abeta_Kout = sample[5] / (365*24)
        model.Abeta_slope = sample[6]
        model.Abeta_Kin = model.Abeta_Kout * model.baseline_Abeta
        model.baseline_tau = sample[7]
        model.tau_Kout = sample[8] / (365*24)
        model.tau_slope = sample[9]
        model.tau_Kin = model.tau_Kout * model.baseline_tau
        # need to add the rest for the other parameters after have implemented those equations.

    # Solve model

    solver_bw = Solution(model_bw, 0, 12960, 1) # 12960 hours = 18 months; 1080 half days
    solver_m = Solution(model_m, 0, 12960, 1) # 12960 hours = 18 months; 1080 half days

    # Save values

    solutions, t_list = solver_bw.solve()

    if i == 0:
        time = t_list
        results['time'] = time

    biomarker = solutions[2]
    central_L = solutions[0]
    peripheral_L = solutions[1]

    dict_of_cols['{}{}'.format(model_type, i)] = biomarker

    solutions_m, t_list_m = solver_m.solve()

    if i == 0:
        time = t_list_m
        results_m['time'] = time

    biomarker_m = solutions_m[2]
    central_L_m = solutions_m[0]
    peripheral_L_m = solutions_m[1]

    dict_of_cols_m['{}{}'.format(model_type, i)] = biomarker_m


# Process data

results = pd.concat([results, pd.DataFrame(dict_of_cols)], axis=1)
results_m = pd.concat([results_m, pd.DataFrame(dict_of_cols_m)], axis=1)

total_df = \
    results[list(results.filter(regex=model_type))]
results["Average"] = total_df.median(axis=1)
results["fifth"] = total_df.quantile(q=0.05, axis=1)
results["ninety-fifth"] = total_df.quantile(q=0.95, axis=1)

results.to_csv("output/biweekly_{}_param_{}.csv".format(n, model_type), index=False)

total_df_m = \
    results_m[list(results_m.filter(regex=model_type))]
results_m["Average"] = total_df_m.median(axis=1)
results_m["fifth"] = total_df_m.quantile(q=0.05, axis=1)
results_m["ninety-fifth"] = total_df_m.quantile(q=0.95, axis=1)

results_m.to_csv("output/biweekly_{}_m_param_{}.csv".format(n, model_type), index=False)

# Visualisation

time = results["time"]
mean_SUVr = results['Average']
fifth_SUVr = results['fifth']
ninety_fifth = results['ninety-fifth']

mean_SUVr_m = results_m['Average']
fifth_SUVr_m = results_m['fifth']
ninety_fifth_m = results_m['ninety-fifth']

plt.plot(time, mean_SUVr, color = 'tab:pink', label='10 mg/kg bi-weekly')
plt.fill_between(time, fifth_SUVr, 
                 ninety_fifth, 
                 color='tab:pink', alpha=0.2, label = 'percentile 95% CI')

plt.plot(time, mean_SUVr_m, color = 'blue', label='10 mg/kg monthly')
plt.fill_between(time, fifth_SUVr_m, 
                 ninety_fifth_m, 
                 color='blue', alpha=0.2, label = 'percentile 95% CI')

plt.xlabel("Time, hours")
plt.ylabel(model_type)
plt.legend()
if model_type == 'suvr':
    plt.ylim((1, 1.4))
    plt.axhline(y=1.17, linestyle='dashed', color = 'black')
    plt.axvline(x = 10950, linestyle='dashed', color = 'red')
plt.suptitle("Change in {} over 18 months treatment".format(model_type), y=1.05, fontsize=18)
plt.title("{} profiles for {} PD parameter samples".format(model_type, n), fontsize=10)

#plt.show()
plt.savefig('plots/{}_biweekly_monthly_{}_param.png'.format(model_type, n))
