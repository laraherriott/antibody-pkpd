#
# example flow - here sampling for 1000 different patients with covariates, inter-individual variability (and, to be added, residual variability)
#

import pandas as pd
import numpy as np
import random
import matplotlib.pyplot as plt

from dose_switch import LecanemabModel
from solution import Solution

#random.seed(1)

n = 1
model_type = 'suvr'

doses = [1+x for x in range(20)]

for dose in doses:

    results = pd.DataFrame()
    dict_of_cols = {}

    model = LecanemabModel(model_type, 10, (14*24), (24*(365*15)), dose_change=[dose, (28*24), (12*1080)], median_patient=True, param_type='noiiv')
    model()


    # order: baseline SUVr, Kin, Emax, EC50; baseline AB, Kout, slope; baseline ptau, Kout, slope
    sample = [1.34, 0.232, 1.54, 75.0, 0.0842, 0.367, 0.00155, 4.06, 0.468, 0.00313]

    for i in range(n):
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
        # Solve model

        solver = Solution(model, 0, (24*(365*15)), 1)

        # Save values

        solutions = solver.solve()

        if i == 0:
            time = solutions.t
            results['time'] = time

        biomarker = solutions.y[2]
        central_L = solutions.y[0]
        peripheral_L = solutions.y[1]

        dict_of_cols['{}{}'.format(model_type, i)] = biomarker


    # Process data
    results = pd.concat([results, pd.DataFrame(dict_of_cols)], axis=1)

    total_df = \
        results[list(results.filter(regex=model_type))]
    results["Average"] = total_df.median(axis=1)
    results["fifth"] = total_df.quantile(q=0.05, axis=1)
    results["ninety-fifth"] = total_df.quantile(q=0.95, axis=1)

    results.to_csv("output/15y_monthly/{}mg_{}_15ycont_param_{}.csv".format(dose, model_type, n), index=False)

# results = pd.read_csv('output/biweekly_1_15y_off_param_suvr.csv')
# results_m = pd.read_csv('output/biweekly_1_15ycont_param_suvr.csv')
# results_old0 = pd.read_csv('output/mid1_biweekly_1_15y_param_suvr.csv')
# results_old1 = pd.read_csv('output/mid2_biweekly_1_15y_param_suvr.csv')
# results_2old0 = pd.read_csv('output/mid3_biweekly_1_15y_param_suvr.csv')
# results_2old1 = pd.read_csv('output/mid4_biweekly_1_15y_param_suvr.csv')

# # Visualisation

# time = results["time"]
# mean_SUVr = results['Average']
# fifth_SUVr = results['fifth']
# ninety_fifth = results['ninety-fifth']

# mean_SUVr_m = results_m['Average']
# fifth_SUVr_m = results_m['fifth']
# ninety_fifth_m = results_m['ninety-fifth']

# mean_SUVr_old0 = results_old0['Average']
# fifth_SUVr_old0 = results_old0['fifth']
# ninety_fifth_old0 = results_old0['ninety-fifth']

# mean_SUVr_old1 = results_old1['Average']
# fifth_SUVr_old1 = results_old1['fifth']
# ninety_fifth_old1 = results_old1['ninety-fifth']

# mean_SUVr_2old0 = results_2old0['Average']
# fifth_SUVr_2old0 = results_2old0['fifth']
# ninety_fifth_2old0 = results_2old0['ninety-fifth']

# mean_SUVr_2old1 = results_2old1['Average']
# fifth_SUVr_2old1 = results_2old1['fifth']
# ninety_fifth_2old1 = results_2old1['ninety-fifth']

# mean_SUVr_0 = results_0['Average']
# fifth_SUVr_0 = results_0['fifth']
# ninety_fifth_0 = results_0['ninety-fifth']

# mean_SUVr_1 = results_1['Average']
# fifth_SUVr_1 = results_1['fifth']
# ninety_fifth_1 = results_1['ninety-fifth']

# fig, ax = plt.subplots()

# ax.plot(time, mean_SUVr, color = 'black', label='10 mg/kg bi-weekly 18 months')

# ax.plot(time, mean_SUVr_m, color = 'black', linestyle='dashed', label='10 mg/kg bi-weekly continuous')

# ax.plot(time, mean_SUVr_old0, color = '#1c9099', label='+ 10 mg/kg monthly')
# ax.plot(time, mean_SUVr_old1, color = '#feb24c', label='+ 5 mg/kg monthly')
# ax.plot(time, mean_SUVr_2old0, color = '#ce56d2', label='+ 2 mg/kg monthly')
# ax.plot(time, mean_SUVr_2old1, color = '#f03b20', label='+ 5 mg/kg bi-monthly')
# ax.plot(time, mean_SUVr_0, color = '#43a2ca', label='+ 10 mg/kg bi-monthly')
# ax.plot(time, mean_SUVr_1, color = '#a8ddb5', label='+ 10 mg/kg quarterly')


# ax.xlabel("Time, years")
# ax.ylabel(model_type)
# no_years = 15
# xticks = [i*(365*24) for i in range(no_years+1)]
# xtick_labels = [i for i in range(no_years+1)]
# ax.xticks(xticks, xtick_labels)
# if model_type == 'suvr':
#     #plt.ylim((1, 1.4))
#     ax.axhline(y=1.17, linestyle='dashed', color = '#808080')
# ax.axvline(x=(24*(365*1.5)), linestyle='dashed', color = '#808080')
# #plt.suptitle("Change in {} over 15 years".format(model_type), y=1.05, fontsize=18)
# ax.title("{} profiles for 18 month vs continuous treatment".format(model_type), fontsize=10)

# pos = ax.get_position()
# ax.set_position([pos.x0, pos.y0, pos.width*0.9, pos.height])
# ax.legend(loc='center right', bbox_to_anchor=(1.25, 0.5))
# #plt.show()
# plt.savefig('plots/{}_biweekly_cont_vs_off_4mid_{}_param.png'.format(model_type, n))

