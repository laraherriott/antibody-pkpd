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
model_type = 'tau'


results_0 = pd.DataFrame()
dict_of_cols_0 = {}
results_1 = pd.DataFrame()
dict_of_cols_1 = {}
# results_2 = pd.DataFrame()
# dict_of_cols_2 = {}
# results_3 = pd.DataFrame()
# dict_of_cols_3 = {}
# results_4 = pd.DataFrame()
# dict_of_cols_4 = {}
# results_5 = pd.DataFrame()
# dict_of_cols_5 = {}

model_0 = LecanemabModel(model_type, 10, (14*24), (24*(365*15)), dose_change=[10, (14*24), (12*1080)], median_patient=True, param_type='noiiv')
model_1 = LecanemabModel(model_type, 10, (14*24), (24*(365*15)), dose_change=[0, (28*2*24), (12*1080)], median_patient=True, param_type='noiiv')

model_0()
model_1()

# model_2 = LecanemabModel(model_type, 10, (14*24), (24*(365*15)), dose_change=[10, (28*3*24), (12*1080)], median_patient=True, param_type='noiiv')
# model_3 = LecanemabModel(model_type, 10, (14*24), (24*(365*15)), dose_change=[5, (28*24), (12*1080)], median_patient=True, param_type='noiiv')

# model_2()
# model_3()

# model_4 = LecanemabModel(model_type, 10, (14*24), (24*(365*15)), dose_change=[5, (28*2*24), (12*1080)], median_patient=True, param_type='noiiv')
# model_5 = LecanemabModel(model_type, 10, (14*24), (24*(365*15)), dose_change=[2, (28*24), (12*1080)], median_patient=True, param_type='noiiv')

# model_4()
# model_5()

# order: baseline SUVr, Kin, Emax, EC50; baseline AB, Kout, slope; baseline ptau, Kout, slope
sample = [1.34, 0.232, 1.54, 75.0, 0.0842, 0.367, 0.00155, 4.06, 0.468, 0.00313]

for i in range(n):
    models = [model_0, model_1]#, model_2, model_3, model_4, model_5]
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
    # Solve model

    solver_0 = Solution(model_0, 0, (24*(365*15)), 1)
    solver_1 = Solution(model_1, 0, (24*(365*15)), 1)
    # solver_2 = Solution(model_2, 0, (24*(365*15)), 1)
    # solver_3 = Solution(model_3, 0, (24*(365*15)), 1)
    # solver_4 = Solution(model_4, 0, (24*(365*15)), 1)
    # solver_5 = Solution(model_5, 0, (24*(365*15)), 1)
    # Save values


    solutions_0 = solver_0.solve()

    if i == 0:
        time = solutions_0.t
        results_0['time'] = time

    biomarker_0 = solutions_0.y[2]
    central_L_0 = solutions_0.y[0]
    peripheral_L_0 = solutions_0.y[1]

    dict_of_cols_0['{}{}'.format(model_type, i)] = biomarker_0

    solutions_1 = solver_1.solve()

    if i == 0:
        time = solutions_1.t
        results_1['time'] = time

    biomarker_1 = solutions_1.y[2]
    central_L_1 = solutions_1.y[0]
    peripheral_L_1 = solutions_1.y[1]

    dict_of_cols_1['{}{}'.format(model_type, i)] = biomarker_1

    # solutions_2 = solver_2.solve()

    # if i == 0:
    #     time = solutions_2.t
    #     results_2['time'] = time

    # biomarker_2 = solutions_2.y[2]
    # central_L_2 = solutions_2.y[0]
    # peripheral_L_2 = solutions_2.y[1]

    # dict_of_cols_2['{}{}'.format(model_type, i)] = biomarker_2

    # solutions_3 = solver_3.solve()

    # if i == 0:
    #     time = solutions_3.t
    #     results_3['time'] = time

    # biomarker_3 = solutions_3.y[2]
    # central_L_3 = solutions_3.y[0]
    # peripheral_L_3 = solutions_3.y[1]

    # dict_of_cols_3['{}{}'.format(model_type, i)] = biomarker_3

    # solutions_4 = solver_4.solve()

    # if i == 0:
    #     time = solutions_4.t
    #     results_4['time'] = time

    # biomarker_4 = solutions_4.y[2]
    # central_L_4 = solutions_4.y[0]
    # peripheral_L_4 = solutions_4.y[1]

    # dict_of_cols_4['{}{}'.format(model_type, i)] = biomarker_4

    # solutions_5 = solver_5.solve()

    # if i == 0:
    #     time = solutions_5.t
    #     results_5['time'] = time

    # biomarker_5 = solutions_5.y[2]
    # central_L_5 = solutions_5.y[0]
    # peripheral_L_5 = solutions_5.y[1]

    # dict_of_cols_5['{}{}'.format(model_type, i)] = biomarker_5

# Process data
results_0 = pd.concat([results_0, pd.DataFrame(dict_of_cols_0)], axis=1)
results_1 = pd.concat([results_1, pd.DataFrame(dict_of_cols_1)], axis=1)

# results_2 = pd.concat([results_2, pd.DataFrame(dict_of_cols_2)], axis=1)
# results_3 = pd.concat([results_3, pd.DataFrame(dict_of_cols_3)], axis=1)

# results_4 = pd.concat([results_4, pd.DataFrame(dict_of_cols_4)], axis=1)
# results_5 = pd.concat([results_5, pd.DataFrame(dict_of_cols_5)], axis=1)

total_df_0 = \
    results_0[list(results_0.filter(regex=model_type))]
results_0["Average"] = total_df_0.median(axis=1)
results_0["fifth"] = total_df_0.quantile(q=0.05, axis=1)
results_0["ninety-fifth"] = total_df_0.quantile(q=0.95, axis=1)

results_0.to_csv("output/15y_replication/biweekly_{}_15ycont_param_{}.csv".format(n, model_type), index=False)

total_df_1 = \
    results_1[list(results_1.filter(regex=model_type))]
results_1["Average"] = total_df_1.median(axis=1)
results_1["fifth"] = total_df_1.quantile(q=0.05, axis=1)
results_1["ninety-fifth"] = total_df_1.quantile(q=0.95, axis=1)

results_1.to_csv("output/15y_replication/biweekly_{}_15y_off_param_{}.csv".format(n, model_type), index=False)

# total_df_2 = \
#     results_2[list(results_2.filter(regex=model_type))]
# results_2["Average"] = total_df_2.median(axis=1)
# results_2["fifth"] = total_df_2.quantile(q=0.05, axis=1)
# results_2["ninety-fifth"] = total_df_2.quantile(q=0.95, axis=1)

# results_2.to_csv("output/mid3_biweekly_{}_15y_param_{}.csv".format(n, model_type), index=False)

# total_df_3 = \
#     results_3[list(results_3.filter(regex=model_type))]
# results_3["Average"] = total_df_3.median(axis=1)
# results_3["fifth"] = total_df_3.quantile(q=0.05, axis=1)
# results_3["ninety-fifth"] = total_df_3.quantile(q=0.95, axis=1)

# results_3.to_csv("output/mid4_biweekly_{}_15y_param_{}.csv".format(n, model_type), index=False)

# total_df_4 = \
#     results_4[list(results_4.filter(regex=model_type))]
# results_4["Average"] = total_df_4.median(axis=1)
# results_4["fifth"] = total_df_4.quantile(q=0.05, axis=1)
# results_4["ninety-fifth"] = total_df_4.quantile(q=0.95, axis=1)

# results_4.to_csv("output/mid5_biweekly_{}_15y_param_{}.csv".format(n, model_type), index=False)

# total_df_5 = \
#     results_5[list(results_5.filter(regex=model_type))]
# results_5["Average"] = total_df_5.median(axis=1)
# results_5["fifth"] = total_df_5.quantile(q=0.05, axis=1)
# results_5["ninety-fifth"] = total_df_5.quantile(q=0.95, axis=1)

# results_5.to_csv("output/mid6_biweekly_{}_15y_param_{}.csv".format(n, model_type), index=False)

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

