#
# example flow - here sampling for 1000 different patients with covariates, inter-individual variability (and, to be added, residual variability)
#

import pandas as pd
import random
import matplotlib.pyplot as plt

from lecanemab_model import LecanemabModel
from solution import Solution

random.seed(1)
n = 2

# bi-weekly
results = pd.DataFrame()
metadata = pd.DataFrame(columns = ['age', 'weight', 'sex', 'ada', 'apoe4', 'japanese', 'albumin', 'baseline_suvr'])
dict_of_cols = {}

# monthly
results_m = pd.DataFrame()
metadata_m = pd.DataFrame(columns = ['age', 'weight', 'sex', 'ada', 'apoe4', 'japanese', 'albumin', 'baseline_suvr'])
dict_of_cols_m = {}


for i in range (n):
    # Set up model - bi-weekly dosing

    model_type = 'suvr'

    model = LecanemabModel(model_type, 10, 3.5, 540)
    model()

    metadata.loc[len(metadata)] = [model.patient.age_SUVr, model.patient.weight_SUVr, model.patient.sex, model.patient.ADA,
                                   model.patient.APOE, model.patient.race, model.patient.albumin, model.baseline_SUVr]

    # Solve model

    solver = Solution(model, 0, 540, 1) # 12960 hours = 18 months; 1080 half days

    # Save values

    solutions = solver.solve()

    if i == 0:
        time = solutions.t
        results['time'] = time

    SUVr = solutions.y[2]
    central_L = solutions.y[0]
    peripheral_L = solutions.y[1]

    dict_of_cols['SUVr_{}'.format(i)] = SUVr

    # Set up model - monthly dosing

    model_m = LecanemabModel(model_type, 10, 30, 540)
    model_m()

    metadata_m.loc[len(metadata_m)] = [model_m.patient.age_SUVr, model_m.patient.weight_SUVr, model_m.patient.sex, model_m.patient.ADA,
                                       model_m.patient.APOE, model_m.patient.race, model_m.patient.albumin, model_m.baseline_SUVr]

    # Solve model

    solver_m = Solution(model_m, 0, 540, 1) # 12960 hours = 18 months; 1080 half days

    # Save values

    solutions_m = solver_m.solve()

    if i == 0:
        time = solutions_m.t
        results_m['time'] = time

    SUVr_m = solutions_m.y[2]
    central_L_m = solutions_m.y[0]
    peripheral_L_m = solutions_m.y[1]

    dict_of_cols_m['SUVr_{}'.format(i)] = SUVr_m


# Process data

results = pd.concat([results, pd.DataFrame(dict_of_cols)], axis=1)
results_m = pd.concat([results_m, pd.DataFrame(dict_of_cols_m)], axis=1)

total_df = \
    results[list(results.filter(regex='SUVr'))]
results["Average"] = total_df.median(axis=1)
results["Std"] = total_df.std(axis=1)

results.to_csv("output/biweekly_1000_individuals.csv", index=False)
metadata.to_csv("output/patient_characteristics.csv")

total_df_m = \
    results_m[list(results_m.filter(regex='SUVr'))]
results_m["Average"] = total_df_m.median(axis=1)
results_m["Std"] = total_df_m.std(axis=1)

results_m.to_csv("output/biweekly_1000_m_individuals.csv", index=False)
metadata_m.to_csv("output/patient_characteristics_m.csv")

# Visualisation

time = results["time"]
mean_SUVr = results['Average']
sd_SUVr = results['Std']

mean_SUVr_m = results_m['Average']
sd_SUVr_m = results_m['Std']

plt.plot(time, mean_SUVr, color = 'tab:pink', label='10 mg/kg bi-weekly')
plt.fill_between(time, mean_SUVr - (2*sd_SUVr), 
                 mean_SUVr + (2*sd_SUVr), 
                 color='tab:pink', alpha=0.2, label = "2 s.d.")

plt.plot(time, mean_SUVr_m, color = 'blue', label='10 mg/kg monthly')
plt.fill_between(time, mean_SUVr_m - (2*sd_SUVr_m), 
                 mean_SUVr_m + (2*sd_SUVr_m), 
                 color='blue', alpha=0.2, label = '2 s.d.')

plt.xlabel("Time, days")
plt.ylabel("SUVr")
plt.legend()
plt.axhline(y=1.17, linestyle='dashed', color = 'black')
plt.suptitle("Change in SUVr over 18 months treatment", y=1.05, fontsize=18)
plt.title("Profiles for 1000 individuals", fontsize=10)
plt.ylim((0.8, 1.8))
#plt.show()
plt.savefig('plots/SUVr_biweekly_monthly_1000_individuals.png')

