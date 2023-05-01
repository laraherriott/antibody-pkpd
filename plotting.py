import pandas as pd
import matplotlib.pyplot as plt

results = pd.read_csv("output/biweekly_100_param_tau.csv")
results_m = pd.read_csv("output/biweekly_100_m_param_tau.csv")
n = 100
model_type = 'Plasma p-tau181 (pg/mL)'

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

plt.xlabel("Time, months")
plt.ylabel(model_type)
no_months = 18
xticks = [i*(28*24) for i in range(0, no_months+1, 3)]
xtick_labels = [i for i in range(0, no_months+1, 3)]
plt.xticks(xticks, xtick_labels)
plt.legend()
if model_type == 'SUVr':
    plt.ylim((1, 1.4))
    plt.axhline(y=1.17, linestyle='dashed', color = 'black')
    plt.axvline(x=(24*(28*15)), linestyle='dashed', color = 'red')
#plt.suptitle("Change in {} over 15 years".format(model_type), y=1.05, fontsize=18)
plt.title("Change in {} over 18 month treatment".format(model_type), fontsize=10)

#plt.show()
plt.savefig('plots/tau_biweekly_monthly_{}_param.png'.format(n))