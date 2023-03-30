import pandas as pd
import matplotlib.pyplot as plt

results = pd.read_csv("output/biweekly_500_param_ab.csv")
results_m = pd.read_csv("output/biweekly_500_m_param_ab.csv")

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

plt.plot(time, mean_SUVr_m, color = 'blue', label='10 mg/kg bi-monthly')
plt.fill_between(time, fifth_SUVr_m, 
                 ninety_fifth_m, 
                 color='blue', alpha=0.2, label = 'percentile 95% CI')

plt.xlabel("Time, days")
plt.ylabel("ab")
plt.legend()
#plt.axhline(y=1.17, linestyle='dashed', color = 'black')
plt.suptitle("Change in ab over 18 months treatment", y=1.05, fontsize=18)
plt.title("ab profiles for 500 PD parameter samples", fontsize=10)
#plt.show()
plt.savefig('plots/ab_biweekly_monthly_500_param.png')