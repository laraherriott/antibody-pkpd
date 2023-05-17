import pandas as pd
import matplotlib.pyplot as plt

results = pd.read_csv("final_data/biweekly_100_param_suvr.csv")
results_m = pd.read_csv("final_data/biweekly_100_m_param_suvr.csv")
results_5bw = pd.read_csv("final_data/5biweekly_100_param_suvr.csv")
results_5m = pd.read_csv("final_data/5monthly_100_param_suvr.csv")
results_25 = pd.read_csv("final_data/2_5100_param_suvr.csv")
#real = pd.read_csv('final_data/suvr_biweekly_vs_monthly_digitized.csv')

# suvr_bw = real[real['type']=='suvr-biweekly']
# suvr_m = real[real['type']=='suvr-monthly']
# tau_bw = real[real['type']=='tau-biweekly']
# tau_m = real[real['type']=='tau-monthly']
# ab_bw = real[real['type']=='ab-biweekly']
# ab_m = real[real['type']=='ab-monthly']

time = results["time"]
mean_SUVr = results['Average']
fifth_SUVr = results['fifth']
ninety_fifth = results['ninety-fifth']

mean_SUVr_m = results_m['Average']
fifth_SUVr_m = results_m['fifth']
ninety_fifth_m = results_m['ninety-fifth']

mean_SUVr_5bw = results_5bw['Average']
fifth_SUVr_5bw = results_5bw['fifth']
ninety_fifth_5bw = results_5bw['ninety-fifth']

mean_SUVr_5m = results_5m['Average']
fifth_SUVr_5m = results_5m['fifth']
ninety_fifth_5m = results_5m['ninety-fifth']

mean_SUVr_25 = results_25['Average']
fifth_SUVr_25 = results_25['fifth']
ninety_fifth_25 = results_25['ninety-fifth']

#time_bars = suvr_bw['time']

# times = []
# for n in range(len(time_bars)):
#     times.append(time_bars[n]*(28*24))
# length = len(time_bars)

# min_bw = tau_bw['lowerCI'].values
# min_m = tau_m['lowerCI'].values


# max_bw = tau_bw['upperCI'].values
# max_m = tau_m['upperCI'].values


plt.plot(time, mean_SUVr, color = 'tab:pink', label='10 mg/kg bi-weekly')
# plt.fill_between(time, fifth_SUVr, 
#                  ninety_fifth, 
#                  color='tab:pink', alpha=0.2, label = 'percentile 95% CI')
plt.plot(time, mean_SUVr_m, color = 'blue', label='10 mg/kg monthly')
# plt.fill_between(time, fifth_SUVr_m, 
#                  ninety_fifth_m, 
#                  color='blue', alpha=0.2, label = 'percentile 95% CI')
plt.plot(time, mean_SUVr_5bw, color = 'yellow', label='5 mg/kg bi-weekly')
# plt.fill_between(time, fifth_SUVr_5bw, 
#                  ninety_fifth_5bw, 
#                  color='yellow', alpha=0.4, label = 'percentile 95% CI')
plt.plot(time, mean_SUVr_5m, color = 'orange', label='5 mg/kg monthly')
# plt.fill_between(time, fifth_SUVr_5m, 
#                  ninety_fifth_5m, 
#                  color='orange', alpha=0.2, label = 'percentile 95% CI')
plt.plot(time, mean_SUVr_25, color = 'red', label='2.5 mg/kg bi-weekly')
# plt.fill_between(time, fifth_SUVr_25, 
#                  ninety_fifth_25, 
#                  color='red', alpha=0.2, label = 'percentile 95% CI')
# plt.plot(times, tau_bw['mean'], color = 'tab:pink', marker='o', linestyle='None')
# plt.plot(times, tau_m['mean'], color = 'blue', marker='o', linestyle='None')
# for i in range(length):
#     plt.vlines(times[i], min_bw[i], max_bw[i], color = 'tab:pink')
# for k in range(length):
#     plt.vlines(times[k], min_m[k], max_m[k], color = 'blue')

plt.xlabel("Time, months")
plt.ylabel("SUVr")
no_months = 18
xticks = [i*(28*24) for i in range(0, no_months+1, 3)]
xtick_labels = [i for i in range(0, no_months+1, 3)]
plt.xticks(xticks, xtick_labels)
plt.legend()
plt.axhline(y=1.17, linestyle='dashed', color = 'black')
plt.ylim(1.0, 1.4)
#plt.suptitle("Change in SUVr over 18 months treatment", y=1.05, fontsize=18)
plt.title("Change in SUVr over 18 month treatment", fontsize=10)
#plt.show()
plt.savefig('final_data/SUVr_all_noerrorbars_biweekly_monthly_100.png', dpi=300)
