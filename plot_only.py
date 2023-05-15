import pandas as pd
import matplotlib.pyplot as plt

results = pd.read_csv("final_data/biweekly_100_param_suvr.csv")
results_m = pd.read_csv("final_data/biweekly_100_m_param_suvr.csv")

real = pd.read_csv('final_data/suvr_biweekly_vs_monthly_digitized.csv')

suvr_bw = real[real['type']=='suvr-biweekly']
suvr_m = real[real['type']=='suvr-monthly']
tau_bw = real[real['type']=='tau-biweekly']
tau_m = real[real['type']=='tau-monthly']
ab_bw = real[real['type']=='ab-biweekly']
ab_m = real[real['type']=='ab-monthly']

time = results["time"]
mean_SUVr = results['Average']
fifth_SUVr = results['fifth']
ninety_fifth = results['ninety-fifth']

mean_SUVr_m = results_m['Average']
fifth_SUVr_m = results_m['fifth']
ninety_fifth_m = results_m['ninety-fifth']

time_bars = suvr_bw['time']

times = []
for n in range(len(time_bars)):
    times.append(time_bars[n]*(28*24))
length = len(time_bars)

print(suvr_m['lowerCI'])

min_bw = suvr_bw['lowerCI']
min_m = suvr_m['lowerCI']


max_bw = suvr_bw['upperCI']
max_m = suvr_m['upperCI']
print(max_bw[0])


plt.plot(time, mean_SUVr, color = 'tab:pink', label='10 mg/kg bi-weekly')
plt.fill_between(time, fifth_SUVr, 
                 ninety_fifth, 
                 color='tab:pink', alpha=0.2, label = 'percentile 95% CI')

plt.plot(time, mean_SUVr_m, color = 'blue', label='10 mg/kg monthly')
plt.fill_between(time, fifth_SUVr_m, 
                 ninety_fifth_m, 
                 color='blue', alpha=0.2, label = 'percentile 95% CI')
plt.plot(times, suvr_bw['mean'], color = 'tab:pink', marker='o', linestyle='None')
#plt.plot(times, suvr_m['mean'], color = 'blue', marker='o', linestyle='None')
for i in range(length):
    plt.vlines(times[i], min_bw[i], max_bw[i], color = 'tab:pink')
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
plt.ylim(1, 1.4)
#plt.suptitle("Change in SUVr over 18 months treatment", y=1.05, fontsize=18)
plt.title("Change in SUVr over 18 month treatment", fontsize=10)
#plt.show()
plt.savefig('final_data/suvr_biweekly_monthly_100.png')
