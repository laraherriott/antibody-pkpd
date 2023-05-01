import pandas as pd
import matplotlib.pyplot as plt

results = pd.read_csv('output/15y_replication/biweekly_1_15y_off_param_tau.csv')
results_m = pd.read_csv('output/15y_replication/biweekly_1_15ycont_param_tau.csv')
results_old0 = pd.read_csv('output/tau_regimens/mid1_biweekly_1_15y_param_tau.csv')
results_old1 = pd.read_csv('output/tau_regimens/mid2_biweekly_1_15y_param_tau.csv')
results_2old0 = pd.read_csv('output/tau_regimens/mid3_biweekly_1_15y_param_tau.csv')
results_2old1 = pd.read_csv('output/tau_regimens/mid4_biweekly_1_15y_param_tau.csv')

results_0 = pd.read_csv('output/tau_regimens/mid5_biweekly_1_15y_param_tau.csv')
results_1 = pd.read_csv('output/tau_regimens/mid6_biweekly_1_15y_param_tau.csv')

# Visualisation

time = results["time"]
mean_SUVr = results['Average']
fifth_SUVr = results['fifth']
ninety_fifth = results['ninety-fifth']

mean_SUVr_m = results_m['Average']
fifth_SUVr_m = results_m['fifth']
ninety_fifth_m = results_m['ninety-fifth']

mean_SUVr_old0 = results_old0['Average']
fifth_SUVr_old0 = results_old0['fifth']
ninety_fifth_old0 = results_old0['ninety-fifth']

mean_SUVr_old1 = results_old1['Average']
fifth_SUVr_old1 = results_old1['fifth']
ninety_fifth_old1 = results_old1['ninety-fifth']

mean_SUVr_2old0 = results_2old0['Average']
fifth_SUVr_2old0 = results_2old0['fifth']
ninety_fifth_2old0 = results_2old0['ninety-fifth']

mean_SUVr_2old1 = results_2old1['Average']
fifth_SUVr_2old1 = results_2old1['fifth']
ninety_fifth_2old1 = results_2old1['ninety-fifth']

mean_SUVr_0 = results_0['Average']
fifth_SUVr_0 = results_0['fifth']
ninety_fifth_0 = results_0['ninety-fifth']

mean_SUVr_1 = results_1['Average']
fifth_SUVr_1 = results_1['fifth']
ninety_fifth_1 = results_1['ninety-fifth']

fig, ax = plt.subplots()



ax.plot(time, mean_SUVr_m, color = 'black', linestyle='dashed', label='10 mg/kg bi-weekly continuous')

ax.plot(time, mean_SUVr_old0, color = '#1c9099', label='+ 10 mg/kg monthly')
ax.plot(time, mean_SUVr_old1, color = '#43a2ca', label='+ 10 mg/kg bi-monthly')
ax.plot(time, mean_SUVr_2old0, color = '#a8ddb5', label='+ 10 mg/kg quarterly')
ax.plot(time, mean_SUVr_2old1, color = '#feb24c', label='+ 5 mg/kg monthly')
ax.plot(time, mean_SUVr_0, color = '#f03b20', label='+ 5 mg/kg bi-monthly')
ax.plot(time, mean_SUVr_1, color = '#ce56d2', label='+ 2 mg/kg monthly')

ax.plot(time, mean_SUVr, color = 'black', label='10 mg/kg bi-weekly 18 months')

ax.set_xlabel("Time, years")
#ax.set_ylabel("Plasma Aβ42/40 Ratio")
ax.set_ylabel("Plasma p-tau181 (pg/mL)")
#ax.set_ylabel("SUVr")
no_years = 15
xticks = [i*(365*24) for i in range(no_years+1)]
xtick_labels = [i for i in range(no_years+1)]
ax.set_xticks(xticks, xtick_labels)
#ax.axhline(y=1.17, linestyle='dashed', color = '#808080')
ax.axvline(x=(24*(365*1.5)), linestyle='dashed', color = '#808080')
#plt.suptitle("Change in {} over 15 years".format(model_type), y=1.05, fontsize=18)
#ax.set_title("Plasma Aβ42/40 Ratio profiles for different treatment regimens")
ax.set_title("Plasma p-tau181 (pg/mL) profiles for different treatment regimens")
#ax.set_title("SUVr profiles for different treatment regimens")

# pos = ax.get_position()
# ax.set_position([pos.x0, pos.y0, pos.width*0.5, pos.height])
# ax.legend(loc='center right', bbox_to_anchor=(2.25, 0.5))
#plt.show()
plt.savefig('plots/tau_biweekly_cont_vs_off_6mid_param.png')
