#
# example flow - for 2 antibodies
#

import random
import pandas as pd
import scipy
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt

from aducanumab_model import MITModel
from antibodies_model import AbModel
from solution import Solution

# data = pd.read_csv('data/covid_ab_level_profile.csv')

# over_80 = data[data['age'] == 80]

# over_80_single_no_prev = over_80[over_80['model'] == 'pf_one_negative']

# spline = CubicSpline(over_80_single_no_prev['time'], over_80_single_no_prev['prediction'])

random.seed(1)

#pathway = MITModel('adu_path')
model = MITModel('brain')
#model = AbModel()
#model()

# solver2 = Solution(pathway, 0, 31536000, 1)
solver = Solution(model, 0, int(24*360*364*1.5), 360)

solutions= solver.solve()
# solutions2 = solver2.solve()

monomer = solutions.y[0]
oligomer = solutions.y[1]
plaque = solutions.y[2]
receptor = solutions.y[3]
antibody = solutions.y[4]
mon_ab = solutions.y[5]
olig_ab = solutions.y[6]
pl_ab = solutions.y[7]
olig_ab_fcr = solutions.y[8]
pl_ab_fcr = solutions.y[9]
lec_plasma = solutions.y[10]
consumed = solutions.y[11]
new_agg = solutions.y[12]
time = solutions.t

results = pd.DataFrame()
# results['time'] = time
# results['PK'] = lec_plasma
#results['Drug'] = lec_plasma
# results['Ab'] = antibody
# results['M'] = monomer
# results['O'] = oligomer
# results['P'] = plaque
# results['FcR'] = receptor
# results['M_Ab'] = mon_ab
# results['O_Ab'] = olig_ab
# results['P_Ab'] = pl_ab
# results['O_Ab_FcR'] = olig_ab_fcr
# results['P_Ab_FcR'] = pl_ab_fcr



#results.to_csv('PK_results.csv')
total_ab = []
per_ab = []
plaque_total = []
olig_total = []
mon_total = []
for i in range(len(antibody)):
    add = antibody[i] + mon_ab[i] + olig_ab[i] + pl_ab[i] + olig_ab_fcr[i] + pl_ab_fcr[i]
    total_ab.append(add)
    percent = (add/lec_plasma[i])*100
    per_ab.append(percent)
    plaque_add = plaque[i] + pl_ab[i] + pl_ab_fcr[i]
    plaque_total.append(plaque_add)
    olig_add = oligomer[i] + olig_ab[i] + olig_ab_fcr[i]
    olig_total.append(olig_add)
    mon_add = monomer[i] + mon_ab[i]
    mon_total.append(mon_add)

print(plaque_total[0], plaque_total[-1])
per_dec = ((plaque[0]-plaque[-1])/plaque[0])*100
print(per_dec)

# ab_1 = []
# ab_2= []

# for i in time:
#     ab = ((spline(i) * (7.375*1000))/1000)
#     ab_1.append(ab*0.5)
#     ab_2.append(ab*0.5)


fig1 = plt.figure(1)
plt.plot(time, plaque, label = 'Plaque (free)')
plt.plot(time, oligomer, label = 'Oligomer (free)')
plt.plot(time, monomer, label = 'Monomer (free)')
no_months = 18
xticks = [i*(28*24*360) for i in range(0, no_months+1, 3)]
xtick_labels = [i for i in range(0, no_months+1, 3)]
plt.xticks(xticks, xtick_labels)
plt.legend()
plt.xlabel("Time, months")
plt.ylabel("Species, nM")
plt.title("Change in species over 18 months lecanemab treatment")
plt.show()
fig1.savefig('detailed_model/per_second/species_adcp_00000036.png')

fig2 = plt.figure(2)
plt.plot(time, antibody, label = 'mAb', color='black')
no_months = 18
xticks = [i*(28*24*360) for i in range(0, no_months+1, 3)]
xtick_labels = [i for i in range(0, no_months+1, 3)]
plt.xticks(xticks, xtick_labels)
#plt.legend()
plt.xlabel("Time, months")
plt.ylabel("Species, nM")
plt.title("Change in free lecanemab concentration in brain")
plt.show()
fig2.savefig('detailed_model/per_second/brain_ab_adcp_00000036.png')

fig3 = plt.figure(3)
plt.plot(time, lec_plasma, label = 'mAb', color='black')
no_months = 18
xticks = [i*(28*24*360) for i in range(0, no_months+1, 3)]
xtick_labels = [i for i in range(0, no_months+1, 3)]
plt.xticks(xticks, xtick_labels)
#plt.legend()
plt.xlabel("Time, months")
plt.ylabel("Species, nM")
plt.title("Change in lecanemab concentration in plasma")
plt.show()
fig3.savefig('detailed_model/per_second/plasma_ab_adcp_00000036.png')

fig4 = plt.figure(4)
plt.plot(time, pl_ab, label = 'Plaque-Ab')
plt.plot(time, olig_ab, label = 'Oligomer-Ab')
plt.plot(time, mon_ab, label = 'Monomer-Ab')
no_months = 18
xticks = [i*(28*24*360) for i in range(0, no_months+1, 3)]
xtick_labels = [i for i in range(0, no_months+1, 3)]
plt.xticks(xticks, xtick_labels)
plt.legend()
plt.xlabel("Time, months")
plt.ylabel("Complex, nM")
plt.title("Change in antibody-species complex concentrations")
plt.show()
fig4.savefig('detailed_model/per_second/species_ab_adcp_00000036.png')

fig5 = plt.figure(5)
plt.plot(time, pl_ab_fcr, label = 'Plaque-mAb-FcR')
plt.plot(time, olig_ab_fcr, label = 'Oligomer-mAb-FcR')
no_months = 18
xticks = [i*(28*24*360) for i in range(0, no_months+1, 3)]
xtick_labels = [i for i in range(0, no_months+1, 3)]
plt.xticks(xticks, xtick_labels)
plt.legend()
plt.xlabel("Time, months")
plt.ylabel("Comples, nM")
plt.title("Change in receptor-bound complex concentration")
plt.show()
fig5.savefig('detailed_model/per_second/species_ab_FcR_adcp_00000036.png')

fig6 = plt.figure(6)
plt.plot(time, receptor, label = 'FcR', color='r')
no_months = 18
xticks = [i*(28*24*360) for i in range(0, no_months+1, 3)]
xtick_labels = [i for i in range(0, no_months+1, 3)]
plt.xticks(xticks, xtick_labels)
plt.xlabel("Time, months")
plt.ylabel("Species, nM")
plt.title("Change in FcR concentration")
plt.legend()
plt.show()
fig6.savefig('detailed_model/per_second/receptor_adcp_00000036.png')

fig7 = plt.figure(7)
plt.plot(time, total_ab, label = 'total brain Ab', color='k')
no_months = 18
xticks = [i*(28*24*360) for i in range(0, no_months+1, 3)]
xtick_labels = [i for i in range(0, no_months+1, 3)]
plt.xticks(xticks, xtick_labels)
#plt.legend()
plt.title('Change in total lecanemab concentration in brain (bound and unbound)')
plt.xlabel("Time, months")
plt.ylabel("Species, nM")
plt.show()
fig7.savefig('detailed_model/per_second/brainAb_adcp_00000036.png')

fig8 = plt.figure(8)
plt.plot(time, per_ab, label = 'brain Ab as percentage of plasma', color='k')
no_months = 18
xticks = [i*(28*24*360) for i in range(0, no_months+1, 3)]
xtick_labels = [i for i in range(0, no_months+1, 3)]
plt.xticks(xticks, xtick_labels)
#plt.legend()
plt.title('Brain lecanemab concentration as a percentage of plasma concentration')
plt.xlabel("Time, months")
plt.ylabel("Brain lecanemab/plasma lecanemab, %")
plt.show()
fig8.savefig('detailed_model/per_second/percent_brainAb_adcp_00000036.png')

fig9 = plt.figure(9)
plt.plot(time, plaque_total, label = 'Plaque')
plt.plot(time, olig_total, label = 'Oligomer')
plt.plot(time, mon_total, label = 'Monomer')
no_months = 18
xticks = [i*(28*24*360) for i in range(0, no_months+1, 3)]
xtick_labels = [i for i in range(0, no_months+1, 3)]
plt.xticks(xticks, xtick_labels)
plt.title('Total plaque (bound and unbound)')
plt.xlabel("Time, months")
plt.ylabel("Plaque, nM")
plt.legend()
plt.title('Change in total concentration of plaque (bound and unbound)')
plt.show()
fig9.savefig('detailed_model/per_second/totalPlaque_adcp_00000036.png')

fig10 = plt.figure(10)
plt.plot(time, consumed, label = 'Total plaque phagocytosed')
no_months = 18
xticks = [i*(28*24*360) for i in range(0, no_months+1, 3)]
xtick_labels = [i for i in range(0, no_months+1, 3)]
plt.xticks(xticks, xtick_labels)
#plt.legend()
plt.xlabel("Time, months")
plt.ylabel("Phagocytosed plaque, nM")
plt.title('Cumulative amount of plaque consumed by phagocytosis')
plt.show()
fig10.savefig('detailed_model/per_second/phagocytosedPlaque_adcp_00000036.png')

fig11 = plt.figure(11)
plt.plot(time, new_agg, label = 'Total plaque aggregating')
plt.show()
# fig11.savefig('detailed_model/aggregating_plaque_adcp_000036_old.png')