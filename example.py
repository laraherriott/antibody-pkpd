#
# example flow - currently set up for 18 month continuous bi-weekly treatment
#

import random
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import uniform_filter1d

#from lecanemab_model import LecanemabModel
from dose_switch import LecanemabModel
from solution import Solution
from patient import Patient
from dosage import DoseFn

random.seed(1)

model_type = 'suvr'

model = LecanemabModel(model_type, 10, (14*24), (24*365*3), dose_change=[2, (28*24), (12*1080)], median_patient=True, param_type = 'noiiv')
model()

# patient = Patient()

# dose = DoseFn(10, (14*24), (24*540), patient)

# regimen = dose.get_profile(28*24)

# t = list(range(0, len(regimen), 1))

solver = Solution(model, 0, (24*365*3), 1) # 12960 hours = 18 months; 1080 half days

# solutions, time = solver.solve()

# SUVr = solutions[2]
# central_L = solutions[0]
# peripheral_L = solutions[1]

solutions = solver.solve()

SUVr = solutions.y[2]
central_L = solutions.y[0]
peripheral_L = solutions.y[1]
time = solutions.t

# half_max_lec = max(central_L) / 2

running_mean = uniform_filter1d(central_L[(1080*12):], size=len(central_L[(1080*12):]))



plt.plot(time, central_L, label="Central compartment")
plt.axhline(y=running_mean[-1], color = 'orange', linestyle = 'dashed', label = "Mean level during second 18 months")
#ax.plot(time, peripheral_L, label='peripheral')
#plt.axhline(y=half_max_lec, linestyle='dashed', color = 'black')
plt.xlabel("Time, years")
no_years = 3
xticks = [i*(365*24) for i in range(no_years+1)]
xtick_labels = [i for i in range(no_years+1)]
plt.xticks(xticks, xtick_labels)
plt.ylabel("Lecanemab level, mg")
plt.title("Lecanemab levels following 18 months 10 mg/kg biweekly \n followed by 2 mg/kg monthly")
print(running_mean[-1])
# pos = ax.get_position()
# ax.set_position([pos.x0, pos.y0, pos.width, pos.height*0.9])
plt.legend(loc='lower center', bbox_to_anchor=(0.5, -0.35))
plt.tight_layout()

#plt.show()
plt.savefig('plots/lec_levels_with_switch_2mg.png')

# plt.plot(time, SUVr, color = 'tab:pink')
# plt.xlabel("Time, hours")
# plt.ylabel("SUVr")
# plt.axhline(y=1.17, linestyle='dashed', color = 'black')
# plt.title("Change in SUVr over 18 months bi-weekly treatment")
# plt.ylim((1.0, 1.4))
# plt.show()

# plt.plot(t, regimen, 'x')
# plt.xlabel("Time, hours")
# plt.ylabel("Lecanemab level, mg")
# plt.title("Dosage regimen")
# #plt.savefig('plots/dosage_regimen_2.png')
# plt.show()


