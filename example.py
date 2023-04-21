#
# example flow - currently set up for 18 month continuous bi-weekly treatment
#

import random
import matplotlib.pyplot as plt

from lecanemab_model import LecanemabModel
from solution import Solution
from patient import Patient
from dosage import DoseFn

random.seed(1)

model_type = 'suvr'

model = LecanemabModel(model_type, 10, (14*24), (24*540), param_type = 'hour')
model()

patient = Patient()

dose = DoseFn(10, (14*24), (24*540), patient)

regimen = dose.get_profile(28*24)

t = list(range(0, len(regimen), 1))

solver = Solution(model, 0, 12960, 1) # 12960 hours = 18 months; 1080 half days

solutions, time = solver.solve()

SUVr = solutions[2]
central_L = solutions[0]
peripheral_L = solutions[1]

# solutions = solver.solve()

# SUVr = solutions.y[2]
# central_L = solutions.y[0]
# peripheral_L = solutions.y[1]
# time = solutions.t

half_max_lec = max(central_L) / 2


plt.plot(time, central_L, label="central")
plt.plot(time, peripheral_L, label='peripheral')
#plt.axhline(y=half_max_lec, linestyle='dashed', color = 'black')
plt.xlabel("Time, hours")
plt.ylabel("Lecanemab level, mg")
plt.title("Decay in lecanemab levels following treatment")
plt.legend()
plt.show()
#plt.savefig('plots/lec_levels_decay.png')

plt.plot(time, SUVr, color = 'tab:pink')
plt.xlabel("Time, hours")
plt.ylabel("SUVr")
plt.axhline(y=1.17, linestyle='dashed', color = 'black')
plt.title("Change in SUVr over 18 months bi-weekly treatment")
plt.ylim((1.0, 1.4))
plt.show()

plt.plot(t, regimen, 'x')
plt.xlabel("Time, hours")
plt.ylabel("Lecanemab level, mg")
plt.title("Dosage regimen")
#plt.savefig('plots/dosage_regimen_2.png')
plt.show()


