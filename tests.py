#
# tests against analytical solutions
# conditions: per hour, for 100 hours
# plot the concentration of lecanemab in the central compartment
#

import random
import matplotlib.pyplot as plt
import pandas as pd

from lecanemab_model import Pharmacokinetic
from solution import Solution

random.seed(1)

constant = Pharmacokinetic(10, 10, 0.5, 'constant')

solver = Solution(constant, 0, 100, 1) # 12960 hours = 18 months; 1080 half days

solutions_c, time = solver.solve()

central_L_c = solutions_c[0]

# solutions_c = solver.solve()

# central_L_c = solutions_c.y[0]

sine = Pharmacokinetic(10, 10, 0.5, 'sine')

solver = Solution(sine, 0, 100, 1) # 12960 hours = 18 months; 1080 half days

solutions_s, time = solver.solve()

central_L_s = solutions_s[0]

# solutions_s = solver.solve()

# central_L_s = solutions_s.y[0]

narrow_sine = Pharmacokinetic(10, 10, 0.5, 'narrow_sine')

solver = Solution(narrow_sine, 0, 100, 1) # 12960 hours = 18 months; 1080 half days

solutions_ns, time = solver.solve()

central_L_ns = solutions_ns[0]

# solutions_ns = solver.solve()

# central_L_ns = solutions_ns.y[0]

# time = solutions_c.t

plt.figure(1)
plt.plot(time, central_L_c)
#plt.axhline(y=half_max_lec, linestyle='dashed', color = 'black')
plt.xlabel("Time, hours")
plt.ylabel("Lecanemab level, mg")
plt.title("Lecanemab profile")
#plt.show()
plt.savefig('analytical/dosefn/constant_euler.png')

plt.figure(2)
plt.plot(time, central_L_s)
#plt.axhline(y=half_max_lec, linestyle='dashed', color = 'black')
plt.xlabel("Time, hours")
plt.ylabel("Lecanemab level, mg")
plt.title("Lecanemab profile")
#plt.show()
plt.savefig('analytical/dosefn/sine_euler.png')

plt.figure(3)
plt.plot(time, central_L_ns)
#plt.axhline(y=half_max_lec, linestyle='dashed', color = 'black')
plt.xlabel("Time, hours")
plt.ylabel("Lecanemab level, mg")
plt.title("Lecanemab profile")
#plt.show()
plt.savefig('analytical/dosefn/narrower_sine_euler.png')

results = pd.DataFrame()
results['time'] = time
results['constant'] = central_L_c
results['sine'] = central_L_s
results['narrow_sine'] = central_L_ns
results.to_csv('analytical/dosefn/euler_solutions.csv')