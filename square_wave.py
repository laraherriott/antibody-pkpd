# plotting square wave

import math
import random
import matplotlib.pyplot as plt
import pandas as pd
from lecanemab_model import Pharmacokinetic
from solution import Solution

A = 10
f = 0.5
delta = 0.1

dose_list = [0, 24, 48, 72]

t = []
y = []

i = -0.5
while i < 1.5:
    sol = 0
    for n in dose_list:
        if i >= n-0.5 and i <= (n+1.5):
            sol =  ((A/2)/math.atan(1/delta))*(math.atan(math.sin(2*math.pi*i*f)/delta)) + A/2
    y.append(sol)
    t.append(i)
    i += 0.1


plt.plot(t, y)
#plt.show()
plt.savefig('plots/square_wave_single.png')

# random.seed(1)

# square = Pharmacokinetic(10, 10, 0.5, 'square')

# solver = Solution(square, 0, 150, 1) # 12960 hours = 18 months; 1080 half days

# solutions_sq, time = solver.solve()

# central_L_sq = solutions_sq[0]

# plt.plot(time, central_L_sq)
# #plt.axhline(y=half_max_lec, linestyle='dashed', color = 'black')
# plt.xlabel("Time, hours")
# plt.ylabel("Lecanemab level, mg")
# plt.title("Lecanemab profile")
# plt.show()
# #plt.savefig('analytical/dosefn/square_euler.png')