# analytic solutions to test against

import math
import matplotlib.pyplot as plt
import pandas as pd

CL = 0.0181
Q = 0.0349
V1 = 3.22
V2 = 2.19

alpha1 = Q / V1
alpha2 = Q / V2
beta = CL / V1

co_a = 1
co_b = (alpha1 + alpha2 + beta)
co_c = alpha2*beta

f_constant = 10
f_1 = 10
epsilon = 0.5

root = math.sqrt((co_b)**2 - 4*co_a*co_c)

lambda1 = (-co_b + root)/2
lambda2 = (-co_b - root)/2

# these are both negative hence the baseline scenario is exponential decay

A_constant = ((1+(lambda2/beta))/(lambda1-lambda2))*f_constant
B_constant = -((1+(lambda1/beta))/(lambda1-lambda2))*f_constant

C_sin = (f_1*((alpha2*beta)-1)-alpha2*f_1*(alpha1+alpha2+beta))/((alpha2*beta-1)**2 + (alpha1+alpha2+beta)**2)
D_sin = (f_1*(alpha1+alpha2+beta)+alpha2*f_1*(alpha2*beta - 1))/((alpha2*beta-1)**2 + (alpha1+alpha2+beta)**2)
E_sin = (f_constant/beta)
A_sin = (-D_sin + (lambda2*C_sin)+(lambda2*E_sin)+f_constant)/(lambda1-lambda2)
B_sin = (-(lambda1*C_sin)-(lambda1*E_sin)+D_sin-f_constant)/(lambda1-lambda2)

C_sin_n = (1/epsilon)*((f_1*((alpha2*beta)-(1/epsilon**2))-(alpha1+alpha2+beta)*alpha2*f_1)/((alpha2*beta-(1/epsilon**2))**2 + (((alpha1+alpha2+beta)**2)/epsilon**2)))
D_sin_n = (((alpha1+alpha2+beta)/epsilon)*C_sin_n + alpha2*f_1)/((alpha2*beta-(1/epsilon**2)))
E_sin_n = (f_constant/beta)
A_sin_n = (-(D_sin_n/epsilon) + (lambda2*C_sin_n)+(lambda2*E_sin_n)+f_constant)/(lambda1-lambda2)
B_sin_n = (-(lambda1*C_sin_n)-(lambda1*E_sin_n)+(D_sin_n/epsilon)-f_constant)/(lambda1-lambda2)

y_constant = []
t_constant = []  
dosage_constant = f_constant

for i in range(100):
    L1 = A_constant*math.exp(lambda1*i) + B_constant*math.exp(lambda2*i) + (f_constant/beta)
    t_constant.append(i)
    y_constant.append(L1)

y_sin = []
t_sin = []

def dosage_sin(t):
    return f_constant + f_1*math.sin(t)

for i in range(100):
    L1 = A_sin*math.exp(lambda1*i) + B_sin*math.exp(lambda2*i) + C_sin*math.cos(i) + D_sin*math.sin(i) + (f_constant/beta)
    t_sin.append(i)
    y_sin.append(L1)

y_sin_n = []
t_sin_n = []

def dosage_sin_n(t):
    return f_constant + (f_1)*math.sin(t/epsilon)

#for i in range(100):
i = 0
while i < 100:
    L1 = A_sin_n*math.exp(lambda1*i) + B_sin_n*math.exp(lambda2*i) + C_sin_n*math.cos(i/epsilon) + D_sin_n*math.sin(i/epsilon) + (f_constant/beta)
    t_sin_n.append(i)
    y_sin_n.append(L1)
    i += 0.1

y_sin_nest = []
t_sin_nest = []

epsilon_2 = 0.1

C_sin_nest = (1/epsilon_2)*((f_1*((alpha2*beta)-(1/epsilon_2**2))-(alpha1+alpha2+beta)*alpha2*f_1)/((alpha2*beta-(1/epsilon_2**2))**2 + (((alpha1+alpha2+beta)**2)/epsilon_2**2)))
D_sin_nest = (((alpha1+alpha2+beta)/epsilon_2)*C_sin_nest + alpha2*f_1)/((alpha2*beta-(1/epsilon_2**2)))
E_sin_nest = (f_constant/beta)
A_sin_nest = (-(D_sin_nest/epsilon_2) + (lambda2*C_sin_nest)+(lambda2*E_sin_nest)+f_constant)/(lambda1-lambda2)
B_sin_nest = (-(lambda1*C_sin_nest)-(lambda1*E_sin_nest)+(D_sin_nest/epsilon_2)-f_constant)/(lambda1-lambda2)

for i in range(100):
    L1 = A_sin_nest*math.exp(lambda1*i) + B_sin_nest*math.exp(lambda2*i) + C_sin_nest*math.cos(i/epsilon_2) + D_sin_nest*math.sin(i/epsilon_2) + (f_constant/beta)
    t_sin_nest.append(i)
    y_sin_nest.append(L1)

plt.figure(1)
plt.plot(t_constant, y_constant)
plt.xlabel('time, hours')
plt.ylabel('Lecanemab, mg')
plt.title("Analytical solution with constant dose")
plt.savefig('analytical/constant.png')

plt.figure(2)
plt.plot(t_sin, y_sin)
plt.xlabel('time, hours')
plt.ylabel('Lecanemab, mg')
plt.title("Analytical solution with sinusoidal dose")
plt.savefig('analytical/sine.png')

plt.figure(3)
plt.plot(t_sin_n, y_sin_n)
plt.xlabel('time, hours')
plt.ylabel('Lecanemab, mg')
plt.title("Analytical solution with narrow sinusoidal dose")
plt.savefig('analytical/narrow_sine.png')

plt.figure(4)
plt.plot(t_sin_nest, y_sin_nest)
plt.xlabel('time, hours')
plt.ylabel('Lecanemab, mg')
plt.title("Analytical solution with narrower sinusoidal dose")
plt.savefig('analytical/narrower_sine.png')

# dosage_y1 = []
# dosage_y2 = []
# dosage_y3 = []
# times = []
# i=0
# while i <= 10:
#     dosage_y1.append(dosage_constant)
#     dosage_y2.append(dosage_sin(i))
#     dosage_y3.append(dosage_sin_n(i))
#     times.append(i)
#     i += 0.1

# plt.figure(4)
# plt.plot(times, dosage_y1, 'g', label = 'constant')
# plt.plot(times, dosage_y2, 'b', label = 'sine')
# plt.plot(times, dosage_y3, 'r', label = 'small sine')
# plt.legend()
# plt.title('Test Dosage Regimen')
# plt.savefig('analytical/dosage')

# results = pd.DataFrame()
# results['time'] = t_constant
# results['constant'] = y_constant
# results['sine'] = y_sin
# results['narrow_sine'] = y_sin_n
# results.to_csv('analytical/100_h_solutions.csv')
