import numpy as np
means = [1.34, 0.232, 1.54, 75.0, 0.0842, 0.367, 0.00155, 4.06, 0.468, 0.00313]
rse = [0.873, 11.1, 11.8, 19.6, 4.28, 1.97, 9.32, 1.61, 20.7, 15.6]
se = [estimate*(error/100) for estimate,error in zip(means,rse)]
var = [value**2 for value in se]
cov = np.zeros((len(var), len(var)))
np.fill_diagonal(cov, var)

parameter_sample = np.random.multivariate_normal(means, cov, 2)

print(parameter_sample)