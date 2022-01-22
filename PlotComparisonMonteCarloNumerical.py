# -*- coding: utf-8 -*-
"""
Created on Sat Jan 22 13:08:25 2022

@author: leand
"""

#%% IMPORTS

import numpy as np
import matplotlib.pyplot as plt

import GroundTruthRawMonteCarlo as GTRMC
import CompareToGroundTruthNumerical as CTGTN
import CompareToGroundTruthMatrix as CTGTM


#%% GENERAL PARAMETERS
mu_skip = 50
l_slide = 5

D = np.arange(1, 125)[::1][1:]

#%% MONTE CARLO

n_walkers = 100

D_MC_10000, n_returns_MC_10000 = GTRMC.calc_n_returns(n_walkers, mu_skip, l_slide, D)

#%%

n_walkers = 100000
D_MC_100000, n_returns_MC_100000 = GTRMC.calc_n_returns(n_walkers, mu_skip, l_slide, D)

#%% NUMERICAL

D_NUM , n_returns_NUM = CTGTN.calc_n_returns(mu_skip, l_slide, D)


#%% MATRIX

D_MAT, n_returns_MAT = CTGTM.calc_n_returns(mu_skip, l_slide, D)

#%% PLOTTING

plt.title('Comparison: Monte Carlo vs. Numerical vs. Matrix Based')
plt.ylabel('Average Number of Returns to the Initial Trap')
plt.xlabel('Distance Between Traps L')


plt.plot(D_MC_10000, n_returns_MC_10000, color = 'grey', linestyle='solid', linewidth = 0.5, label='Monte Carlo Simulation, N_particles=100')
plt.plot(D_MC_100000, n_returns_MC_100000, color = 'grey', linestyle='solid', linewidth = 3, label='Monte Carlo Simulation, N_particles=100000')
plt.plot(D_NUM, n_returns_NUM, linestyle='solid', label='Numerical Simulation')
plt.plot(D_MAT, n_returns_MAT, linestyle='dashed', label='Matrix Based Approach')

plt.legend()

plt.savefig('ComparisonSimulationMethods', dpi=400)