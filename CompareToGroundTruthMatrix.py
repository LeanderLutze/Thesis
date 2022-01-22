# -*- coding: utf-8 -*-
"""
Created on Sat Jan 22 13:41:57 2022

@author: leand
"""

import numpy as np

import matplotlib.pyplot as plt

#%% UNIQUE CODE


def create_matrix_nReturns_largeFlank_bothMove(D_max, params, D_actual):
    mu_skip, l_slide = params
    
    p_skip = 1 - (l_slide**2 / (1 + l_slide**2))
    
    p_slide = 1 - p_skip
    
    nDmu = 0
    
    D = nDmu*D_max + nDmu*mu_skip + D_actual
    
    A = np.zeros((D, D))
    
    
    #MODIFIED
    locationStartingTrap = [int(0.5 * (D - D_actual))]
    
    locationFinalTraps = [int(0.5 * (D - D_actual)) + D_actual - 1]
    
    absorbingTraps = locationStartingTrap + locationFinalTraps
    
    trapLocations = [locationStartingTrap, locationFinalTraps]
    
    for i in range(len(A)):
        if i - 1 >= 0:
            A[i][i - 1] += 0.5 * p_slide
        else:
            A[i][i] += 0.5 * p_slide
        
        if i + 1 < D:
            A[i][i + 1] += 0.5 * p_slide
        else:
            A[i][i] += 0.5 * p_slide
        
        if i - mu_skip >= 0:
            A[i][i - mu_skip] += 0.5 * (1 - p_slide)
        else:
            A[i][i] += 0.5 * (1 - p_slide)
        
        if i + mu_skip < D:
            A[i][i + mu_skip] += 0.5 * (1 - p_slide)
        else:
            A[i][i] += 0.5 * (1 - p_slide)
    
    prev = [A[a][a] for a in absorbingTraps]
    
    for a in absorbingTraps:
        for j in range(D):
            A[j][a] = 0
        
        A[a][a] = 1
    
    return A, trapLocations, prev

#%% COMMON CODE

def computeLimitMatrix(A):
    B = np.array([[a for a in aa] for aa in A])
    for i in range(100):
        B = np.matmul(B, B)
    return B


def calcAverageNumberVisits(p_nextToTrap):
    return p_nextToTrap / (1 - p_nextToTrap)


def extractP_nextToTrap(A, B, s, prev):
    # print(A)
    
    row = A[s]
    row[s] = prev[0]
    return sum(row * B[s])
    

def runAverageNumberVisits(D_max, D_actual, params, flank = 0):    
    # A, trapLocations, prev = create_matrix_nReturns_noFlank_notSymmetric_absorbingTrapMoves(D_max, params)
    # A, trapLocations, prev = create_matrix_nReturns_noFlank_notSymmetric_startingPointMoves(D_max, params, D_actual)
    A, trapLocations, prev = create_matrix_nReturns_largeFlank_bothMove(D_max, params, D_actual)

    
    
    B = computeLimitMatrix(A)
    
    # print(len(A), trapLocations)
    
    locationStartingTrap, locationFinalTraps = trapLocations
    
    
    
    
    s = locationStartingTrap[0]
    
    p_return = extractP_nextToTrap(A, B, s, prev)
    

    # p_return = B[s][s+1]

    
    N_returns = calcAverageNumberVisits(p_return)
    
    return N_returns
    

def executeOverDifferentValuesAndCollect(D_max, DDD_actual, MU_SKIP, L_SLIDE, flank = 0):
    N_RETURNS = [0 for i in range(len(DDD_actual))]
    
    for i in range(len(DDD_actual)):
        mu_skip = MU_SKIP[i]
        l_slide = L_SLIDE[i]
        params = [mu_skip, l_slide]
        
        DD_actual = DDD_actual[i]
        N_returns = np.zeros(len(DD_actual))
        
        for j,D_actual in enumerate(DD_actual):
            print("RUN",i,", ",j,"/",len(DD_actual) - 1, ", ", D_actual, "/", DD_actual[-1])
            N_returns[j] = runAverageNumberVisits(D_max, D_actual, params)
        plt.plot(DD_actual, N_returns)
        plt.show()
        
        N_RETURNS[i] = N_returns
    
    for i, N in enumerate(N_RETURNS):
        plt.plot(DDD_actual[i], N)
    
    plt.show()
    return N_RETURNS

#%%

def calc_n_returns(mu_skip, l_slide, D = np.arange(0, 125)[::1][1:]):
    DD_actual = D
    params = [mu_skip, l_slide]
    DD_max = max(DD_actual + 2)
    N_returns = np.zeros(len(DD_actual))
    for i,D_actual in enumerate(DD_actual):
        print(i, '/', len(DD_actual)-1,', ', D_actual, '/', DD_actual[-1])
        N_returns[i] = runAverageNumberVisits(DD_max, D_actual, params)
        
    return DD_actual, N_returns
    


#%% EXECUTION
# ______________________________________________________________________________

if __name__ == '__main__':
    
    mu_skip = 50
    l_slide = 5
    
    DD_actual = np.arange(0, 125)[::1][1:]
    
    params = [mu_skip, l_slide]
    
    
    # DD_actual = np.arange(2, 125)[::5]
    DD_max = max(DD_actual + 2)
    N_returns = np.zeros(len(DD_actual))
    
    for i,D_actual in enumerate(DD_actual):
        print(i, '/', len(DD_actual)-1,', ', D_actual, '/', DD_actual[-1])
        N_returns[i] = runAverageNumberVisits(DD_max, D_actual, params)
        
    plt.plot(DD_actual, N_returns)
