# -*- coding: utf-8 -*-
"""
Created on Fri Aug 13 11:45:33 2021

@author: leand
"""

import numpy as np
# from random import random as rnd

import matplotlib.pyplot as plt


def isAllowed(D, a):    
    return (a >= 0) * (a < D) 


def diffuseUntilEnd(D, params):
    n_walkers_init, mu_skip, l_slide, p_skip = params
    walkers = np.zeros(n_walkers_init)
    
    n_walkers_current = n_walkers_init
    
    visits_to_start = 0
    
    while len(walkers) > 0:
        skip_or_slide = np.where(np.random.rand(n_walkers_current) < p_skip, mu_skip, 1)
        forwards_or_backwards = np.where(np.random.rand(n_walkers_current) < 0.5, 1, -1)
        
        effect = skip_or_slide * forwards_or_backwards
        
        changedWalkers = walkers + effect
        
        
        newWalkers = np.where(isAllowed(D, changedWalkers), changedWalkers, walkers)
        
        walkers = newWalkers
        
        visits_to_start += np.count_nonzero(walkers == 0)
        walkers = walkers[np.where(walkers != (D - 1))]
        
        n_walkers_current = len(walkers)
        
    return visits_to_start / n_walkers_init
            

def calc_n_returns(n_walkers, mu_skip, l_slide, D = np.arange(0, 125)[::1][1:]):
    p_skip = 1 - (l_slide**2 / (1 + l_slide**2))
    
    params = (n_walkers, mu_skip, l_slide, p_skip)
    
    
    D = np.arange(0, 125)[::1][1:]
    n_returns = np.zeros(len(D))
    
    for i, d in enumerate(D):
        print(i, '/', len(D)-1,', ', d)
        n_returns[i] = diffuseUntilEnd(d, params)
        
    return D, n_returns
    

if __name__ == '__main__':
        
    n_walkers = 10000
    
    mu_skip = 50
    l_slide = 5
    
    p_skip = 1 - (l_slide**2 / (1 + l_slide**2))
    
    params = (n_walkers, mu_skip, l_slide, p_skip)
    
    
    D = np.arange(0, 125)[::1][1:]
    n_returns = np.zeros(len(D))
    
    for i, d in enumerate(D):
        print(i, '/', len(D)-1,', ', d)
        n_returns[i] = diffuseUntilEnd(d, params)
        
    plt.plot(D, n_returns)
    
