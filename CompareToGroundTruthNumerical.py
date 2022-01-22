# -*- coding: utf-8 -*-
"""
Created on Fri Aug 13 12:49:22 2021

@author: leand
"""

import numpy as np
import matplotlib.pyplot as plt

def isAllowed(D, a):    
    return (a >= 0) * (a < D) 


def diffuseUntilEnd(D, params):
    mu_skip, l_slide, p_skip = params
    
    p_slide = 1 - p_skip
    
    p_dist = np.zeros(D)
    p_dist[0] = 1
    
    cum_p_start = 0
    
    
    for i in range(100000):
        new_p_dist = np.zeros(D)
        
        '''
        Difference may be coming from the fact that I am forcing particles
        to remain near the edges, which would count them at zero
        '''
        
        new_p_dist[1:] += p_slide * 0.5 * p_dist[:-1]
        new_p_dist[0] += p_slide * 0.5 * p_dist[0]
        
        new_p_dist[:-1] += p_slide * 0.5 * p_dist[1:]
        new_p_dist[-1] += p_slide * 0.5 * p_dist[-1]
        
        new_p_dist[mu_skip:] += p_skip * 0.5 * p_dist[:-mu_skip]
        new_p_dist[:mu_skip] += p_skip * 0.5 * p_dist[:mu_skip]
        
        new_p_dist[:-mu_skip] += p_skip * 0.5 * p_dist[mu_skip:]
        new_p_dist[-mu_skip:] += p_skip * 0.5 * p_dist[-mu_skip:]
        
        new_p_dist[D - 1] = 0
        
        cum_p_start += new_p_dist[0]
        
        # print(sum(p_dist))
        
        p_dist = new_p_dist
        
        
    return cum_p_start


def calc_n_returns(mu_skip, l_slide, D = np.arange(0, 125)[::1][1:]):
    p_skip = 1 - (l_slide**2 / (1 + l_slide**2))
    
    params = (mu_skip, l_slide, p_skip)
    n_returns = np.zeros(len(D))
    
    for i, d in enumerate(D):
        print(i, '/', len(D)-1,', ', d)
        n_returns[i] = diffuseUntilEnd(d, params)
        
    return D, n_returns


if __name__ == '__main__':
    
    mu_skip = 50
    l_slide = 5
    
    p_skip = 1 - (l_slide**2 / (1 + l_slide**2))
    
    params = (mu_skip, l_slide, p_skip)
    
    
    D = np.arange(0, 125)[::1][1:]
    n_returns = np.zeros(len(D))
    
    for i, d in enumerate(D):
        print(i, '/', len(D)-1,', ', d)
        n_returns[i] = diffuseUntilEnd(d, params)
        
    plt.plot(D, n_returns)