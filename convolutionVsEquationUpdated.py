# -*- coding: utf-8 -*-
"""
Created on Fri Dec 31 10:47:14 2021

@author: leand
"""



import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


#%%
# # Set up formatting for the movie files
# import matplotlib.animation as animation
# Writer = animation.writers['ffmpeg']
# writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
#%%


def Gaussian(x, sigma, mu):
    return 1/(np.sqrt(2*np.pi) * sigma) * np.exp(-0.5 * ((x - mu)/sigma)^2)


def single_term(x, k, N_a, mu, p_skip, p_slide):
    return 1/(2 * np.pi * N_a * np.sqrt(p_skip * p_slide)) \
            * np.exp(-0.5 * (\
                             x**2/(p_skip * N_a * mu**2)
                             + (x - k * mu)**2 / (p_slide * N_a)
                             ))

def P(N_a, x, p_skip, p_slide, L, mu):
    k_max = int(2 * L / mu)
    k = np.arange(-k_max, k_max + 1)
    X, K = np.meshgrid(x, k)
    Y = single_term(X, K, N_a, mu, p_skip, p_slide)
    return np.sum(Y, 0)


def calcEquation(L, mu, p_slide, simrange):
    x = np.arange(-L, L+1)
    
    lattice = np.zeros(2*L+1)
    lattice[L] = 1

    p_skip = 1-p_slide
    
    result_equation = [lattice.copy()]
    
    for i in range(1, simrange):
        result_equation.append(P(i, x, p_skip, p_slide, L, mu))
        print(i)
        print(sum(P(i, x, p_skip, p_slide, L, mu)))

    result_equation = np.array(result_equation)
    return result_equation


def single_termUpdated(x, k, N_a, mu, p_skip, p_slide):
    return 1/(2 * np.pi * N_a * np.sqrt(p_skip * p_slide)) \
            * np.exp(-0.5 * (\
                             k**2/(p_skip * N_a)
                             + (x - k * mu)**2 / (p_slide * N_a)
                             ))

def PUpdated(N_a, x, p_skip, p_slide, L, mu):
    k_max = int(2 * L / mu)
    k = np.arange(-k_max, k_max + 1)
    X, K = np.meshgrid(x, k)
    Y = single_termUpdated(X, K, N_a, mu, p_skip, p_slide)
    return np.sum(Y, 0)


def calcEquationUpdated(L, mu, p_slide, simrange):
    x = np.arange(-L, L+1)
    
    lattice = np.zeros(2*L+1)
    lattice[L] = 1

    p_skip = 1-p_slide
    
    result_equation = [lattice.copy()]
    
    for i in range(1, simrange):
        result_equation.append(PUpdated(i, x, p_skip, p_slide, L, mu))
        print(i)
        print(sum(PUpdated(i, x, p_skip, p_slide, L, mu)))

    result_equation = np.array(result_equation)
    return result_equation


def convolutionFancySplit(L, mu, p_slide, simrange):

    lattice = np.zeros(2*L+1)
    lattice[L] = 1
    
    p_skip = 1-p_slide
    
    
    PMF_mod = np.zeros(2*mu+3)
    PMF_mod[0] = 0.5 * 0.5 * p_skip * p_slide
    PMF_mod[1] = 0.5 * p_skip * (1-p_slide)
    PMF_mod[2] = 0.5 * 0.5 * p_skip * p_slide
    
    PMF_mod[mu] = 0.5 * p_slide * (1-p_skip)
    PMF_mod[mu+1] = (1-p_slide) * (1-p_skip)
    PMF_mod[mu+2] = 0.5 * p_slide * (1-p_skip)
    
    PMF_mod[-1] = 0.5 * 0.5 * p_skip * p_slide
    PMF_mod[-2] = 0.5 * p_skip * (1-p_slide)
    PMF_mod[-3] = 0.5 * 0.5 * p_skip * p_slide
    
    
    result_mod = [lattice.copy()]
    
    old_mod = lattice.copy()
    
    for i in range(1, simrange):
        new_mod = np.convolve(old_mod, PMF_mod)
        new_mod = new_mod[mu+1:-(mu+1)]
        result_mod.append(1*new_mod)
        old_mod = new_mod.copy()

    
    result_mod = np.array(result_mod)
    return result_mod

def convolutionRegularSplit(L, mu, p_slide, simrange):
        

    lattice = np.zeros(2*L+1)
    lattice[L] = 1
    
    p_skip = 1-p_slide
    
    
    PMF_mod = np.zeros(2*mu+3)
    PMF_mod[1] = 0.5 * p_skip
    
    PMF_mod[mu] = 0.5 * p_slide
    PMF_mod[mu+2] = 0.5 * p_slide
    
    PMF_mod[-2] = 0.5 * p_skip
    
    
    result_mod = [lattice.copy()]
    
    old_mod = lattice.copy()

    
    for i in range(1, simrange):
        new_mod = np.convolve(old_mod, PMF_mod)
        new_mod = new_mod[mu+1:-(mu+1)]
        result_mod.append(1*new_mod)
        old_mod = new_mod.copy()

    
    result_mod = np.array(result_mod)
    return result_mod




#%% PARAMETERS

L = 3000
mu = 50
p_slide = 0.99
simrange = 1000
x = np.arange(-L, L+1)

#%% SIMULATION

result_mod = convolutionFancySplit(L, mu, p_slide, simrange)
result_basic = convolutionRegularSplit(L, mu, p_slide, simrange)
# result_equation = calcEquation(L, mu, p_slide, simrange)
result_equationUpdated = calcEquationUpdated(L, mu, p_slide, simrange)



#%% PRELIMINARY PLOTTING

start = 8
stop = simrange
# stop = int(simrange/2)

N_a = np.arange(start, stop)

result_mod_cumsum = np.cumsum(result_mod[start:stop], 0)
result_basic_cumsum = np.cumsum(result_basic[start:stop], 0)
# result_equation_cumsum= np.cumsum(result_equation[start:stop], 0)
result_equationUpdated_cumsum = np.cumsum(result_equationUpdated[start:stop], 0)

fig_cumsum, ax_cumsum = plt.subplots()

ax_cumsum.set_title('L='+str(L)+', mu='+str(mu)+ ', p_slide='+str(p_slide))

adjustment = -3

ax_cumsum.plot(N_a, result_mod_cumsum[:, L], label = 'Discrete Result')
# ax_cumsum.plot(N_a, result_basic_cumsum[:, L], label = 'Discrete Result')
# ax_cumsum.plot(N_a, result_equation_cumsum[:, L], label = 'Continuous Approximation')
ax_cumsum.plot(N_a, result_equationUpdated_cumsum[:, L] , label = 'Continuous Approximation')
# ax_cumsum.plot(N_a, result_equationUpdated_cumsum[:, L] + adjustment, label = 'Continuous Approximation + ' + str(adjustment))

ax_cumsum.set_ylim(0, ax_cumsum.get_ylim()[1])

ax_cumsum.set_xlabel('N_a')
ax_cumsum.set_ylabel('Sum from ' + str(start)+' to N_a over p(N_a, 0)')

ax_cumsum.legend()
    


#%% EFFECT OF CHANGING INTEGRATION LOWER BOUND

start = 0
stop = 100

# start_2 = 2
# stop_2 = simrange
# stop = int(simrange/2)

N_a = np.arange(start, stop)
# N_a_2 = np.arange(start_2, stop_2)

# result_diff = (result_mod[start:stop] - result_basic[start:stop])**2
result_diff = (result_mod[start:stop] - result_basic[start:stop])
# result_diff_2 = (result_mod[start_2:stop_2] - result_basic[start_2:stop_2])**2


fig_cumsum, ax_cumsum = plt.subplots()

ax_cumsum.set_title('L='+str(L)+', mu='+str(mu)+ ', p_slide='+str(p_slide))

ax_cumsum.plot(N_a[1:], result_diff[1:, L], label='difference')
ax_cumsum.plot(N_a[1:], 0.4/N_a[1:], label='0.4/N_a')

# ax_cumsum.plot(N_a, result_diff[:, L], label = 'Integration from N_a = ' + str(start))
# ax_cumsum.plot(N_a_2, result_diff_2[:, L], label = 'Integration from N_a = ' + str(start_2))
# ax_cumsum.plot(N_a, result_equation_cumsum[:, L], label = 'Continuous Approximation')
# ax_cumsum.plot(result_equation_cumsum[:, L] + adjustment, label = 'Continuous Approximation + ' + str(adjustment))

# ax_cumsum.set_ylim(- ax_cumsum.get_ylim()[1],  ax_cumsum.get_ylim()[1])
# ax_cumsum.set_ylim([0,  2 * result_diff[2]])

ax_cumsum.set_xlabel('N_a')
ax_cumsum.set_ylabel('Difference between strict and independent approach')


ax_cumsum.legend()

#%% ANIMATION CONTROL BLOCK

labels = []
X = []
Y = []

X.append(x)
Y.append(result_mod)
labels.append('Exact Result')

# X.append(x)
# Y.append(result_equation)
# labels.append('Continuous Approximation')

X.append(x)
Y.append(result_equationUpdated)
labels.append('Continuous Approximation')
#%% ANIMATION

fig, ax = plt.subplots()
# ax.set_xlim(-L, L)
ax.set_xlim(-180, 180)
peak = 5

# ax.set_xlim(peak - 30, peak + 30)
ax.set_ylim(0, 1)

LINES = [ax.plot(X[i], Y[i][0])[0] for i in range(len(Y))]
LINES[0].set_linewidth(0.5)
LINES[1].set_linewidth(0.5)

for i, l in enumerate(LINES):
    l.set_label(labels[i])

ax.legend()
ax.set_xlabel('x')
ax.set_ylabel('p(N_a, x)')
ax.set_title('N_a from 0 to '+str(simrange))

counter = 0
def animation_frame(i):
    counter = int(i)
    ax.set_ylim(0, 1.5*max(Y[0][counter]))
    for j, line in enumerate(LINES):
        line.set_ydata(Y[j][counter])
    
animation = FuncAnimation(fig, func=animation_frame, frames=np.arange(0, simrange-1), interval=50)
plt.show()

# animation.to_html5_video()

# SAVE ANIMATION (you need ffmpeg in you path for this to work)

# animation.save('convoVS.mp4', writer = 'ffmpeg', fps = 60, dpi=300)
print("done")

