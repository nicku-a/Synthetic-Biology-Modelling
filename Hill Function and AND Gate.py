#!/usr/bin/env python
# coding: utf-8

# In[286]:


from scipy.integrate import odeint
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import pandas as pd
import os
from matplotlib.colors import LogNorm
from mpl_toolkits.mplot3d import Axes3D


# In[246]:


# System of equations
def hill_model(params, log):
    
    k, n, K1, a = params
    
    f = [k*(a+i**n/((K1*10**(-3))**n+i**n)) for i in log]

    return [f]


# In[276]:


parameters


# In[277]:


def find_R_S(I, A):
    
    k = [2071,12900]
    n = [1.282,1.323]
    K1 = [0.287,0.513]
    a = [0.0037,0.0013]
    
    R = k[0]*(a[0]+I**n[0]/((K1[0]*10**(-3))**n[0]+I**n[0]))
    S = k[1]*(a[0]+A**n[0]/((K1[0]*10**(-3))**n[0]+A**n[0]))
    
    return [R,S]


# In[278]:


# System of equations
def AND_gate(params, r, s):
    
    kr, ks, nr, ns = params
    
    f = (r/kr)**nr*(s/ks)**ns/((1+(r/kr)**nr)*(1+(s/ks)**ns))

    return [f]


# In[267]:


base_dir = 'C:/Users/Nicholas/OneDrive/Biomedical Engineering MEng/Year 3/Group Project'
filename = 'HillFunctionData.csv'

os.chdir(base_dir)
parameters = pd.read_csv(filename, sep = ',', index_col=False, skiprows = 0)

k_new = parameters["k"].str.split("+", n = 1, expand = True)
parameters["k"]= k_new[0]
parameters["k+"]= k_new[1]

n_new = parameters["n"].str.split("+", n = 1, expand = True)
parameters["n"]= n_new[0]
parameters["n+"]= n_new[1]

k1_new = parameters["K1"].str.split("+", n = 1, expand = True)
parameters["K1"]= k1_new[0]
parameters["K1+"]= k1_new[1]

a_new = parameters["alpha"].str.split("+", n = 1, expand = True)
parameters["alpha"]= a_new[0]
parameters["alpha+"]= a_new[1]


# In[355]:


log = np.logspace(-12,-2, 100)

# Parameters
params = np.zeros([len(parameters),4])
for i in range(len(parameters)):
    params[i] = [parameters.at[i,'k'],parameters.at[i,'n'],parameters.at[i,'K1'],parameters.at[i,'alpha']]

# Create list of solutions for different parameters
solutions = [hill_model(p, log) for p in params]

# Parameters
selected_params = [206.1,3135,2.381,1.835]

iptg = arab = np.logspace(-7,-2, 100)
IPTG, ARAB = np.meshgrid(iptg, arab)

R, S = find_R_S(IPTG, ARAB)

# Create list of solutions for different rates
solution = AND_gate(selected_params, R, S)[0]


# In[356]:


np.savetxt('IPTG.csv', IPTG, delimiter=',')
np.savetxt('ARAB.csv', ARAB, delimiter=',')
np.savetxt('AND_gate_solution.csv', solution, delimiter=',')


# In[357]:


fig, ax = plt.subplots(nrows=2, ncols=3,  figsize=(25, 9), squeeze=False)

# Figure d
for i in range(6):
    ax[0,0].plot(log, solutions[i][:][0], label=parameters.at[i,'promoter/RBS'])
ax[0,0].set_xscale('log')
ax[0,0].set_xlabel('[IPTG] (M)')
ax[0,0].set_ylabel('Fluo/OD600 (a.u.)')
ax[0,0].legend(loc=0)
ax[0,0].set_xlim([10**(-7),10**(-2)])

# Figure e
for i in range(6,12): 
    ax[0,1].plot(log, solutions[i][:][0], label=parameters.at[i,'promoter/RBS'])
#     plt.plot(log, solutions[i][:][0], label=parameters.at[i,'promoter/RBS'])
    ax[0,1].set_xscale('log')
    ax[0,1].set_xlabel('[Arabinose] (M)')
    ax[0,1].set_ylabel('Fluo/OD600 (a.u.)')
    ax[0,1].legend(loc=0)
    ax[0,1].set_xlim([10**(-7),10**(-2)])

# Figure f
for i in range(12,18):
    ax[0,2].plot(log, solutions[i][:][0], label=parameters.at[i,'promoter/RBS'])
ax[0,2].set_xscale('log')
ax[0,2].set_xlabel('[AHL] (M)')
ax[0,2].set_ylabel('Fluo/OD600 (a.u.)')
ax[0,2].legend(loc=0)
ax[0,2].set_xlim([10**(-12),10**(-7)])

# Figure g
ax[1,0].plot(log, solutions[5][:][0], label=parameters.at[5,'promoter/RBS'])
ax[1,0].plot(log, solutions[18][:][0], label=parameters.at[18,'promoter/RBS'])
ax[1,0].set_xscale('log')
ax[1,0].set_xlabel('[IPTG] (M)')
ax[1,0].set_ylabel('Fluo/OD600 (a.u.)')
ax[1,0].legend(loc=0)
ax[1,0].set_xlim([10**(-7),10**(-2)])

# Figure h
ax[1,1].plot(log, solutions[9][:][0], label=parameters.at[9,'promoter/RBS'])
ax[1,1].plot(log, solutions[19][:][0], label=parameters.at[19,'promoter/RBS'])
ax[1,1].set_xscale('log')
ax[1,1].set_xlabel('[Arabinose] (M)')
ax[1,1].set_ylabel('Fluo/OD600 (a.u.)')
ax[1,1].legend(loc=0)
ax[1,1].set_xlim([10**(-7),10**(-2)])

# Figure i
ax[1,2].plot(log, solutions[15][:][0], label=parameters.at[15,'promoter/RBS'])
ax[1,2].plot(log, solutions[20][:][0], label=parameters.at[20,'promoter/RBS'])
ax[1,2].set_xscale('log')
ax[1,2].set_xlabel('[AHL] (M)')
ax[1,2].set_ylabel('Fluo/OD600 (a.u.)')
ax[1,2].legend(loc=0)
ax[1,2].set_xlim([10**(-12),10**(-6)])

plt.show()


# In[358]:


fig, ax = plt.subplots()
im = ax.imshow(solution, cmap=cm.plasma, norm=LogNorm(vmin=0.01, vmax=1), extent=[10**(-7), 10**(-2), 10**(-7), 10**(-2)])
ax.set_xlabel('[IPTG] (M)')
ax.set_ylabel('[Arabinose] (M)')
plt.show()

fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(IPTG, ARAB, solution, cmap=cm.plasma)
# ax.set_xscale('log')
# ax.set_yscale('log')
ax.set_xlabel('[IPTG] (M)')
ax.set_ylabel('[Arabinose] (M)')
plt.show()


# In[ ]:




