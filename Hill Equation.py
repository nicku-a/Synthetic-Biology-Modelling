#!/usr/bin/env python
# coding: utf-8

# In[8]:


from scipy.integrate import odeint
import matplotlib.pyplot as plt
import numpy as np


# In[9]:


# System of equations
def hill_model(state, t, params):
    mRNA, prot = state
    k1, k2, d1, d2, km, n, r = params

    mRNA_d = (k1*km**n)/(km**n + r**n) - d1*mRNA
    prot_d = k2*mRNA - d2*prot

    return [mRNA_d, prot_d]


# In[35]:


t  = np.linspace(0, 400., 1000)

# Initial conditions
mRNA0 = 0
prot0 = 0

# Parameters
d1 = 0.693/5
d2 = 0.693/35
k1 = 2.5*d1
k2 = 1000*d2
km = 100
n = [1,2,3]
r = 100

params = np.zeros([3,7])

for i in range(len(n)):
    params[i] = [k1,k2,d1,d2,km,n[i],r]

# Create list of solutions for different rates
solutions = [odeint(hill_model, [mRNA0, prot0], t, args=(p,)) for p in params]


# In[37]:


# Choose individual solution from list
for sn in enumerate(solutions):

    solution = sn[1]
    
    # Concentrations of mRNA and protein
    mRNA_sol = solution[:, 0]
    prot_sol = solution[:, 1]

    # Plot solution
    plt.figure()
    plt.plot(t, mRNA_sol, 'b', label='[mRNA]')
    plt.title('Simulation {}: n={}'.format(sn[0]+1,params[sn[0]][5]))
    plt.xlabel('Time')
    plt.ylabel('Concentration')
    plt.legend(loc=0)
    
    plt.figure()
    plt.plot(t, prot_sol, 'y', label='[protein]')
    plt.title('Simulation {}: n={}'.format(sn[0]+1,params[sn[0]][5]))
    plt.xlabel('Time')
    plt.ylabel('Concentration')
    plt.legend(loc=0)


# In[ ]:




