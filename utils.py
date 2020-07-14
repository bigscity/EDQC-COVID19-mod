#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
import pickle
import json
import csv
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import datetime
import numpy as np
import pandas as pd
import datetime as dt
from pymc import *
from time import time


# In[2]:


def str_to_dt(datestr):
    return dt.datetime.strptime(datestr, '%Y/%m/%d').date()
def dt_to_str(datedt):
    return datedt.strftime('%Y/%m/%d')


# In[3]:


def save_obj(obj, name ):
    with open('objs/'+ name + '.pkl', 'wb') as f:
        pickle.dump(obj, f)

def load_obj(name ):
    with open('objs/' + name + '.pkl', 'rb') as f:
        return pickle.load(f)


# In[4]:


def get_3dayavg(x):
    res = [x[0],]
    for t in range(1, len(x)-1):
        res.append(np.mean([x[t-1], x[t], x[t+1]]))
    res.append(x[-1])
    return np.array(res)


# In[5]:


# Data of cumulative E
Edata = pd.read_csv('modeldata/Edata.csv')


# In[6]:


# Data of cumulative I
Idata = pd.read_csv('modeldata/Idata.csv')


# In[7]:


# Data of cumulative R
Rdata = pd.read_csv('modeldata/Rdata.csv')


# In[8]:


# Average duration from symptom onset to report
confirmrate_data = pd.read_csv('modeldata/confirm_duration.csv')


# In[9]:


# Intensity of human mobility
bmi_data = pd.read_csv('modeldata/mobility.csv')


# In[10]:


# Intensity of human mobility without intervention
bmi_data_nocontrol = pd.read_csv('modeldata/mobility_noctrl.csv')


# In[11]:


# eta1: close-contact tracing, eta2: community-based NAT
eta_data = pd.read_csv('modeldata/eta12.csv')


# In[12]:


# eta220: community-based NAT stopped after June 20
eta220_data = pd.read_csv('modeldata/eta220.csv')


# In[ ]:





# In[13]:


# starting from June 4 (first case onset)
start_t = 2
population = np.array([22000000,])
Etrue = Edata.drop(['ID', 'city'], axis= 1)
Etrue = np.array(Etrue)[:, start_t:]
Itrue = Idata.drop(['ID', 'city',], axis= 1)
Itrue = np.array(Itrue)[:, start_t:]
Rtrue = Rdata.drop(['ID', 'city',], axis= 1)
Rtrue = np.array(Rtrue)[:, start_t:]

Etrue = Etrue - Itrue 
Itrue = Itrue - Rtrue
confirmrate = np.array(confirmrate_data['confirm_duration']).flatten()[start_t:]
confirmrate = get_3dayavg(confirmrate)
bmi = np.array(bmi_data['bmi']).flatten()[start_t:]
bmi_noctrl = np.array(bmi_data_nocontrol['bmi']).flatten()[start_t:]


# In[14]:


alldates = list(Edata.columns)[2+start_t:]
alldates_dt = [dt.datetime.strptime(d, '%Y/%m/%d').date() for d in alldates]
allday_t = len(alldates)


# In[15]:


eta1 = (np.array(eta_data['eta1']))[start_t:]
eta2 = (np.array(eta_data['eta2']))[start_t:]
eta220 = (np.array(eta220_data['eta220']))[start_t:]


# In[16]:


eta1_noctrl = np.zeros(eta1.shape)
eta2_noctrl = np.zeros(eta2.shape)


# In[17]:


'''
N: total population of each region, order by valid_names
initI: first day's infection num, order by valid_names
'''
N = population
initE = [x for x in Etrue[:, 0]]
initI = [x for x in Itrue[:, 0]]
initR = [x for x in Rtrue[:, 0]]


# In[ ]:





# In[18]:


def get_sigma(t, sigma):
    if t < len(sigma):
        return sigma[t]
    return sigma[-1]


# In[19]:


def get_bmi(t, bmi):
    if t < len(bmi):
        return bmi[t]
    while t >= len(bmi):
        t -= 7
    return bmi[t]


# In[20]:


def get_arr(t, arr):
    if t < len(arr):
        return arr[t]
    return arr[-1]


# In[21]:


def sim_seiisr(INPUT, alpha1, bmi, beta1, sigma, t_range):
    T = len(t_range)
    S = np.zeros(T)
    E = np.zeros(T)
    I = np.zeros(T)
    IS = np.zeros(T)
    R = np.zeros(T)
    S[0] = INPUT[0]
    E[0] = INPUT[1]
    I[0] = INPUT[2]
    IS[0] = INPUT[3]
    R[0] = INPUT[4]
    print(INPUT)
    beta1 = 1 / beta1
    for j in range(1, T):
        eta1j = get_arr(j-1, eta1)
        eta2j = get_arr(j-1, eta2)
        bmij = get_bmi(j-1, bmi)
        sigmaj = get_sigma(j-1, sigma)
        S[j] = S[j - 1] - alpha1 * S[j - 1] * bmij * I[j - 1] 
        E[j] = E[j - 1] + alpha1 * S[j - 1] * bmij * I[j - 1]   - beta1 * E[j - 1]
        I[j] = I[j - 1] + (1-eta1j-eta2j) * beta1 * E[j - 1] - sigmaj * I[j - 1]
        IS[j] = IS[j - 1] + (eta1j+eta2j) * beta1 * E[j - 1] - sigmaj * IS[j - 1]
        R[j] = R[j - 1] + sigmaj * (I[j - 1] + IS[j - 1])

    return np.array([S, E, I, IS, R]).T


# In[22]:


def simulate_single(alpha1, e2, beta1, sigma):

    if alpha1 < 0 or beta1 < 0:
        print('Warning ', alpha1, beta1)
    
    INPUT1 = np.zeros((5))
    INPUT1[0] = population[0]
    INPUT1[1] = initE[streetid]
    INPUT1[2] = initI[streetid] * (1 - eta1[0] - eta2[0])
    INPUT1[3] = initI[streetid] * (eta1[0] + eta2[0])
    INPUT1[4] = initR[streetid]
    
    t_range = np.arange(0.0, t, 1.0)
    RES1 = sim_seiisr(INPUT1, alpha1, e2, beta1, sigma, t_range)
    
    return RES1


# In[ ]:





# In[167]:


import numpy as np

def SEIR_MODEL1(Ecumm, Icumm, Rcumm, init_S, init_E, init_I, init_R, bmi, sigma):
    T = len(Ecumm)
    print(T)
    case_data = np.concatenate([Ecumm, Icumm, Rcumm])
    
    alpha1 = Uniform('alpha1', 1e-8, 5e-4, value = 7.9e-8)
    beta1 = Weibull('beta1', alpha = 1.681228, beta = 6.687700 )
    
    @deterministic
    def sim(alpha1 = alpha1, beta1 = beta1):
        S = np.zeros(T)
        E = np.zeros(T)
        I = np.zeros(T)
        IS = np.zeros(T)
        S[0] = init_S
        E[0] = init_E
        I[0] = init_I * (1 - eta1[0] - eta2[0])
        IS[0] = init_I * (eta1[0] + eta2[0])
        cumulative_cases = np.zeros(3*T)
        cumulative_cases[0] = Ecumm[0]
        cumulative_cases[T] = Icumm[0]
        cumulative_cases[2*T] = Rcumm[0]

        for j in range(1, T):
            eta1j = get_arr(j-1, eta1)
            eta2j = get_arr(j-1, eta2)
            bmij = get_bmi(j-1, bmi)
            sigmaj = get_sigma(j-1, sigma)
            
            S[j] = S[j - 1] - alpha1 * S[j - 1] * bmij * I[j - 1] 
            E[j] = E[j - 1] + alpha1 * S[j - 1] * bmij * I[j - 1] - E[j - 1] / beta1
            I[j] = I[j - 1] + (1-eta1j-eta2j) * E[j - 1] / beta1 - sigmaj * I[j - 1]
            IS[j] = IS[j - 1] + (eta1j+eta2j) * E[j - 1] / beta1 - sigmaj * IS[j - 1]
            cumulative_cases[j] = E[j] 
            cumulative_cases[j + T] = cumulative_cases[j + T - 1] + E[j - 1] / beta1
            cumulative_cases[j + 2*T] = cumulative_cases[j + 2*T - 1] + sigmaj * (I[j - 1] + IS[j - 1])
            
        return cumulative_cases[:]
    cases = Lambda('cases', lambda sim = sim : sim)
    A = Poisson('A', mu = cases, value = case_data, observed = True)
    return locals()


# In[186]:


# sample parameters using MCMC methods and save sampled result to Beijing_tracetest.pkl.
# The pymc package seems to achieve different simulation results in each run, though the seed is set.
# Beijing_trace.pkl saves the result shown in the paper.
def sample_from_MCMC():
    streetname = 'Beijing'
    streetid = 0
    e2 = bmi
    sigma = np.divide(1, confirmrate)
    Ecumm = Etrue[streetid, :]
    Icumm = Itrue[streetid, :]
    Rcumm = Rtrue[streetid, :]

    print(N[streetid], initE[streetid], initI[streetid], initR[streetid], streetid, streetname)
    np.random.seed(202018)
    mod1 = SEIR_MODEL1(Ecumm, Icumm + Rcumm, Rcumm, N[streetid], initE[streetid], initI[streetid], initR[streetid], 
                      e2, sigma)
    mc = MCMC(mod1)
    mc.use_step_method(AdaptiveMetropolis, [mod1['alpha1'], mod1['beta1']])#, mod1['s0'], mod1['e0']
    mc.sample(iter = 45000, burn = 1000, thin = 44, verbose = 0) #23000 35000 number of saved parameters simulated = (item - burn) // thin.

    plt.figure(figsize=(8, 6))
    plt.title('SEIR Model for {}'.format(streetname),)
    plt.plot(mod1['case_data'], 's', mec = 'black', color = 'black')
    plt.plot(mod1['cases'].stats()['mean'], color = 'blue', linewidth = 2)
    plt.plot(mod1['cases'].stats()['95% HPD interval'][0], 
             color = 'blue', linewidth = 1, linestyle = 'dotted')
    plt.plot(mod1['cases'].stats()['95% HPD interval'][1],
             color = 'blue', linewidth = 1, linestyle = 'dotted')
    plt.ylabel('Cumulative Number of Cases')
    plt.xlabel('Time (days)')
    plt.grid(True)
    plt.show()
    
    street_params = {}
    street_params[streetname] = {'alpha_trace': mod1['alpha1'].trace(), 'beta_trace': mod1['beta1'].trace()}
    save_obj(street_params, 'Beijing_tracetest')
    

