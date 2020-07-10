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


print('Loading data...')
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
print('Initialization')
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





# ## Define Simulation Function

# In[23]:


def sim_seiisr_exp(INPUT, alpha1, bmi, beta1, sigma, t_range, eta1, eta2):
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


# In[24]:


def simulate_exp(alpha1, bmi, beta1, sigma, eta1, eta2):

#     k = 30
    
    if alpha1 < 0 or beta1 < 0:
        print('Warning ', alpha1, beta1)
    
    INPUT1 = np.zeros((5))
    INPUT1[0] = population[0]
    INPUT1[1] = initE[streetid]
    INPUT1[2] = initI[streetid] * (1 - eta1[0] - eta2[0])
    INPUT1[3] = initI[streetid] * (eta1[0] + eta2[0])
    INPUT1[4] = initR[streetid]
    
    t_range = np.arange(0.0, t, 1.0)
    RES1 = sim_seiisr_exp(INPUT1, alpha1, bmi, beta1, sigma, t_range, eta1, eta2)
    
    return RES1


# ## Start Simulation

# In[25]:


street_params = load_obj('Beijing_trace')
streetname = 'Beijing'
streetid = 0
e2 = bmi
sigma = np.divide(1, confirmrate)


# In[26]:


params = street_params['Beijing']

t = len(alldates) + 0
based_time = str_to_dt(alldates[0])
t_range_subdt = [based_time + dt.timedelta(days = x) for x in range(t)]


# In[27]:


def get_result_with_CI(params, bmi, eta1, eta2):
    allres = []
    alphas = params['alpha_trace']
    betas = params['beta_trace']
    
    for alpha, beta in zip(alphas, betas):
        allres.append(simulate_exp( 
                          alpha, bmi, 
                          beta,  sigma, eta1, eta2
                            )
                     )
#     print(len(allres))
    result = np.zeros(allres[0].shape)
    resultup = np.zeros(allres[0].shape)
    resultdown = np.zeros(allres[0].shape)
    for i in range(allres[0].shape[0]):
        for j in range(allres[0].shape[1]):
            tmp = [x[i, j] for x in allres]
            result[i, j] = np.percentile(tmp, 50)
            resultup[i, j] = np.percentile(tmp, 75)
            resultdown[i, j] = np.percentile(tmp, 25)
    return result, resultup, resultdown


# In[28]:


def save_resultforsimu(result, resultup, resultdown, orgres, orgresup, orgresdown, filename, dt_range):
    df = pd.DataFrame(data = {
        'Date': dt_range,
        'conbined_intervention_I': orgres[:, 2] + orgres[:, 3] + orgres[:, 4],
        'conbined_intervention_I_up': orgresup[:, 2] + orgresup[:, 3] + orgresup[:, 4],
        'conbined_intervention_I_down': orgresdown[:, 2] + orgresdown[:, 3] + orgresdown[:, 4],
        'reported_I': Itrue[streetid, :] + Rtrue[streetid, :],
        'simulation_I': result[:, 2] + result[:, 3] + result[:, 4],
        'simulation_I_up': resultup[:, 2] + resultup[:, 3] + resultup[:, 4],
        'simulation_I_down': resultdown[:, 2] + resultdown[:, 3] + resultdown[:, 4],
    })
    df.to_csv(filename, index=False)


# ## 1. Under all 3 TIs (actual scenario, Fig.4A)

# In[29]:


print('Simulating scenario under all TIs (actual scenario) and generating graphs...')
orgres, orgresup, orgresdown = get_result_with_CI(params, bmi, eta1, eta2)
save_resultforsimu(orgres, orgresup, orgresdown, orgres, orgresup, orgresdown, './simuresult/resultsimu111.csv', t_range_subdt)


# In[30]:


from matplotlib.dates import DateFormatter
formatter = DateFormatter("%b %d")
plt.figure(figsize=(10, 5))
plt.plot(t_range_subdt, orgres[:, 2] + orgres[:, 3] + orgres[:, 4], 'b', label = 'Estimated')
plt.fill_between(t_range_subdt, orgresdown[:, 2] + orgresdown[:, 3] + orgresdown[:, 4], orgresup[:, 2] + orgresup[:, 3] + orgresup[:, 4], alpha=0.5, color='lightblue')
plt.plot(t_range_subdt, Itrue[streetid, :] + Rtrue[streetid, :], 'ko', label = 'Reported')
plt.xlabel('Date')
plt.ylabel('Number of cases')
plt.gcf().axes[0].xaxis.set_major_formatter(formatter)
plt.gcf().autofmt_xdate()
plt.legend(loc='upper left')
plt.savefig('./simuresult/Fig4A.pdf', bbox_inches='tight')


# ## 2. Under no TIs (Fig.4B)

# In[31]:


print('Simulating scenario under no TIs and generating graphs...')
plt.figure(figsize=(10, 5))
result, resultup, resultdown = get_result_with_CI(params, bmi_noctrl, eta1_noctrl, eta2_noctrl)
save_resultforsimu(result, resultup, resultdown, orgres, orgresup, orgresdown, './simuresult/resultsimu000.csv', t_range_subdt)
plt.plot(t_range_subdt, orgres[:, 2] + orgres[:, 3] + orgres[:, 4], 'b', label = 'Combined interventions')
plt.fill_between(t_range_subdt, orgresdown[:, 2] + orgresdown[:, 3] + orgresdown[:, 4], orgresup[:, 2] + orgresup[:, 3] + orgresup[:, 4], alpha=0.5, color='lightblue')
plt.plot(t_range_subdt, result[:, 2] + result[:, 3] + result[:, 4], 'orchid', label = 'No TIs')
plt.fill_between(t_range_subdt, resultdown[:, 2] + resultdown[:, 3] + resultdown[:, 4], resultup[:, 2] + resultup[:, 3] + resultup[:, 4], alpha=0.7, color='lavender')
plt.xlabel('Date')
plt.ylabel('Number of cases')
plt.gcf().axes[0].xaxis.set_major_formatter(formatter)
plt.gcf().autofmt_xdate()
# plt.grid()
plt.legend(loc='upper left')
# plt.savefig('./NO_intervention.png', dpi=500, bbox_inches='tight')
plt.savefig('./simuresult/Fig4B.pdf', bbox_inches='tight')


# ## 3. Under 1 or 2 TIs (Fig.4C-4H)

# In[32]:


print('Simulating scenarios under 1 or 2 TIs and generating graphs...')
fig, axes = plt.subplots(2, 3, figsize=(20, 12))

# 231 : NO Localized lockdown
ax = plt.subplot(231)
result, resultup, resultdown = get_result_with_CI(params, bmi_noctrl, eta1, eta2)
save_resultforsimu(result, resultup, resultdown, orgres, orgresup, orgresdown, './simuresult/resultsimu011.csv', t_range_subdt)
plt.plot(t_range_subdt, orgres[:, 2] + orgres[:, 3] + orgres[:, 4], 'b', label = 'Combined interventions')
plt.fill_between(t_range_subdt, orgresdown[:, 2] + orgresdown[:, 3] + orgresdown[:, 4], orgresup[:, 2] + orgresup[:, 3] + orgresup[:, 4], alpha=0.5, color='lightblue')
plt.plot(t_range_subdt, result[:, 2] + result[:, 3] + result[:, 4], 'r', label = 'No localized lockdown')
plt.fill_between(t_range_subdt, resultdown[:, 2] + resultdown[:, 3] + resultdown[:, 4], resultup[:, 2] + resultup[:, 3] + resultup[:, 4], alpha=0.5, color='pink')
plt.xlabel('Date')
plt.ylabel('Number of cases')
ax.xaxis.set_major_formatter(formatter)
plt.gcf().autofmt_xdate()
# plt.grid()
plt.legend(loc='upper left')

# 232: NO close-contact tracing
ax = plt.subplot(232)
result, resultup, resultdown = get_result_with_CI(params, bmi, eta1_noctrl, eta2)
save_resultforsimu(result, resultup, resultdown, orgres, orgresup, orgresdown, './simuresult/resultsimu101.csv', t_range_subdt)
plt.plot(t_range_subdt, orgres[:, 2] + orgres[:, 3] + orgres[:, 4], 'b', label = 'Combined interventions')
plt.fill_between(t_range_subdt, orgresdown[:, 2] + orgresdown[:, 3] + orgresdown[:, 4], orgresup[:, 2] + orgresup[:, 3] + orgresup[:, 4], alpha=0.5, color='lightblue')
plt.plot(t_range_subdt, result[:, 2] + result[:, 3] + result[:, 4], 'orange', label = 'No close-contact tracing')
plt.fill_between(t_range_subdt, resultdown[:, 2] + resultdown[:, 3] + resultdown[:, 4], resultup[:, 2] + resultup[:, 3] + resultup[:, 4], alpha=0.5, color='wheat')
plt.xlabel('Date')
plt.ylabel('Number of cases')
ax.xaxis.set_major_formatter(formatter)
plt.gcf().autofmt_xdate()
# plt.grid()
plt.legend(loc='upper left')

# 233: NO NAT (all time or after 20 June)
ax = plt.subplot(233)
result, resultup, resultdown = get_result_with_CI(params, bmi, eta1, eta2_noctrl)
save_resultforsimu(result, resultup, resultdown, orgres, orgresup, orgresdown, './simuresult/resultsimu110.csv', t_range_subdt)
plt.plot(t_range_subdt, orgres[:, 2] + orgres[:, 3] + orgres[:, 4], 'b', label = 'Combined interventions')
plt.fill_between(t_range_subdt, orgresdown[:, 2] + orgresdown[:, 3] + orgresdown[:, 4], orgresup[:, 2] + orgresup[:, 3] + orgresup[:, 4], alpha=0.5, color='lightblue')
plt.plot(t_range_subdt, result[:, 2] + result[:, 3] + result[:, 4], 'orchid', label = 'No community-based NAT all time')
plt.fill_between(t_range_subdt, resultdown[:, 2] + resultdown[:, 3] + resultdown[:, 4], resultup[:, 2] + resultup[:, 3] + resultup[:, 4], alpha=0.7, color='lavender')
result, resultup, resultdown = get_result_with_CI(params, bmi, eta1, eta220)
save_resultforsimu(result, resultup, resultdown, orgres, orgresup, orgresdown, './simuresult/resultsimu110Jun20.csv', t_range_subdt)
plt.plot(t_range_subdt, result[:, 2] + result[:, 3] + result[:, 4], 'tomato', label = 'No community-based NAT after 20 June')
plt.fill_between(t_range_subdt, resultdown[:, 2] + resultdown[:, 3] + resultdown[:, 4], resultup[:, 2] + resultup[:, 3] + resultup[:, 4], alpha=0.3, color='salmon')
plt.xlabel('Date')
plt.ylabel('Number of cases')
ax.xaxis.set_major_formatter(formatter)
plt.gcf().autofmt_xdate()
# plt.grid()
plt.legend(loc='upper left')

# 234: Localized lockdown only
ax = plt.subplot(234)
result, resultup, resultdown = get_result_with_CI(params, bmi, eta1_noctrl, eta2_noctrl)
save_resultforsimu(result, resultup, resultdown, orgres, orgresup, orgresdown, './simuresult/resultsimu100.csv', t_range_subdt)
plt.plot(t_range_subdt, orgres[:, 2] + orgres[:, 3] + orgres[:, 4], 'b', label = 'Combined interventions')
plt.fill_between(t_range_subdt, orgresdown[:, 2] + orgresdown[:, 3] + orgresdown[:, 4], orgresup[:, 2] + orgresup[:, 3] + orgresup[:, 4], alpha=0.5, color='lightblue')
plt.plot(t_range_subdt, result[:, 2] + result[:, 3] + result[:, 4], 'r', label = 'Localized lockdown only')
plt.fill_between(t_range_subdt, resultdown[:, 2] + resultdown[:, 3] + resultdown[:, 4], resultup[:, 2] + resultup[:, 3] + resultup[:, 4], alpha=0.5, color='pink')
plt.xlabel('Date')
plt.ylabel('Number of cases')
ax.xaxis.set_major_formatter(formatter)
plt.gcf().autofmt_xdate()
# plt.grid()
plt.legend(loc='upper left')

# 235: Close-contact tracing only
ax = plt.subplot(235)
result, resultup, resultdown = get_result_with_CI(params, bmi_noctrl, eta1, eta2_noctrl)
save_resultforsimu(result, resultup, resultdown, orgres, orgresup, orgresdown, './simuresult/resultsimu010.csv', t_range_subdt)
plt.plot(t_range_subdt, orgres[:, 2] + orgres[:, 3] + orgres[:, 4], 'b', label = 'Combined interventions')
plt.fill_between(t_range_subdt, orgresdown[:, 2] + orgresdown[:, 3] + orgresdown[:, 4], orgresup[:, 2] + orgresup[:, 3] + orgresup[:, 4], alpha=0.5, color='lightblue')
plt.plot(t_range_subdt, result[:, 2] + result[:, 3] + result[:, 4], 'orange', label = 'Close-contact tracing only')
plt.fill_between(t_range_subdt, resultdown[:, 2] + resultdown[:, 3] + resultdown[:, 4], resultup[:, 2] + resultup[:, 3] + resultup[:, 4], alpha=0.5, color='wheat')
plt.xlabel('Date')
plt.ylabel('Number of cases')
ax.xaxis.set_major_formatter(formatter)
plt.gcf().autofmt_xdate()
# plt.grid()
plt.legend(loc='upper left')

# 236: community-based NAT only
ax = plt.subplot(236)
result, resultup, resultdown = get_result_with_CI(params, bmi_noctrl, eta1_noctrl, eta2)
save_resultforsimu(result, resultup, resultdown, orgres, orgresup, orgresdown, './simuresult/resultsimu001.csv', t_range_subdt)
plt.plot(t_range_subdt, orgres[:, 2] + orgres[:, 3] + orgres[:, 4], 'b', label = 'Combined interventions')
plt.fill_between(t_range_subdt, orgresdown[:, 2] + orgresdown[:, 3] + orgresdown[:, 4], orgresup[:, 2] + orgresup[:, 3] + orgresup[:, 4], alpha=0.5, color='lightblue')
plt.plot(t_range_subdt, result[:, 2] + result[:, 3] + result[:, 4], 'orchid', label = 'Community-based NAT only')
plt.fill_between(t_range_subdt, resultdown[:, 2] + resultdown[:, 3] + resultdown[:, 4], resultup[:, 2] + resultup[:, 3] + resultup[:, 4], alpha=0.7, color='lavender')
plt.xlabel('Date')
plt.ylabel('Number of cases')
ax.xaxis.set_major_formatter(formatter)
plt.gcf().autofmt_xdate()
plt.legend(loc='upper left')
plt.savefig('./simuresult/intervention_simulation.pdf', bbox_inches='tight')


# In[33]:


# community-based NAT only (before June 21)
result, resultup, resultdown = get_result_with_CI(params, bmi_noctrl, eta1_noctrl, eta220)
save_resultforsimu(result, resultup, resultdown, orgres, orgresup, orgresdown, './simuresult/resultsim001beforejun21.csv', t_range_subdt)


# In[35]:


print('Simulation Finished. Results are saved in ./simuresult/')


# In[ ]:




