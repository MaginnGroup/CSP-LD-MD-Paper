
#-----------------------------------------------------------------------------------------
# PACOTES

import mics as mx
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import figstyle

#-----------------------------------------------------------------------------------------
# ENTRADAS

T = 350.0 #temperatura
kb = 1.9872e-3 #kcal/mol/K
conv = 1e-27*24.21726e-3 #atm A3 --> atm L --> kcal
Navogadro = 6.02214076e23 #atomos/mol
nstates = 21
kT = T*kb

# -5.400710251047948 0.0006876612244249315 kJ/mol

#-----------------------------------------------------------------------------------------
# POTENCIAL E ENERGIA LIVRE DOS ESTADOS AMOSTRADOS
Pressure = np.zeros(nstates)
Volume = np.zeros(nstates)
samples = mx.pooledsample()
for state in range(21):
    print(f'state {state}')
    data = pd.read_csv(f'{state}_output_{state}.txt', sep=" ")   
    for i in range(21):
        data_new = pd.read_csv(f'{state}_output_{i}.txt', sep=" ") 
        data[f'Upot_{i}'] = data_new['PotEng'] 
    data.index = np.arange(0, len(data))
    Volume[state] = data['Volume'][0]
    kwargs = {'beta': 1.0/kT}
    kwargs['N'] = 162.0
    kwargs['V'] = Volume[state]
    samples.append(mx.sample(data, f'-N*ln(V)+beta*(Upot_{state})', acfun='PotEng', **kwargs)) 
    Pressure[state] = np.mean(data['Press'].values)
    
samples.subsampling(integratedACF=True)

mixture = mx.mixture(samples, engine=mx.MBAR())

results = mixture.free_energies()
results['F'] = results['f']*kT
results['dF'] = results['df']*kT

#-----------------------------------------------------------------------------------------
# REWEIGHTING

reweighting = pd.DataFrame()

for state in range(nstates):
    variables = dict(V=[])
    variables['V'].append(Volume[state])
    properties = dict(
        P = 'Press'
        ) 
    rw = mixture.reweighting(
        potential=f'-N*ln(V)+beta*(Upot_{state})',
        properties=properties,
        conditions=pd.DataFrame(variables), 
        beta=1/kT,
        N=162.0
        )
    reweighting = reweighting._append(rw, True)

reweighting['F'] = reweighting['f']*kT
reweighting['dF'] = reweighting['df']*kT
#results['Pressure'] = Pressure
#reweighting.to_csv('data_step2_reweighting.out', index=None, sep=' ', mode='w')
#results.to_csv('data_step2.out', index=None, sep=' ', mode='w')

#-----------------------------------------------------------------------------------------
# ENERGIA LIVRE FINAL

Nat = 162.0

deltaG1 = results['F'].iloc[-1]/Nat*4.184
d_deltaG1 = results['dF'].iloc[-1]/Nat*4.184

print(f'Delta G = {deltaG1} {d_deltaG1} kJ/mol')

np.savetxt('data1.txt', (np.array([results['V'], results['F']/Nat*4.184])).T, delimiter=',')
np.savetxt('data2.txt', (np.array([reweighting['V'], reweighting['F']/Nat*4.184])).T, delimiter=',')

#-----------------------------------------------------------------------------------------
# PLOTAR GRÁFICOS

fig, ax = plt.subplots(2, 1, figsize=(2.5, 5.0),  dpi=300)
ax[0].plot(results['V'], results['F']/Nat*4.184, 'k-', label='')
ax[0].plot(results['V'], results['F']/Nat*4.184, 'bo', label='')
ax[0].legend()
ax[1].plot(reweighting['V'], reweighting['P'], 'k-', label='')
ax[1].plot(Volume, Pressure, 'bo', label='')
ax[0].set_ylabel('$A$ (kJ/mol)')
ax[1].set_ylabel('$P$ (atm)')
ax[1].set_xlabel('$V$ (\AA$^3$)')
ax[0].set_xticklabels([])
plt.savefig('figure_MBAR_total.png', dpi=300)
plt.show()


