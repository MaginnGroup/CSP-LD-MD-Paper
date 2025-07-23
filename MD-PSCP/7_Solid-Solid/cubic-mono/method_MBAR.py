
#-----------------------------------------------------------------------------------------

import mics as mx
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import figstyle

#-----------------------------------------------------------------------------------------

conv = 1e-27*24.21726e-3 #atm A3 --> atm L --> kcal
Navogadro = 6.02214076e23 #atomos/mol
kb = 1.9872e-3 #kcal/mol/K
temp = np.array([110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300,310,320,330,340,350,360,370,380,390,400,410,420,430,440,450,460,470,480,490])
nstates = len(temp)
deltaSp = np.zeros(nstates)
deltaHp = np.zeros(nstates)
deltaGp = np.zeros(nstates)
dev_deltaGp = np.zeros(nstates)
rho_a = np.zeros(nstates)
rho_b = np.zeros(nstates)
H_a = np.zeros(nstates)
H_b = np.zeros(nstates)

deltaf_alch = -0.24/4.184/(kb*350.0)
dev = 0.09/4.184/(kb*350.0)

#-----------------------------------------------------------------------------------------

samples_a = mx.pooledsample()
samples_b = mx.pooledsample()
mx.verbose = True
for state in range(nstates):

    print(f'state {state}')

    kwargs = {}
    kwargs['T'] = temp[state]
    kwargs['beta'] = 1/(kb*temp[state])  
    
    dataa = pd.read_csv(f'../cubic/{int(temp[state])}/output.out', sep=" ")
    dataa['U_pot'] = dataa['PotEng']
    dataa['PV'] = 1.0*dataa['Volume']*conv*Navogadro
    dataa.drop(index=range(5000), inplace=True)
    
    dataa.index = np.arange(0, len(dataa))
    samples_a.append(mx.sample(dataa, 'beta*(U_pot + PV)', acfun='Temp', **kwargs))

    datab = pd.read_csv(f'../mono/{int(temp[state])}/output.out', sep=" ")
    datab['U_pot'] = datab['PotEng']
    datab['PV'] = 1.0*datab['Volume']*conv*Navogadro
    datab.drop(index=range(5000), inplace=True)

    datab.index = np.arange(0, len(datab))
    samples_b.append(mx.sample(datab, 'beta*(U_pot + PV)', acfun='Temp', **kwargs))
    
    H_a[state] = np.mean(dataa['Enthalpy'])/162.0
    H_b[state] = np.mean(datab['Enthalpy'])/160.0
    rho_a[state] = np.mean(dataa['Density'])
    rho_b[state] = np.mean(datab['Density'])

samples_a.subsampling(integratedACF=True)
samples_b.subsampling(integratedACF=True)

mixture_a = mx.mixture(samples_a, engine=mx.MBAR())
mixture_b = mx.mixture(samples_b, engine=mx.MBAR())

results_a = mixture_a.free_energies()
results_b = mixture_b.free_energies()

#-----------------------------------------------------------------------------------------

variables = dict(beta=[], T=[])

for temps in temp:

    variables['beta'].append( 1/(kb*temps)  )
    variables['T'].append( temps  )

properties = dict(
    Hp='Enthalpy'
    )

reweighting_a = mixture_a.reweighting(
    potential='beta*(U_pot + PV)',
    properties=properties,
    conditions=pd.DataFrame(variables)
    )

reweighting_b = mixture_b.reweighting(
    potential='beta*(U_pot + PV)',
    properties=properties,
    conditions=pd.DataFrame(variables)
    )

reweighting_a['Hp'] = reweighting_a['Hp']/162.0
reweighting_b['Hp'] = reweighting_b['Hp']/160.0
reweighting_a['dHp'] = reweighting_a['dHp']/162.0
reweighting_b['dHp'] = reweighting_b['dHp']/160.0
results_a['f'] = results_a['f']/162.0
results_b['f'] = results_b['f']/160.0
results_a['df'] = results_a['df']/162.0
results_b['df'] = results_b['df']/160.0

#-----------------------------------------------------------------------------------------

for state in range(nstates):
    deltaGp[state] = (results_b['f'][state] - results_a['f'][state] - \
                      results_b['f'][24] + results_a['f'][24] + \
                      deltaf_alch)*kb*temp[state]
    dev_deltaGp[state] = np.sqrt( (results_b['df'][state])**2 + (results_a['df'][state])**2 + \
                     (results_b['df'][24])**2 + (results_a['df'][24])**2 + \
                         dev**2)*kb*temp[state]


lala = (reweighting_b['Hp']-reweighting_a['Hp']).values
dev_deltaHp = 2*np.sqrt(reweighting_b['dHp']**2+reweighting_a['dHp']**2).values
from scipy.ndimage import gaussian_filter1d
deltaHp = gaussian_filter1d(lala,sigma=2)
deltaSp = (deltaHp - deltaGp)/temp
dev_deltaSp = np.sqrt(dev_deltaHp**2 + dev_deltaGp**2)/temp
#-----------------------------------------------------------------------------------------

np.savetxt('data_T_H_cubic.txt', (np.array([temp, H_a*4.184])).T, delimiter=',')
np.savetxt('data_T_rho_cubic.txt', (np.array([temp, rho_a])).T, delimiter=',')
np.savetxt('data_T_H.txt', (np.array([temp, H_b*4.184])).T, delimiter=',')
np.savetxt('data_T_rho.txt', (np.array([temp, rho_b])).T, delimiter=',')
np.savetxt('data_T_deltaGp.txt', (np.array([temp, deltaGp*4.184, dev_deltaGp*4.184])).T, delimiter=',')
np.savetxt('data_T_deltaHp.txt', (np.array([temp, deltaHp*4.184, dev_deltaHp*4.184])).T, delimiter=',')
np.savetxt('data_T_TdeltaSp.txt', (np.array([temp, deltaSp*temp*4.184, dev_deltaSp*temp*4.184])).T, delimiter=',')

#-----------------------------------------------------------------------------------------
fig = plt.figure(figsize=(3.0, 4.0),  dpi=300)
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)

ax1.set_ylabel(r'$\beta G - \beta_\text{ref} G_\text{ref}$')
ax2.set_ylabel(r'$\Delta G$ (kJ/mol)')
ax2.set_xlabel(r'$T$ (K)')
ax1.set_xticklabels([])

ax1.plot(results_a['T'], results_a['f'], 'co', label='')
ax1.plot(results_b['T'], results_b['f'], 'mv', label='')

ax2.plot(temp, deltaGp*4.184, 'ko', label='')
ax2.fill_between(temp.flatten(), (deltaGp*4.184 - dev_deltaGp*4.184).flatten(), (deltaGp*4.184 + dev_deltaGp*4.184).flatten(), color='black', alpha=0.2)

fig.align_ylabels([ax1, ax2])
plt.tight_layout()

plt.savefig('figura_mbar_1.png', dpi=300)
plt.show()

fig = plt.figure(figsize=(3.0, 2.0),  dpi=300)
plt.ylabel(r'$\Delta H$ (kJ/mol)')
plt.xlabel(r'$T$ (K)')
plt.plot(temp, deltaHp*4.184, 'ko-', label='')
plt.fill_between(temp.flatten(), ((deltaHp - dev_deltaHp)*4.184).flatten(), ((deltaHp + dev_deltaHp)*4.184).flatten(), color='black', alpha=0.2)
plt.ylim(-2,10)
plt.savefig('figura_mbar_2.png', dpi=300)
plt.show()

fig = plt.figure(figsize=(3.0, 2.0),  dpi=300)
plt.ylabel(r'$\rho$ (kg/m$³$)')
plt.xlabel(r'$T$ (K)')
plt.plot(temp, rho_b*1000, 'ko', label='')
plt.savefig('figura_mbar_3.png', dpi=300)
plt.show()

fig = plt.figure(figsize=(3.0, 2.0),  dpi=300)
plt.ylabel(r'$-T\Delta S$ (kJ/mol)')
plt.xlabel(r'$T$ (K)')
plt.plot(temp, -deltaSp*temp*4.184, 'ko-', label='')
plt.fill_between(temp.flatten(), (-(deltaSp - dev_deltaSp)*temp*4.184).flatten(), (-(deltaSp + dev_deltaSp)*temp*4.184).flatten(), color='black', alpha=0.2)
plt.ylim(-8,2)
plt.savefig('figura_mbar_6.png', dpi=300)
plt.show()



