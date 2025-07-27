
#-----------------------------------------------------------------------------------------

import mics as mx
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

#-----------------------------------------------------------------------------------------
# ENTRADAS

conv = 1e-27*24.21726e-3 #atm A3 --> atm L --> kcal
Navogadro = 6.02214076e23 #atomos/mol
kb = 1.9872e-3 #kcal/mol/K
temp = np.array([110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300,310,320,330,340,350,360,370,380,390,400,410,420,430,440,450,460,470,480,490])
nstates = len(temp)
deltaG = np.zeros(npoints)
deltaG_upper = np.zeros(npoints)
deltaG_lower = np.zeros(npoints)
rho_a = np.zeros(nstates)
rho_b = np.zeros(nstates)
H_a = np.zeros(nstates)
H_b = np.zeros(nstates)
deltaSp = np.zeros(nstates)
deltaHp = np.zeros(nstates)
deltaGp = np.zeros(nstates)
dev_deltaGp = np.zeros(nstates)

deltaf_alch = 2.33/4.184/(kb*350.0)
dev = 0.09/4.184/(kb*350.0)

#-----------------------------------------------------------------------------------------

samples_liq = mx.pooledsample()
samples_sol = mx.pooledsample()
for state in range(nstates):

    print(f'state {state}')

    kwargs = {}
    kwargs['T'] = temp[state]
    kwargs['beta'] = 1/(kb*temp[state])  
    
    data = pd.read_csv(f'../../../1_NPT_equilibration/cubic/Outputs/{int(temp[state])}_output.out', sep=" ")
    data['U_pot'] = data['PotEng']
    data['PV'] = 1.0*data['Volume']*conv*Navogadro
    data.drop(index=range(5000), inplace=True)
    data.index = np.arange(0, len(data))
    samples_liq.append(mx.sample(data, 'beta*(U_pot + PV)', acfun='Temp', **kwargs))
    
    rho_a[state] = np.mean(data['Density'])
    H_a[state] = np.mean(data['Enthalpy'])/162.0

    data = pd.read_csv(f'../../../1_NPT_equilibration/index-9/Outputs/{int(temp[state])}_output.out', sep=" ")
    data['U_pot'] = data['PotEng']
    data['PV'] = 1.0*data['Volume']*conv*Navogadro
    data.drop(index=range(5000), inplace=True)
    data.index = np.arange(0, len(data))
    samples_sol.append(mx.sample(data, 'beta*(U_pot + PV)', acfun='Temp', **kwargs))
    
    rho_b[state] = np.mean(data['Density'])
    H_b[state] = np.mean(data['Enthalpy'])/180.0
    

samples_liq.subsampling(integratedACF=True)
samples_sol.subsampling(integratedACF=True)

mixture_liq = mx.mixture(samples_liq, engine=mx.MBAR())
mixture_sol = mx.mixture(samples_sol, engine=mx.MBAR())

results_liq = mixture_liq.free_energies()
results_sol = mixture_sol.free_energies()

#-----------------------------------------------------------------------------------------

variables = dict(beta=[], T=[])

for temps in temp:

    variables['beta'].append( 1/(kb*temps)  )
    variables['T'].append( temps  )

properties = dict(
    Hp='Enthalpy'
    )

reweighting_liq = mixture_liq.reweighting(
    potential='beta*(U_pot + PV)',
    properties=properties,
    conditions=pd.DataFrame(variables)
    )

reweighting_sol = mixture_sol.reweighting(
    potential='beta*(U_pot + PV)',
    properties=properties,
    conditions=pd.DataFrame(variables)
    )

reweighting_liq['Hp'] = reweighting_liq['Hp']/162.0
reweighting_sol['Hp'] = reweighting_sol['Hp']/180.0
reweighting_liq['dHp'] = reweighting_liq['dHp']/162.0
reweighting_sol['dHp'] = reweighting_sol['dHp']/180.0
results_liq['f'] = results_liq['f']/162.0
results_sol['f'] = results_sol['f']/180.0
results_liq['df'] = results_liq['df']/162.0
results_sol['df'] = results_sol['df']/180.0

#-----------------------------------------------------------------------------------------

for point in range(nstates):
    deltaGp[point] = (results_sol['f'][point] - results_liq['f'][point] - \
                     results_sol['f'][24] + results_liq['f'][24] + \
                         deltaf_alch)*kb*results_liq['T'][point]
    dev_deltaGp[point] = np.sqrt( (results_sol['df'][point])**2 + (results_liq['df'][point])**2 + \
                     (results_sol['df'][24])**2 + (results_liq['df'][24])**2 + \
                         dev**2)*kb*results_liq['T'][point]

lala = (reweighting_sol['Hp']-reweighting_liq['Hp']).values
dev_deltaHp = 2*np.sqrt(reweighting_sol['dHp']**2+reweighting_liq['dHp']**2).values
from scipy.ndimage import gaussian_filter1d
deltaHp = gaussian_filter1d(lala,sigma=2)
deltaSp = (deltaHp - deltaGp)/temp
dev_deltaSp = np.sqrt(dev_deltaHp**2 + dev_deltaGp**2)/temp

#-----------------------------------------------------------------------------------------

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
ax2.set_ylabel(r'$\Delta G_\textsc{l,s}$ (kJ/mol)')
ax2.set_xlabel(r'$T$ (K)')
ax1.set_xticklabels([])

p2, = ax1.plot(results_liq['T'], results_liq['f'], 'co', label='')
p4, = ax1.plot(results_sol['T'], results_sol['f'], 'mv', label='')

ax2.plot(temp, deltaGp*4.184, 'ko', label='')
ax2.fill_between(temp.flatten(), (deltaGp*4.184 - dev_deltaGp*4.184).flatten(), (deltaGp*4.184 + dev_deltaGp*4.184).flatten(), color='black', alpha=0.2)

fig.align_ylabels([ax1, ax2])
plt.tight_layout()

plt.savefig('figura_mbar.png', dpi=300)
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


