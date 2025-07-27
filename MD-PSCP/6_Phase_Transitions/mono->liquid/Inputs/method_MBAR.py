
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
temp = np.array([350,360,370,380,390,400,410,420,430,440,450,460,470,480,490])
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
deltacp = np.zeros(nstates)
cps = np.zeros(nstates)
cpl = np.zeros(nstates)

deltaf_alch = -5.22/4.184/(kb*350.0)

#-----------------------------------------------------------------------------------------

samples_liq = mx.pooledsample()
samples_sol = mx.pooledsample()
for state in range(nstates):

    print(f'state {state}')

    kwargs = {}
    kwargs['T'] = temp[state]
    kwargs['beta'] = 1/(kb*temp[state])  
    
    data = pd.read_csv(f'../../../1_NPT_equilibration/liquid/Outputs/{int(temp[state])}_output.out', sep=" ")
    data['U_pot'] = data['PotEng']
    data['PV'] = 1.0*data['Volume']*conv*Navogadro
    data.drop(index=range(5000), inplace=True)
    data.index = np.arange(0, len(data))
    samples_liq.append(mx.sample(data, 'beta*(U_pot + PV)', acfun='Temp', **kwargs))
    
    rho_a[state] = np.mean(data['Density'])
    H_a[state] = np.mean(data['Enthalpy'])/162.0
    a = (data['U_pot'] - data['E_mol'])
    b = (data['U_pot'] - data['E_mol'] + data['PV'])
    cpl[state] = (np.mean(a*b) - np.mean(a)*np.mean(b))/162.0/(temp[state])**2/kb +  (np.mean(data['PV']*b) - np.mean(data['PV'])*np.mean(b))/162.0/(temp[state])**2/kb

    data = pd.read_csv(f'../../../1_NPT_equilibration/mono/Outputs/{int(temp[state])}_output.out', sep=" ")
    data['U_pot'] = data['PotEng']
    data['PV'] = 1.0*data['Volume']*conv*Navogadro
    data.drop(index=range(5000), inplace=True)
    data.index = np.arange(0, len(data))
    samples_sol.append(mx.sample(data, 'beta*(U_pot + PV)', acfun='Temp', **kwargs))
    
    rho_b[state] = np.mean(data['Density'])
    H_b[state] = np.mean(data['Enthalpy'])/160.0
    a = (data['U_pot'] - data['E_mol'])
    b = (data['U_pot'] - data['E_mol'] + data['PV'])
    cps[state] = (np.mean(a*b) - np.mean(a)*np.mean(b))/160.0/(temp[state])**2/kb +  (np.mean(data['PV']*b) - np.mean(data['PV'])*np.mean(b))/160.0/(temp[state])**2/kb
    
    deltaHp[state] = H_a[state] - H_b[state]
    deltacp[state] = cpl[state] - cps[state]

samples_liq.subsampling(integratedACF=True)
samples_sol.subsampling(integratedACF=True)

mixture_liq = mx.mixture(samples_liq, engine=mx.MBAR())
mixture_sol = mx.mixture(samples_sol, engine=mx.MBAR())

results_liq = mixture_liq.free_energies()
results_sol = mixture_sol.free_energies()

#-----------------------------------------------------------------------------------------

variables = dict(beta=[], T=[])
temp_new = np.linspace(temp.min(), temp.max(), npoints)

for point in range(npoints):

    variables['beta'].append( 1/(kb*temp_new[point])  )
    variables['T'].append( temp_new[point]  )


reweighting_liq = mixture_liq.reweighting(
    potential='beta*(U_pot + PV)',
    conditions=pd.DataFrame(variables)
    )

reweighting_sol = mixture_sol.reweighting(
    potential='beta*(U_pot + PV)',
    conditions=pd.DataFrame(variables)
    )

reweighting_liq['f'] = reweighting_liq['f']/162.0
reweighting_sol['f'] = reweighting_sol['f']/160.0
results_liq['f'] = results_liq['f']/162.0
results_sol['f'] = results_sol['f']/160.0

#-----------------------------------------------------------------------------------------

for point in range(npoints):
    deltaG[point] = (reweighting_sol['f'][point] - reweighting_liq['f'][point] - \
                     results_sol['f'][0] + results_liq['f'][0] + \
                         deltaf_alch)*kb*reweighting_liq['T'][point]
                         
#-----------------------------------------------------------------------------------------

indice = np.where(np.diff(np.sign(deltaG)))[0] 
raiz = (reweighting_liq['T'].values[indice[0]] + reweighting_liq['T'].values[indice[0]+1] )/2
print(raiz)

np.savetxt('data_T_H.txt', (np.array([temp, H_b])).T, delimiter=',')
np.savetxt('data_T_rho.txt', (np.array([temp, rho_b*1000])).T, delimiter=',')
np.savetxt('data_T_deltaG.txt', (np.array([reweighting_liq['T'].values, deltaG*4.184])).T, delimiter=',')
np.savetxt('data_T_deltaGp.txt', (np.array([temp, deltaGp*4.184])).T, delimiter=',')
#-----------------------------------------------------------------------------------------
fig = plt.figure(figsize=(3.0, 4.0),  dpi=300)
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)

ax1.set_ylabel(r'$\beta G - \beta_\text{ref} G_\text{ref}$')
ax2.set_ylabel(r'$\Delta G_\textsc{l,s}$ (kJ/mol)')
ax2.set_xlabel(r'$T$ (K)')
ax1.set_xticklabels([])

p1, = ax1.plot(reweighting_liq['T'], reweighting_liq['f'], 'c-', label='Liquid')
p2, = ax1.plot(results_liq['T'], results_liq['f'], 'co', label='')
p3, = ax1.plot(reweighting_sol['T'], reweighting_sol['f'], 'm-', label='Solid')
p4, = ax1.plot(results_sol['T'], results_sol['f'], 'mv', label='')

from matplotlib.legend_handler import HandlerTuple
ax1.legend([(p1, p2), (p3, p4)], ['Liquid', 'Solid'], 
               handler_map={tuple: HandlerTuple(ndivide=None)}, handlelength=2)

ax2.plot(reweighting_liq['T'], deltaG*4.184, 'k-', label='')
ax2.plot(np.linspace(350, 490, 100), np.linspace(0, 0, 100), 'k--', label='')

fig.align_ylabels([ax1, ax2])
plt.tight_layout()

plt.savefig('figura1.png', dpi=300)
plt.show()

fig = plt.figure(figsize=(3.0, 2.0),  dpi=300)
plt.ylabel(r'$\rho$ (kg/m$^3$)')
plt.xlabel(r'$T$ (K)')
plt.plot(temp, rho_a*1000, 'ko', label=r'liquid')
plt.plot(temp, rho_b*1000, 'ro', label=r'solid')
plt.legend()
plt.savefig('figura2.png', dpi=300)
plt.show()

from scipy.optimize import curve_fit
def linear_fit(x, y):
    def linear_function(x, a, b):
        return a * x + b   
    params, _ = curve_fit(linear_function, x, y)
    return params 

params = linear_fit(temp[8:13], deltaHp[8:13])

fig = plt.figure(figsize=(3.0, 2.0),  dpi=300)
plt.ylabel(r'$\Delta H_{s \rightarrow l}$ (kJ/mol)')
plt.xlabel(r'$T$ (K)')
plt.plot(temp, deltaHp*4.184, 'ko', label='')
plt.plot(temp, (np.array(temp) * params[0] + params[1])*4.184, color='red', label="Ajuste Linear")
plt.savefig('figura3.png', dpi=300)
plt.show()

print(params[0]*4.184*1e3)

fig = plt.figure(figsize=(3.0, 2.0),  dpi=300)
plt.ylabel(r'$\Delta C_{s \rightarrow l}$ (J/K/mol)')
plt.xlabel(r'$T$ (K)')
plt.plot(temp, deltacp*4.184*1e3, 'ko', label='')
plt.plot(temp, cpl*4.184*1e3, 'ro', label='')
plt.plot(temp, cps*4.184*1e3, 'bo', label='')
plt.savefig('figura4.png', dpi=300)
plt.show()



