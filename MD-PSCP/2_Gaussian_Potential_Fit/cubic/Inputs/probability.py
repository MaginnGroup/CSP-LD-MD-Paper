
import numpy as np
import matplotlib.pyplot as plt
from numba import njit

freq = 200 
nsteps_total = 1000000 
N = 162 
d = 1000
nframes = int(nsteps_total/freq)

rv = np.zeros((nframes, N, 3), dtype=np.float_)
drv = np.zeros((nframes, N, 3), dtype=np.float_)
rv_mean = np.zeros((N), dtype=np.float_)
P = np.zeros((d), dtype=np.float_)
dr = np.zeros((nframes, N), dtype=np.float_)

f = open("CMTCNE.out", "r")

for n in range(3):
    f.readline()

for frame in range(nframes):

    f.readline()

    for n in range(N):
        
        atom, rx, ry, rz = f.readline().split()       
        atom = int(atom) - 1
        rv[frame, atom, :] = float(rx), float(ry), float(rz)

f.close()

rv_mean = np.mean(rv, axis=0)

for frame in range(nframes):
    drv[frame, :, :] = rv[frame, :, :] - rv_mean
    dr[frame, :] = np.sqrt((rv[frame, :, 0] - rv_mean[:, 0])**2 + (rv[frame, :, 1] - rv_mean[:, 1])**2 + (rv[frame, :, 2] - rv_mean[:, 2])**2)


delta = np.max(dr)/float(d)

dr0 = np.arange(0, d)*delta + delta/2.0

dr_range = np.arange(0, d+1)*delta

@njit
def loop(nframes, N, d, dr, dr_range, P):
    for frame in range(nframes):
        for n in range(N):
            for i in range(d):
                if ((dr[frame, n] < dr_range[i+1]) and (dr[frame, n] >= dr_range[i])):
                    P[i] += 1.0
                    break
    return P

P = loop(nframes, N, d, dr, dr_range, P)  


f = open('dados_dr0_P.txt', 'w')
for i in range(d):
    f.write('%f \t %f\n' % (dr0[i], P[i]/(delta*np.sum(P))))
f.close()

fig = plt.figure(figsize=(2.5, 2.5), dpi=300)
plt.plot(dr0, P/(delta*np.sum(P)), 'k-', label='')
#plt.xticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4])
plt.ylabel('$P$')
plt.xlabel(r'$r$ (A)')
plt.savefig('probability.png', dpi=300)
plt.show()

