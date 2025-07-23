import os, sys, subprocess
import numpy as np

restart = 2800200

ai = np.array([30.025, 30.0889, 30.1528, 30.2167, 30.2806, 30.3445, 30.4084, 30.4723, 30.5362, 30.6001, 30.664, 30.7279, 30.7918, 30.8557, 30.9196, 30.9835, 31.0474, 31.1113, 31.1752, 31.2391, 31.303])
bi = np.array([30.026, 30.08985, 30.1537, 30.21755, 30.2814, 30.34525, 30.4091, 30.47295, 30.5368, 30.60065, 30.6645, 30.72835, 30.7922, 30.85605, 30.9199, 30.98375, 31.0476, 31.11145, 31.1753, 31.23915, 31.303])
ci = np.array([30.026, 30.08985, 30.1537, 30.21755, 30.2814, 30.34525, 30.4091, 30.47295, 30.5368, 30.60065, 30.6645, 30.72835, 30.7922, 30.85605, 30.9199, 30.98375, 31.0476, 31.11145, 31.1753, 31.23915, 31.303])
alphai = np.array([90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0])
betai = np.array([90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0])
gammai = np.array([90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0])

for i in range(21): 
    ft=open(f'16_output_{i}.txt','a') 
    ft.write('Step Volume Press KinEng PotEng Enthalpy E_mol E_vdwl E_tail E_coul E_long \n')
    ft.close()

for j in range(10000):
    for i in range(21): 

        template ='''
variable        NAME index tcne
log             ${NAME}.log

variable	 a equal $variavel1
variable	 b equal $variavel2
variable	 c equal $variavel3
variable	 alpha equal $variavel4
variable	 beta equal $variavel5
variable	 gamma equal $variavel6

variable        xv equal v_a
variable        xyv equal v_b*cos(v_gamma*PI/180)
variable        yv equal sqrt(v_b*v_b-v_xyv*v_xyv)
variable        xzv equal v_c*cos(v_beta*PI/180)
variable        yzv equal (v_b*v_c*cos(v_alpha*PI/180)-v_xyv*v_xzv)/v_yv
variable        zv equal sqrt(v_c*v_c-v_xzv*v_xzv-v_yzv*v_yzv)

units           real
special_bonds   lj 0.0 0.0 0.5 coul 0.0 0.0 0.5
atom_style      full
pair_style	hybrid/overlay lj/cut 12.0 coul/long 12.0 

bond_style	 harmonic
dihedral_style	 opls
angle_style	 harmonic

read_restart    restart.tcne.2.$variavel7
''' 
        template_mod = template*1
        template_mod = template_mod.replace('$variavel1',f'{ai[i]}')
        template_mod = template_mod.replace('$variavel2',f'{bi[i]}')
        template_mod = template_mod.replace('$variavel3',f'{ci[i]}')
        template_mod = template_mod.replace('$variavel4',f'{alphai[i]}')
        template_mod = template_mod.replace('$variavel5',f'{betai[i]}')
        template_mod = template_mod.replace('$variavel6',f'{gammai[i]}')
        template_mod = template_mod.replace('$variavel7',f'{int(restart) + j*200}')

        f=open('in.input','w')
        f.write(template_mod)
        f.close()

        msg  = subprocess.run(['mpirun', '-np', '24', '/afs/crc.nd.edu/user/g/gcorrea2/lammps-29Sep2021/src/lmp_mpi', '-in', 'in.step2.pos'],capture_output=True)

        f = open("output.out", "w")
        msg  = subprocess.run(['/afs/crc.nd.edu/user/g/gcorrea2/postlammps/postlammps', '-in', 'tcne.log', 'print', 'Step', 'Volume', 'Press', 'KinEng', 'PotEng', 'Enthalpy', 'E_mol', 'E_vdwl', 'E_tail', 'E_coul', 'E_long'], stdout=f)
        f.close()

        fi=open('output.out','r') 
        ft=open(f'16_output_{i}.txt','a') 
    
        fi.readline() #descartando primeira linha
        ft.write(fi.readline())
        fi.close()
        ft.close()
    
    
