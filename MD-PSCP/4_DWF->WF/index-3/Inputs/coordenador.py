
import os, sys, subprocess
import numpy as np

restart = 2000200

ai = np.array([32.029, 32.107385, 32.18577, 32.264155, 32.34254, 32.420925, 32.49931, 32.577695, 32.65608, 32.734465, 32.81285, 32.891235, 32.96962, 33.048005, 33.12639, 33.204775, 33.28316, 33.361545, 33.43993, 33.518315, 33.5967])
bi = np.array([32.029, 32.107385, 32.18577, 32.264155, 32.34254, 32.420925, 32.49931, 32.577695, 32.65608, 32.734465, 32.81285, 32.891235, 32.96962, 33.048005, 33.12639, 33.204775, 33.28316, 33.361545, 33.43993, 33.518315, 33.5967])
ci = np.array([32.029, 32.107385, 32.18577, 32.264155, 32.34254, 32.420925, 32.49931, 32.577695, 32.65608, 32.734465, 32.81285, 32.891235, 32.96962, 33.048005, 33.12639, 33.204775, 33.28316, 33.361545, 33.43993, 33.518315, 33.5967])
alphai = np.array([90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0])
betai = np.array([90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0])
gammai = np.array([90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0])

for i in range(21): 
    ft=open(f'0_output_{i}.txt','a') 
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

        os.system("mpirun -np 24 /afs/crc.nd.edu/user/g/gcorrea2/lammps-29Sep2021/src/lmp_mpi -in in.step2.pos")

        f = open("output.out", "w")
        msg  = subprocess.run(['/afs/crc.nd.edu/user/g/gcorrea2/postlammps/postlammps', '-in', 'tcne.log', 'print', 'Step', 'Volume', 'Press', 'KinEng', 'PotEng', 'Enthalpy', 'E_mol', 'E_vdwl', 'E_tail', 'E_coul', 'E_long'], stdout=f)
        f.close()

        fi=open('output.out','r') 
        ft=open(f'0_output_{i}.txt','a') 
    
        fi.readline() #descartando primeira linha
        ft.write(fi.readline())
        fi.close()
        ft.close()
    
    
