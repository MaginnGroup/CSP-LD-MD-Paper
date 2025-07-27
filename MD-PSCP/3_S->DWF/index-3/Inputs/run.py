
import numpy as np
import os

pwd=os.getcwd()

lbda = np.array([0.0001, 0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275, 0.3, 0.325, 0.35, 0.375, 0.4, 0.425, 0.45, 0.475, 0.5, 0.525, 0.55, 0.575, 0.6, 0.625, 0.65, 0.675, 0.7, 0.725, 0.75, 0.775, 0.8, 0.825, 0.85, 0.875, 0.9, 0.925, 0.95, 0.975, 0.9999])
n=len(lbda)

with open('in.step1_model', 'r') as file :
  file_data = file.read()

for i in range(0,n):
    os.system(f"mkdir {i}")
    os.chdir(pwd+f"/{i}")
    os.system("cp ../sub.lmp sub.lmp")
    file_fixed = file_data.replace('xst1', f'{lbda[i]}')
    with open('in.step1', 'w') as file:
        file.write(file_fixed)
    os.chdir(pwd)
    
   
    
