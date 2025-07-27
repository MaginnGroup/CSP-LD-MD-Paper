
import numpy as np
import os

pwd=os.getcwd()

lbda = np.array([0.0, 0.025, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0])
n=np.size(lbda)

with open('in.step3_model', 'r') as file :
  file_data = file.read()

for i in range(0,n):
    os.system(f"mkdir {i}")
    os.chdir(pwd+f"/{i}")
    os.system("cp ../sub.lmp sub.lmp")
    file_fixed = file_data.replace('xst1', f'{lbda[i]}')
    with open('in.step3', 'w') as file:
        file.write(file_fixed)
    os.chdir(pwd)
    
   
   
