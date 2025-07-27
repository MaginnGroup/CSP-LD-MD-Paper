# Electronic Supporting File Repository 

Authors: G. B. Correa, S. Konstantinopoulos, B. I. Tan, Y. Zhang, F. W. Tavares, C. S. Adjiman, E. J. Maginn

This repository contains supporting files for the article 'Assessing Polymorph Stability and Phase Transitions at Finite Temperature: Integrating Crystal Structure Prediction, Lattice Dynamics, and Molecular Dynamics' (DOI: **[*INSERT LATER*]**).

In this work, we propose a workflow to efficiently assess the phase behaviour of organic molecular crystals at finite temperature. The proposed methodology combines zeroth-order Crystal Structure Prediction (CSP), harmonic approximation lattice dynamics (HA-LD), and molecular dynamics based on the pseudo-supercritical path (PSCP) method. 

The assocaited input (and output) files for each of these stages are compiled in this repository and described below. 

## CSP

Zeroth-order CSP is performed to generate potentially observeable crystal structures. This entails two steps: (1) a global search using CrystalPredictor [1-5], followed by (2) a refinement using CrystalOptimizer [6]. 

For each step, the following inputs/outputs have been provided:

### 1_Global_Search
#### Input Files
1. **input.in**: Main CrystalPredictor input file
2. **potential.in**: Repulsion-dispersion (rd) parameter file
3. **rigid_lam_intra**: LAM file with gas-phase minimized geometry and HLYGAt atomic-charges. The target molecule is treated as rigid at this step, so no additional LAMs were generated. 

#### Output Files
These outputs correspond to the clustered global search landscape shown in Figure 4a of the article. 
1. **N/N.res**: Optimized crystal structure corresponding to cluster N, in SHELX format
2. **N/Cluster_details**: Corresponding lattice energy and density for cluster N
3. **All_Energy_Density_CrystPred.dat**: Summary of lattice energies and densities for entire landscape

### 2_Refinement
#### Input Files
1. **CrystOpt.input**: Main CrystalOptimizer input file
2. **bondlengths**: Standard bond lengths file
3. **dmarel.axis**: DMACRYS [7] molecular axis file
4. **pote.dat**: Repulsion-dispersion (rd) parameter file
5. **mol.input**: Molecular definition file
6. **dmaSCF_2**: Input settings for the GDMA program [8,9]
7. **structures.res**: A target crystal structure in SHELX format should be provided to the program as well (not given here)

#### Notes on CrystalOptimizer
1. In **CrystOpt.input** and **mol.input** files, "PATH" should be replaced by local paths
2. CrystalOptimizer has external dependecies on:
- DMACRYS [7]
- GDMA [8,9]

#### Output Files
These outputs correspond to the clustered refinement landscape shown in Figure 4b of the article.
1. **N/N.res**: Optimized crystal structure corresponding to cluster N, in SHELX format
2. **N/summary_CrystalOptimizer.out**: Corresponding lattice energy and density for cluster N
3. **All_Energy_Density_CrystOpt.dat**: Summary of lattice energies and densities for entire landscape

## HA-LD
Lattice dynamics with the harmonic approximation is performed using the CrystalDynamics code [10,11]. This acts as a finite temperature prescreening.

The following inputs/outputs have been provided:

#### Input Files
1. **vibrations.in**: Main CrystalDynamics input file
2. **dmaSCF_2**: Input settings for the GDMA program [8,9]
3. **structures.res**: A target crystal structure in SHELX format should be provided to the program as well. In this case, these correspond to index 1-100 from the CSP refinement step. 

#### Output Files
These outputs correspond to the 100 structures subject to HA-LD in Figure 5 of the article.
1. **N/summary_CrystalDynamics.out**: Relative free energy, with respect to cluster 1, for cluster N between 0-400 K 

## MD-PSCP
Molecular dynamics based on the PSCP method [12,13] is employed as a final assessment of the finite temperature phase behavior. Each folder contains simulation data for multiple polymorphs (cubic, mono, index-3, index-9, index-13).

The following folders have been provided:

### 1_NPT_equilibration
Equilibrate solid and liquid phases at target pressure and temperature (NPT ensemble).
#### Input Files
1. **box.lmp**: Starting configuration and molecular topology (LAMMPS format).
2. **in.solid**: LAMMPS input script.
#### Output Files
1. **_output.out**: LAMMPS outputs for each temperature.
   
### 2_Gaussian_Potential_Fit
Fit Gaussian potential for intermediate restrained states.
#### Input Files
1. **fit.py**, **positions.py**, **probability.py**: Python scripts for potential fitting.
2.  **data.step0.initial**: Starting configuration and molecular topology (LAMMPS format).
3.  **in.gauss_wells**: LAMMPS input script.
#### Output Files
1. **positions_.txt**: Equilibrium coordinates of atoms of types.
2. **dados_dr0_P.txt**: Probability distribution used to fit the parameter k of the Gaussian potential.

### 3_S->DWF
Thermodynamic transformation from the fully interacting solid (S) to the dense weak fluid (DWF).
#### Input Files
1. **run.py**: Python script to automate the MD runs.
2. **data.step1.initial**, **data.out**: Starting configuration and molecular topology (LAMMPS format).
3. **in.step1_model**: LAMMPS input script.
#### Output Files
1. **_output.out**: Output files from MD runs for different λ windows.
2. **method_MBAR.py**: Python script used to compute free energy differences using MBAR.
3. **figure_MBAR_total.png**: Visualization of the MBAR integration for this step.

### 4_DWF->WF
Thermodynamic transformation from the dense weak fluid (DWF) to the weal fluid (WF).
#### Input Files
1. **coordenador.py**: Python script to automate the MD runs.
2.  **data.initial**, **data2.out**: Starting configuration and molecular topology (LAMMPS format).
3. **in.step2.1**, **in.step2.2**, **in.step2.pos**: LAMMPS input scripts.
#### Output Files
1. **_output_.out**: Output files from MD runs for different volume windows.
2. **method_MBAR.py**: Python script used to compute free energy differences using MBAR.
3. **figure_MBAR_total.png**: Visualization of the MBAR integration for this step.

### 5_WF->L

#### Input Files

#### Output Files

### 6_Phase_Transitions

#### Input Files

#### Output Files

## References
[1] Karamertzanis, P. G.; Pantelides, C. C. Ab initio crystal structure prediction—I. Rigid molecules. *Journal of Computational Chemistry* **2005**, *26*, 304-324, DOI: https://doi.org/10.1002/jcc.20165. 

[2] Karamertzanis, P. G.; Pantelides, C. C. Ab initio crystal structure prediction. II. Flexible molecules. *Molecular Physics* **2007**, *105*, 273-291, DOI:
10.1080/00268970601143317.  

[3] Habgood, M.; Sugden, I. J.; Kazantsev, A. V.; Adjiman, C. S.; Pantelides, C. C.
Efficient Handling of Molecular Flexibility in Ab Initio Generation of Crystal Structures. *Journal of Chemical Theory and Computation* **2015**, *11*, 1957-1969, DOI: 10.1021/ct500621v.

[4] Sugden, I.; Adjiman, C. S.; Pantelides, C. C. Accurate and efficient representation
of intramolecular energy in ab initio generation of crystal structures. I. Adaptive local approximate models. *Acta Crystallographica Section B* **2016**, *72*, 864-874, DOI: 10.1107/S2052520616015122.  

[5] Sugden, I. J.; Adjiman, C. S.; Pantelides, C. C. Accurate and efficient representation of intramolecular energy in ab initio generation of crystal structures. II. Smoothed intramolecular potentials. *Acta Crystallographica Section B* **2019**, *75*, 423-433, DOI: 10.1107/S2052520619005778.

[6] Kazantsev, A. V.; Karamertzanis, P. G.; Adjiman, C. S.; Pantelides, C. C. Efficient Handling of Molecular Flexibility in Lattice Energy Minimization of Organic
Crystals. *Journal of Chemical Theory and Computation* **2011**, *7*, 1998–2016, DOI: 10.1021/ct100597e.

[7] Price, S. L.; Leslie, M.; Welch, G. W. A.; Habgood, M.; Price, L. S.; Karamertzanis, P. G.; Day, G. M. Modelling organic crystal structures using distributed multipole
and polarizability-based model intermolecular potentials. *Phys. Chem. Chem. Phys.*
**2010**, *12*, 8478-8490, DOI: 10.1039/C004164E.  

[8] Stone, A. J.; Alderton, M. Distributed multipole analysis Methods and applications.
*Molecular Physics* **2002**, *100*, 221-233, DOI: 10.1080/00268970110089432.

[9] Stone, A. J. Distributed Multipole Analysis: Stability for Large Basis Sets. *Journal of Chemical Theory and Computation* **2005**, *1*, 1128-1132, DOI: 10.1021/ct050190+.

[10] Vasileiadis, M. Calculation of the free energy of crystalline solids. Ph.D. thesis, Imperial College London, 2013.

[11] Konstantinopoulos, S. Free Energy calculations for Crystal Structure Prediction studies. Ph.D. thesis, Imperial College London, 2023.

[12] Eike, D. M.; Brennecke, J. F.; Maginn, E. J. Toward a robust and general molecular
simulation method for computing solid-liquid coexistence. *The Journal of chemical
physics* **2005**, *122*.

[13] Eike, D. M.; Maginn, E. J. Atomistic simulation of solid-liquid coexistence for molecular systems: Application to triazole and benzene. *The Journal of chemical physics* **2006**, *124*.

[14] Correa, G. B.; Zhang, Y.; Abreu, C. R. A.; Tavares, F. W.; Maginn, E. J. Revisiting the pseudo-supercritical path method: An improved formulation for the alchemical calculation of solid-liquid coexistence. *J. Chem. Phys.* **2023**, *159*, 104105, DOI: 10.1063/5.0163564.
