# PWSCF: First run

## :one: Purposes
Performing the very first PWSCF calculation and extracting results.

## :two: Hows
### :large_blue_diamond: Performing a scf calculation
In this folder, there are two files:
```
input.inp
job
```
- Inspect the files
- Submit the `job` file to run this calculation
```
qsub job
```
**Note that, this tutorial was designed for `challenger`. If you are using other machine, you do need to modify `input.inp` and `job` files, accordingly**

After calculation is done, remove `tmp` directory.

## :three: Final products
The following items must be delivered to `LESSON_04` folder before moving to `LESSON_05`

:heavy_check_mark: A text file, name `RESULT.dat`, with the following information obtained from `output.dat`:

- Name of the program used.
- Version of the program used.
- Date and time when program started.
- Number of processors used. 
- unit-cell volume         
- number of atoms   
- number of atomic types   
- number of electrons      
- number of bands 
- kinetic-energy cutoff    
- charge density cutoff    
- convergence threshold    
- number of k-points
- Potential file   
- Number of valence electron per atom
- Number of iterations
- Total energy
- Date and time when program stoped.
- Total running time (Walltime).
