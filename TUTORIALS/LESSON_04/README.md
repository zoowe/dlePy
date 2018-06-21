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

After calculation is done, remove `tmp` directory.

## :three: Final products
The following items must be delivered to `LESSON_04` folder before moving to `LESSON_05`

:heavy_check_mark: A text file, name `RESULT.dat`, with the following information obtained from `output.dat`:

- Name of the program used.
- Version of the program used.
- Date and time when program started.
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