# PWSCF: K-point sampling  

## :one: Purposes
Performing convergence test for k-point sampling. Generally, the denser k-point sampling is the more accurate calculation is but computational cost will increase as well. Similar to kinetic energy cutoff determination (`LESSON_05`), we want to find the best possible k-point sampling mesh for accurate results but without extensive computational time.


## :two: Hows
### :large_blue_diamond: Performing scf calculation for a series of k-point mesh.

- In this step, a series of scf calculations will be performed with k-point mesh (`k x k x k`) with `k` is the number of grid point, ranging from 4 to 18.

- Each calculation must be done inside a folder named `K-k`, where `k` is the number of grid point mentioned above.

- You must copy `input.inp` and `job` files from `LESSON_04` to each `K-k` folder, then edit `input.inp` file with new value of Kinetic Energy Cut-off for wavefunctions and for charge density obtained from `LESSON_05` and new value for k-point sampling..

- Submit your jobs ( 19 jobs ). Refer to `LESSON_04` for job submission.

- Once your jobs are done, remove all `tmp` folder.

### :large_blue_diamond: Collecting data

In each `output.dat` file, the following information need to be extracted:

- Number of grid point used in k-point sampling
- Number of k-points in calculations 
- Total energy
- Total running time (Walltime)

Save the information to a 4-column data file, named `DATA.dat`.

### :large_blue_diamond: Plot the data

Three plots are needed:

- Total energy vs Number of grid point used in k-point sampling 
- Total energy vs Number of k-points in calculations 
- Total running time vs Number of k-points in calculations

You need to write one or two short codes, using `python` and `matplotlib` (refer to `LESSON_02`). The code needs to read data from `DATA.dat`. Use appropriate `label` each axis (check out `matplotlib` manual for instruction)).

## :three: Final products
The following items must be delivered to `LESSON_06` folder before moving to `LESSON_07`

:heavy_check_mark: 19 folders, named `K-k`, where `k=4, 5, 6,..., 18`. In each folder, there are 3 files: `input.inp`, `output,dat`, and `job`.

:heavy_check_mark: `DATA.dat`, one or more python codes for plotting, and three plots for total energy and for total running time as function of  k-point smapling (see above for needed plots). **Use meaningful name for your files**

:heavy_check_mark: A file name `RESULTS.dat` with the following information:
- The smallest value of number of grid point used in k-point sampling that gives total energy converged with less than 1 meV accuracy.
- Explanation of why we need to find the mallest value of number of grid point used in k-point sampling (Hints: Computational time depends on k-point sampling).

**Important note** k-point sampling depends also on the `smearing` method and its parameters. These tutorials do not cover this aspect. Refer elsewhere for details.
