# PWSCF: Kinetic Energy Cut-off

## :one: Purposes
Performing convergence test for Kinetic Energy Cut-off for wave-functions.

In most cases, the larger Kinetic Energy Cut-off is the more accurate, reliable calculation is. However, large Kinetic Energy Cut-off requires more computational time. Thus, we must determine the lowest possible Kinetic Energy Cut-off but give accurate results.

## :two: Hows
### :large_blue_diamond: Performing scf calculation for a series of Kinetic Energy Cut-off.

- In this step, a series of scf calculations will be performed with Kinetic Energy Cut-off for wavefunctions ranging from 20 Ry to 70 Ry with a step of 5 Ry (12 calucations).

- Do make sure that number of bands is large enough to include some empty states. 

- Each calculation must be done inside a folder named `EN-XX`, where `XX` is the `value` of Kinetic Energy Cut-off for wavefunctions.

- You must copy `input.inp` and `job` files from `LESSON_04` to each `EN-XX` folder, then edit `input.inp` file with new value of Kinetic Energy Cut-off for wavefunctions. **You must also use new value Kinetic Energy Cut-off for charge density that is 4 times of value of Kinetic Energy Cut-off for wavefunctions**  

- Submit your jobs ( 12 jobs ). Refer to `LESSON_04` for job submission.

- Once your jobs are done, remove all `tmp` folder.

### :large_blue_diamond: Collecting data

In each `output.dat` file, the following information need to be extracted:

- Kinetic energy cut-off
- Total energy
- Total running time (Walltime)

Save the information to a 3-column data file, named `DATA.dat`.

### :large_blue_diamond: Plot the data

Two plots are needed:

- Total energy vs Kinetic energy cut-off
- Total running time vs Kinetic energy cut-off

You need to write one or two short codes, using `python` and `matplotlib` (refer to `LESSON_02`). The code needs to read data from `DATA.dat`. Use appropriate `label` each axis (check out `matplotlib` manual for instruction).

## :three: Final products
The following items must be delivered to `LESSON_05` folder before moving to `LESSON_06`

:heavy_check_mark: 12 folders, named `EN-XX`, where `XX=20, 25, 30,..., 70`. In each folder, there are 3 files: `input.inp`, `output,dat`, and `job`.

:heavy_check_mark: `DATA.dat`, one or two python codes for plotting, and two plots for total energy and for total running time as function of  Kinetic energy cut-off. **Use meaningful name for your files**

:heavy_check_mark: A file name `RESULTS.dat` with the following information:
- The smallest value of Kinetic energy cut-off that gives total energy converged with less than 1 meV accuracy.
- Explanation of why finding optimal Kinetic energy cut-off is needed.

