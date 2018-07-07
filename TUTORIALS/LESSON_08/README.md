# PWSCF: Optimize lattice parameter of Ge

## :one: Purposes
Use DFT to find/optimize lattice parameter of Ge and compare with experimental data

## :two: Hows
### :large_blue_diamond: Performing scf calculation for bulk Ge with a series of lattice parameters and collect data.
**Calculation for each lattice parameter must be in a folder named `A-X.XX`, where `X.XX` is the lattice parameter. In each folder, there must be a file named `make_input.py` which is used to create `input.inp` file (See [LESSON_07](../LESSON_07)). Use the optimized `kinetice energy cutoff` and `k-point mesh`. Copy `job` file from previous lesson to each folder.**

- Step 1: Perform `scf` calculation for buld Ge with lattice parameters ranging from 5.40 to 5.70 Angstrom with a step of 0.05 Angstrom.
- Step 2: Collect total energy for each value of lattice parameter. Save them (lattice parameter, total energy) to a 2 column data file named `DATA_01.dat`.
- Step 3: Find the value of lattice parameter ( refered as `a1` ) that gives the lowest total energy. 
- Step 4: Perform `scf` calculation for buld Ge with lattice parameters ranging from `a1 - 0.07` to `a1 + 0.07` Angstrom with a step of 0.01 Angstrom.
- Step 5: Collect total energy for each value of lattice parameter. Save them (lattice parameter, total energy) to a 2 column data file named `DATA_02.dat`.
- Step 6: Find the value of lattice parameter ( refered as `a2` ) that gives the lowest total energy. 
- Step 7: Collect 11 data points ( lattice parameter, total energy) for lattice parameters ranging from `a2 - 0.05` to `a2 + 0.05` Angstrom with a step of 0.01 Angstrom. If any data is not available, perform additional scf calcuation to obtain data. Save data to 2 column data file named `DATA.dat`
- Step 8: Plot data from `DATA.dat`. We want a figure like the one in https://wiki.fysik.dtu.dk/ase/tutorials/eos/eos.html (need line and scatter plot). Save the plot as `aVSe.png`

### :large_blue_diamond: Determine the optimized lattice parameter

In this step, we will fit the data obtained above to an equation of states to determine the optimized lattice parameters. Check out this link for reference to [equation of states](https://wiki.fysik.dtu.dk/ase/ase/eos.html)

We will have to write a short code to perform fitting. The code must do the following:
```
Read the data from DATA.dat to two array: a and e
```
Since we need volume for equation of states so:
```
Use simple math to calculate volume from lattice parameters.
v = function of a
```
Now, fit the `v` and `e` to equation of states. Check out this [link](https://wiki.fysik.dtu.dk/ase/tutorials/eos/eos.html) for example how it is done.
```
eos = EquationOfState(volumes, energies)
v0, e0, B = eos.fit()
```
`v0, e0, B` are the optimized volume, total energy, and bulk modulus of Ge.
Finaly, use simple math to calculate to optimized lattice parameter `a0`
```
a0 = function of v0
```

## :three: Final products
The following items must be delivered to `LESSON_08` folder before moving to `LESSON_09`

:heavy_check_mark: All folders, named `A-X.XX`, where `X.XX` is lattice parameter. In each folder, there are 3 files: `make_input.py`, `input.inp`, `output,dat`, and `job`.

:heavy_check_mark: `DATA_01.dat`, `DATA_02.dat`, `DATA.dat`, `aVSe.png`

:heavy_check_mark: A file name `RESULTS.dat` with the following information:
- Optimized lattice parameter of Ge
- Bulk modulus of Ge
- Comparison of the above value with experimental one.

