# PWSCF: Generating input file 

We can create input file manulally. However, the best practise is to use scripts to generate input files to:

- Make input file neat, clear
- Avoid mistakes
- Keep record of how and why all keywords and values are used

There are many ways to do this. `ASE` has its own interface. However, in this tutorial, we use `pwscf` included in `dlePy` package.

## :one: Purposes
Learning how to generate input file for `PWSCF` calculations by using `dlePy`

## :two: Hows
### :large_blue_diamond: Load `dlePy`

To use `pwscf` module of `dlePy`, you need to load it first by:

```
module load dlePy
```

### :large_blue_diamond: Study few examples

In this folder, there are three `write_input_XX.py` scripts:

- `write_input_default.py`: A very simple script to write default input file for an `scf` calculation of `bulk Ge`

The structure of the script is very simple, usually consisting of 3 components: `making structure`, `update or add keywords`, `write input file`.

Run this script, `python write_input_default.py`, a `PWSCF` input file named `input_default.inp` will be generated.

Inside `input_default.inp`, keywords are grouped to several blocks. Pay attention to the lines stated with `!`. They are corresponding to "name of blocks". They will be used to update value of keywords belong to these blocks or to add new keywords (see `write_input_01.py` and `write_input_02.py`

- `write_input_01.py`: Similar to `write_input_default.py` but with some example of how to change keywords' value. Run this, `input_01.inp` will be generated.

- `write_input_02.py`: Similar to `write_input_01.py` but with an example of adding a new keyword. Run this, `input_02.inp` will be generated.

## :three: Final products
The following items must be delivered to `LESSON_07` folder before moving to `LESSON_08`

:heavy_check_mark: A python script named `write_input.py` that can generate an input file for `pwscf` calculation named `input.inp` with the same keywords and values as shown in `LESSION_04/input.inp`, except for kinetic energy cutoff for wave function and for charge density and k-point mesh. Use value obtained in `LESSON_05` and `LESSON_06` accordingly.


