# PWSCF: Electronic band structure of Ge 

## :one: Purposes
- Understand the Brillouin zone of bulk Ge
- Calculate electronic band structure of Ge

## :two: Hows
### :large_blue_diamond: STEP 1: SCF calculation

Run a `scf` calculation for bulk Ge in a folder named `SCF`.

- We need to write a script to generate input file that uses lattice constant, energy cutoff, k-point mesh obtained from previous lessons. In addition, `disk_io = 'low'` (or something other than `'none'`) and `tstress = .true.`. 
- Check the result to see if the total stress is close to zero (if not, we may have to reoptimize lattice parameter again). 
- **This step is for getting charge density**, which will be written in `outdir`. 
- **NOTE:** `prefix` must be the same for all calculations performed in this `LESSON`.

### :large_blue_diamond: STEP 2: BANDS (NSCF) calculation

Run a `bands` calculation for bulk Ge in a folder named `BANDS_01`.

- This step is to calculate eigenvalues for a series of idividual k-points along symmetry paths of the Brillouin Zone. Google for Brillouin Zone of `Ge` or `fcc` structure for finding the paths (we can also use `xcrysden` software to visualize the Brillouin Zone).   

- For simplification, we need to calculate band structure along ![\text{L} \rightarrow  \Gamma  \rightarrow \text{X}](http://www.sciweavers.org/tex2img.php?eq=%5Ctext%7BL%7D%20%5Crightarrow%20%20%5CGamma%20%20%5Crightarrow%20%5Ctext%7BX%7D&bc=White&fc=Black&im=gif&fs=12&ff=arev&edit=0) path. We need to find coordinate of ![\text{L, }\Gamma\text{, and X}](http://www.sciweavers.org/tex2img.php?eq=%5Ctext%7BL%2C%20%7D%5CGamma%5Ctext%7B%2C%20and%20X%7D&bc=White&fc=Black&im=gif&fs=12&ff=arev&edit=0) points (in crystal coordinate) and put in to `gen_k.py` script to generate 31 k-points along the path.

- We need an input file, which is the same as the one in `scf` calculation in `Step 1`. However, `calculation` needs to be `bands` and `K_POINTS` section becomes:
```
K_POINTS crystal
31
 0.500000000000  0.500000000000  0.500000000000     1.000000
 0.490909090909  0.490909090909  0.490909090909     1.000000
 ...

```
The coordinates of 31 k-points were generated from `gen_k.py` above.

- **We also need `outdir` where charge density calculated in `Step 1` is save. We can defind `outdir` so that it points to the `outdir` of `scf` calclation. However, in this `LESSON`, we copy the whole `outdir` of `scf` run to this folder `BANDS_01`.**

- Modify the job file if needed to run `pw.x` for this input.

- Make sure the run is succesful.

### :large_blue_diamond: STEP 3: Get data for plotting

Run a `bands.x` to extract eigenvalue for each k-point.

- Create input file for `bands.x` run (name the file as `bands.inp`) with the following (or with your chosen `prefix`, `outdir`, and `filband`:
```
&BANDS
    prefix  = 'pwscf'
    outdir='./tmp'
    filband = 'bands.dat'
 /
```
Pleas refer to https://www.quantum-espresso.org/Doc/INPUT_BANDS.html for more infor about `bands.x`.

- Modify `job` file so it can run `bands.x` for `bands.inp` input. Name the output file as `bands.out`.

- Once this step is done, we should have `bands.dat` and `bands.dat.gnu` with all desriable data.

### :large_blue_diamond: STEP 4: Plot the band structure

In this step, we will plot band structure using data in `bands.dat.gnu`. We will need to write a python code to do the job, including the following steps:

- Read data from `bands.dat.gnu` to `kpt` and `eigenval`.

- Use scatter plot to plot the data. Horizontal axis shows k-points along the ![\text{L} \rightarrow  \Gamma  \rightarrow \text{X}](http://www.sciweavers.org/tex2img.php?eq=%5Ctext%7BL%7D%20%5Crightarrow%20%20%5CGamma%20%20%5Crightarrow%20%5Ctext%7BX%7D&bc=White&fc=Black&im=gif&fs=12&ff=arev&edit=0) path. Vertical axis shows eigen values for each k-points.

- Draw three vertical lines at the position of ![\text{L, }\Gamma\text{, and X}](http://www.sciweavers.org/tex2img.php?eq=%5Ctext%7BL%2C%20%7D%5CGamma%5Ctext%7B%2C%20and%20X%7D&bc=White&fc=Black&im=gif&fs=12&ff=arev&edit=0) points.

- Draw a horizontal line at to show the Fermi level (the value of it can be found in `output.dat`).

### :large_blue_diamond: STEP 5: More k-points for band structure

Repeat Step 2, 3, and 4 for 301 k-points along the same path in a folder named `BANDS_02`.


## :three: Final products
The following items must be delivered to `LESSON_09` folder before moving to `LESSON_10`

:heavy_check_mark: Three folders: `SCF`, `BANDS_01`, and `BANDS_02` with all files used and preduced in this LESSON.

:heavy_check_mark: A `RESULTS.dat` file with the following information:

- Comparison between the two plots (with 31 and 301 k-points). What do we learn from it?

- Comparison between the calculated band structure with experimental one and with other caculated band structures (from literature)


