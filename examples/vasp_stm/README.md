# How to generate STM image from VASP run

## Step 1: Write a WAVECAR with lot of k-points.
NBANDS must be large. 

## Step 2: Generate PARCHG

- Make a folder 

- Link WAVECAR

- Copy POSCAR, INCAR, KPOINTS, POTCAR to this folder

- In INCAR: Need to make sure

```
NSW = 0; IBRION = -1
ICHARG = 0
ISTART = 1
LPARD = .TRUE.
NBMOD = -3
EINT  = 0 1   #This is for all unoccupied states from the Fermi energy to 1eV
EINT  = -1 0  #This is for all occupied states from the Fermi energy to -1eV
```
## Step 3: Visulize the PARCHG
Since the file is large, you should not scp it to your machine. Instead, reduce the grid by using this script:

```python
from dlePy.vasp.chgcar import reduce_chgcar

INDATA = 'PARCHG'        # CHGCAR or PARCHG...
struc_file = 'CONTCAR'   # POSCAR, CONTCAR
factor = 2               # If = 1, no reduction
                         # Use this option for faster plotting

reduce_chgcar( INDATA, factor, CONTCAR=struc_file)
```

Now, visualize the ```SMALL.PARCHG.vasp``` with iso-surface plot to find the best local density of states value so you can an iso-surface that cover the entire surface. Note this value.

## Step 4: Generate STM image
Use ```make_stm.py``` with approriate ```LDOS``` value and name of ```PARCHG```



