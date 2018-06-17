# Partial Phonopy

These scripts allow to use `phonopy` to perform vibrational calculation, using finit different method, for part of a system.

These scripts also allow to collect data from two different runs so one do not have to redo calculations that have been done.

# How this example works

## Code needed

I just code this one. So there may be bugs. Improments can be done.

```
module load dlePy/alpha
```

**Note:** Only alpha version works.


## make_POSCAR.py

Create a 5 layer (1x1) Cu slab with one H atom adsorbed.

**Important Note:** The code will work only if atoms are sorted as following:

```
All first layer atoms
All second layer atoms
...
All top layer atoms
All adsorbate atoms (can be devided to many layers)
```


## first_run.py

Create displacemnt for top 2 Cu layer and H atoms. Follow instruction printed after running this script for how and where to run vasp.

## second_run.py

Create displacemnt for top 4 Cu layer and H atoms. It also instructs how to collect calculations that have been done in first_run.py so we do not have to redo them. Follow instruction printed after running this script for how and where to run vasp.

## third_run.py

Same as second_run.py, but it creates displacemnt for all atoms and instructs to collect all data from second_run to save time. Follow instruction printed after running this script for how and where to run vasp.

## REFERENCES folder

If every thing is going as planned, you should see the same data/files/foldes in this REFERENCES folder.
