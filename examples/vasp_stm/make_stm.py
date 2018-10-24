from dlePy.vasp.stm import *

# Structure file (POSCAR or CONTCAR)
CONTCAR='CONTCAR'

# PARCHG file
PARCHG='SMALL.PARCHG'

# Prefix for OUTPUT
PREFIX='SOMETHING'

# Additional points to add to the top of cell. This is helpful
# if you have small vacuum.
TOP = 0

# Value of LDOS for scanning iso surface
LDOS=0.0001 

# Reference point, where you want to set to be zero, in crystal coordinate
REF = 0.0

# main command for getting stm image.
tersoff_hamann(CONTCAR,PARCHG,PREFIX,TOP,LDOS,REF)
