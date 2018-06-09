from dlePy.vasp.chgcar import reduce_chgcar

INDATA = 'CHGCAR'        # CHGCAR or PARCHG...
struc_file = 'CONTCAR'   # POSCAR, CONTCAR
factor = 3               # If = 1, no reduction
                         # Use this option for faster plotting

reduce_chgcar( INDATA, factor, CONTCAR=struc_file)
