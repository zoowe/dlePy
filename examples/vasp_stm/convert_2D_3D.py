from ase import *
from ase.io import *
from dlePy.vasp.stm import convert_2D_to_3D

system=read('CONTCAR')
DATA2D='STM.-0.50.10'
DATA3D='SOMETHING.vasp'

convert_2D_to_3D(DATA2D,system,DATA3D)
