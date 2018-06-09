"""
Creating (7 0, 6 6 ) Cu(111) surface

"""

from ase.build import fcc111
from ase.io import write
from dlePy.supercell import create_matrix_surface

primitive_cell = fcc111( 'Au', ( 1, 1, 5 ), vacuum=7.5)
matrix = ( 7, 0, 6, 6 )
system = create_matrix_surface( primitive_cell, matrix  )

write( 'POSCAR', system, format = 'vasp', direct = True, vasp5= True)

