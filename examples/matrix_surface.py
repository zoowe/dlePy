from ase.build       import molecule, fcc111, add_adsorbate
from dlePy.supercell import create_matrix_surface
from ase.visualize   import view

primitive = fcc111( 'Cu', size = ( 1, 1, 5 ), vacuum = 7.5 )
surface   = create_matrix_surface( primitive, matrix = ( 2, 1, -1, 4 ) )
mol       = molecule( 'CO' )
add_adsorbate( surface, mol, 1.5, position = ( 0., 0. ), mol_index = 1 )

view( surface )


