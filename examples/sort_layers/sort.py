from dlePy.strucmod import sort_layers
from ase.io import read, write

system = read( 'CONTCAR' )

L1 = [ at.index for at in system if at.position[ 2 ] > 21 ]
L2 = [ at.index for at in system if ( at.position[ 2 ] > 18 ) and ( at.position[ 2 ] < 21 )  ]
L3 = [ at.index for at in system if ( at.position[ 2 ] > 15 ) and ( at.position[ 2 ] < 18 )  ]
L4 = [ at.index for at in system if ( at.position[ 2 ] > 12 ) and ( at.position[ 2 ] < 15 )  ]
L5 = [ at.index for at in system if at.position[ 2 ] < 12  ]

layer = [ L1 + L5, L2 + L4, L3 ]

new_system = sort_layers( system, layer )
# This will put L1 + L5 to the bottom of POSCAR and L3 to the top of POSCAR

write( 'POSCAR', new_system, format = 'vasp', direct = True, vasp5 = True )

