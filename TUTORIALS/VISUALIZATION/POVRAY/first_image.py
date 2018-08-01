from ase.build.molecule import molecule
from ase.io import write

benzene = molecule( 'C6H6' )

write( 'first_image.pov', benzene, format = 'pov', run_povray = True )
