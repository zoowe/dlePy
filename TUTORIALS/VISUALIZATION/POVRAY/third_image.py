from ase.build.molecule import molecule
from ase.io import write
import numpy as np

benzene = molecule( 'C6H6' )

# Instead of using default value for radius of each atom, which is covalent radius, here we manually set radius for each atom.
# Default values come from ase.data import covalent_radii

r_C  = 0.4    
r_H  = 0.2

r = np.zeros( [ benzene.get_number_of_atoms() ] )
for i, at in enumerate( benzene ):
    if at.symbol == 'C':
        r[ i ] = r_C
    else:
        r[ i ] = r_H


write( 'third_image.pov', benzene, format = 'pov', run_povray = True, 
       canvas_width = 200,   # Set width, in pixel
       radii = r              # Set radius 
      )                 
