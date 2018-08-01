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


# Determine bonds between atoms. get_bondpairs will search all pairs with distance smaller than radius * sum of covalent bonds of the pair.
from ase.io.pov import get_bondpairs

bond_pairs = get_bondpairs( benzene, radius = 1.1 )

write( 'fourth_image.pov', benzene, format = 'pov', run_povray = True,
       canvas_width = 1000,    # Set width, in pixel
       radii = r,              # Set radius 
       bondatoms = bond_pairs  # Display bonds 
      )
