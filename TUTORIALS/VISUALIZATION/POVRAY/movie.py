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

# manually set color for each atom. Default colors are good, but sometime, we need to set our own colors.
# Colors can be defined using RGB format, and they are saved in a tuple or a list
colors = [ ]
for i, at in enumerate( benzene ):
    if at.symbol == 'C':
        # Set color for C: from black to gray then to white
        colors.append( (0. + i* 1/6., 0. + i * 1/6., 0. + i * 1/6. ) )
    else:
        # Set color for all H atoms: green
        colors.append( (0., 1., 0.) )


# For making movie, it is important to put the molecule in a box for easier control the image.
benzene.set_cell( [20, 20, 20] )
benzene.center()

#Boundary for image
center = (10, 10, 10 )
width  = 6
height = 6
xmin = center[ 0 ] - width / 2.
ymin = center[ 1 ] - height / 2.
xmax = center[ 0 ] + width / 2.
ymax = center[ 1 ] + height / 2.

number_of_frames = 18
angle = 360 / float( number_of_frames )
for i in range( number_of_frames ):
    benzene.rotate( angle, 'z', center = center )
    write( 'frame_' + str( i ) + '.pov', benzene, format = 'pov', run_povray = True, 
       canvas_width = 600,    # Set width, in pixel
       radii = r,              # Set radius 
       bondatoms = bond_pairs, # Display bonds
       colors = colors,        # Set colors
       bbox  = ( xmin, ymin, xmax, ymax ) # set boundary for image
      )           

# Make movie
#List images
s = ' '
for i in range( number_of_frames ):
    s +=  ' frame_' + str( i ) + '.png' 

#Delay     
delay = 50  # 50 1/100 s

#loop
loop = 0    # Number of loops. 0 = infinit. 

# Convert *png to gif
import os
os.system( 'convert -delay ' + str(delay ) + ' -loop ' + str( loop ) + s + ' movie.gif' )
