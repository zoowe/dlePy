from ase.build.molecule import molecule
from ase.io import write

benzene = molecule( 'C6H6' )

write( 'second_image.pov', benzene, format = 'pov', run_povray = True, 
        canvas_width = 1000 # Set resolution of image's width to 1000 pixels
     )
