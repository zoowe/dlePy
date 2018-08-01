from ase.build.molecule import molecule
from ase.io import write

benzene = molecule( 'C6H6' )

write( 'second_image.pov', benzene, format = 'pov', run_povray = True, 
        canvas_width = 150 # Set resolution of image's width to 300 pixels
     )
