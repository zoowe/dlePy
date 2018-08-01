from ase.io import write
from ase.cluster.icosahedron import Icosahedron
import matplotlib.cm as cmx
import os

system=Icosahedron('Ag',3)
for j in range(3):
    d=max(system.positions[:,j])-min(system.positions[:,j])
    system.cell[j,j]=d+15.

system.center()

center = (system.cell[0,0]/2., system.cell[1,1]/2., system.cell[2,2]/2.)
width = max(system.positions[:,0])-min(system.positions[:,0])
height = max(system.positions[:,1])-min(system.positions[:,1])
xmin = min(system.positions[:,0])
xmax = max(system.positions[:,0])

colors = [ ]
map = cmx.rainbow
for i, at in enumerate( system ):
    code = int( ( at.position[ 0 ] - xmin ) / width * 255 )
    color_ = map( code )
    color = (color_[0], color_[1], color_[2], 0.5 )
    colors.append( color )

pad = 3
xmin = center[ 0 ] - width / 2. - pad
ymin = center[ 1 ] - height / 2.- pad
xmax = center[ 0 ] + width / 2. + pad
ymax = center[ 1 ] + height / 2.+ pad

number_of_frames = 180
angle = 360 / float( number_of_frames )
for i in range( number_of_frames ):
    print i
    system.rotate( -angle, 'y', center = center )
    write( 'ico_' + str( i ) + '.pov', system, format = 'pov', run_povray = True,
       canvas_width = 100,    # Set width, in pixel
       colors = colors,        # Set colors
       bbox  = ( xmin, ymin, xmax, ymax ) # set boundary for image
      )


# Make movie
#List images
s = ' '
for i in range( number_of_frames ):
    s +=  ' ico_' + str( i ) + '.png' 

#Delay     
delay = 5  # 5 1/100 s

#loop
loop = 0    # Number of loops. 0 = infinit. 

# Convert *png to gif
os.system( 'convert -delay ' + str(delay ) + ' -loop ' + str( loop ) + s + ' ico.gif' )
