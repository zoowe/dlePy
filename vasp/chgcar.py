"""
****************************************
Author: Duy Le
Department of Physics
University of Central Florida

Email: duy.le@ucf.edu
Website: http://www.physics.ufc.edu/~dle
****************************************
"""

import numpy as np
import datetime
from ase.atoms import Atoms
from ase.units import Bohr
from ase.io import read, write

def read_chgcar(INDATA, CONTCAR='CONTCAR'):

    print ( 'Input data: ',INDATA )

    #Read CONTCAR to define the system.
    system=read(CONTCAR)

    #We don't need this, thus we need to count how many lines
    #we need to ignore
    nline     = 9
    startline = nline + system.get_number_of_atoms()

    #Open file 
    file=open( INDATA,'r')

    #Read all ignored lines
    for i in range( startline ):
        dump = file.readline( )

    #Read the dimensions of the array
    ngr = file.readline( ).split( )

    #Make it becomes integers 
    ng  = ( int( ngr[ 0 ] ), int( ngr[ 1 ] ), int( ngr[ 2 ] ) )

    print ( 'Read ',INDATA )
    print ( datetime.datetime.now( ) )

    total = np.empty( ng[ 0 ] * ng[ 1 ] * ng[ 2 ] )
    total = np.fromfile( file, count = ng[ 0 ] * ng[ 1 ] * ng[ 2 ], sep=' ')
    total = total.reshape( ng[ 2 ], ng[ 1 ], ng[ 0 ] ).T

    print ( datetime.datetime.now( ) )

    return total / system.get_volume( )

def write_chgcar( fobj, atoms, data = None):

    del atoms.constraints 
    write( fobj, atoms, format='vasp', direct = True, vasp5 = True )
    
    if isinstance( fobj, str ):
        fobj = open( fobj, 'a' )
        
    fobj.write('\n')
    
    volume = atoms.get_volume()
    if data is None:
        data = np.ones( ( 2, 2, 2 ) )
    
    #CHGCAR format uses DATA*volume
    data = np.asarray( data ) * volume

    if data.dtype == complex:
        data = np.abs( data )


    shape = np.array( data.shape )
    
    fobj.write(' %5i %5i %5i \n' % ( shape[ 0 ], shape[ 1 ], shape[ 2 ] ) )


    # Make a 1D copy of chg, must take transpose to get ordering right
    chgtmp = data.T.ravel( )

    # Must be a tuple to pass to string conversion
    chgtmp = tuple( chgtmp )

    # Other formats - 5 columns
    # Write all but the last row
    print ( datetime.datetime.now( ) )

    for ii in range( ( len( chgtmp ) - 1 ) / 5 ):
         fobj.write(' %17.10E %17.10E %17.10E %17.10E %17.10E\n' % chgtmp[ ii * 5 : ( ii + 1 ) * 5])

    # If the last row contains 5 values then write them without a newline
    if len( chgtmp ) % 5 == 0:
         fobj.write(' %17.10E %17.10E %17.10E %17.10E %17.10E' % chgtmp[ len( chgtmp ) - 5 : len( chgtmp ) ] )

    # Otherwise write fewer columns without a newline
    else:
         for ii in range( len( chgtmp ) % 5 ):
              fobj.write( (' %17.10E') % chgtmp[ len( chgtmp ) - len( chgtmp ) % 5 + ii ] )

    # Write a newline whatever format it is
    fobj.write( '\n' )

    print ( datetime.datetime.now() )
    # Clean up
    del chgtmp


def reduce_array( data, factor ):
    ng = data.shape
    if np.sum ( np.array( ng ) % factor ) != 0:
        raise RuntimeError ( """ERROR: Remainder of ng / factor must be zero. ng = { }""".format( ng ) )

    ng_    = np.array( ng ) / factor
    data_  = np.empty( ng_ )
    data_  = data[ ::factor, ::factor, ::factor ]
   
    return data_

def reduce_chgcar( INDATA, factor, CONTCAR = 'CONTCAR' ):
    OUTPRE = 'Reduced.' + INDATA + '.vasp'
    print ( 'Input data: ', INDATA )

    #Read CONTCAR to define the system.
    system = read( CONTCAR )

    total  = read_chgcar( INDATA, CONTCAR = CONTCAR )
    total_ = reduce_array( total, factor ) 

    print ( 'Write Charge' )
    write_chgcar( OUTPRE, system, data = total_ )

def reduce_spin_chgcar(INDATA, factor,CONTCAR='CONTCAR'):

    OUTPRE=INDATA+'.vasp'
    if factor > 1:
        OUTPRE='Reduced.'+INDATA+'.vasp'

    print ( 'Input data: ',INDATA )
    
    #Read CONTCAR to define the system.
    system = read( CONTCAR )
    
    
    #We don't need this, thus we need to count how many lines
    #we need to ignore
    nline = 9
    startline = nline + system.get_number_of_atoms( )
    
    #Open file 
    file=open( INDATA, 'r' )
    
    #Read all ignored lines
    for i in range( startline ):
        dump = file.readline( )
    
    #Read the dimensions of the array
    ngr = file.readline( ).split( )
    #Make it becomes integers 
    ng = ( int( ngr[ 0 ] ), int( ngr[ 1 ] ), int( ngr[ 2 ] ) )

    print ( 'Read total charge' )

    total = np.empty( ng[ 0 ] * ng[ 1 ] * ng[ 2 ] )
    total = np.fromfile( file, count = ng[ 0 ] * ng[ 1 ] * ng[ 2 ], sep=' ')
    total = total.reshape( ng[ 2 ], ng[ 1 ], ng[ 0 ] ).T

    total_ = reduce_array( total, factor ) 
    del (total)

    print ( 'Write Charge' )
    write_chgcar( OUTPRE, system, data = total_ )
    
    for i in range(10000):
        dump = file.readline( )
        ngr_ = dump.split( )  
        if ngr[ 0 ] == ngr_[ 0 ]:
           break
    print ( 'Read SPIN' )

    spin = np.empty(ng[0]*ng[1]*ng[2])
    spin = np.fromfile(file, count = ng[0]*ng[1]*ng[2], sep=' ')
    spin = spin.reshape (ng[2], ng[1],ng[0]).T

    spin_ = reduce_array( spin, factor ) 
    del (spin)
    
    file.close()

    print ( 'Write Spin' )
    write_chgcar('SPIN' + OUTPRE, system, data = spin_ )

    print ( 'Now, working on SPIN_UP' )
    spin_up = np.empty( spin_.shape )

    print ( 'Now, calculating  SPIN_UP' )
    spin_up = ( total_ + spin_ ) / 2.

    print ( 'Now, write  SPIN_UP' )
    write_chgcar('SPIN_UP.' + OUTPRE, system, data = spin_up )
    del ( spin_up )

    print ( 'Now, working on SPIN_DN' )
    spin_dn=np.empty( spin_.shape )
    print ( 'Now, calculating  SPIN_DN' )
    spin_dn = ( total_ - spin_ ) / 2.

    print ( 'Now, writing  SPIN_DN' )
    write_chgcar('SPIN_DN.' + OUTPRE, system, data = spin_dn )
    del ( spin_dn )
    del ( total_ )
    del ( spin_ )


def write_vtk( fobj, atoms, data = None):
    del atoms.constraints 
    
    if isinstance( fobj, str ):
        fobj = open( fobj, 'w' )
        
    shape = np.array( data.shape )    
    
    fobj.write('# vtk DataFile Version 3.0 \n')
    fobj.write('vtk output \n')
    fobj.write('ASCII \n')
    fobj.write('DATASET STRUCTURED_GRID \n')
    fobj.write('DIMENSIONS  %5i %5i %5i \n' % ( shape[0], shape[1], shape[2]))
    #fobj.write('SPACING 1 1 1 \n')
    #fobj.write('ORIGIN 0 0 0 \n')
    fobj.write('POINTS %10i float \n' %( shape[0]*shape[1]*shape[2]))
    """ Write coordinate of points """
    ix,iy,iz = np.mgrid [0:shape[0], 0:shape[1], 0:shape[2]]
    
    x = ix / float (shape[0]) * atoms.cell[0,0] + \
        iy / float (shape[1]) * atoms.cell[1,0] + \
        iz / float (shape[2]) * atoms.cell[2,0] 
    y = ix / float (shape[0]) * atoms.cell[0,1] + \
        iy / float (shape[1]) * atoms.cell[1,1] + \
        iz / float (shape[2]) * atoms.cell[2,1] 
    z = ix / float (shape[0]) * atoms.cell[0,2] + \
        iy / float (shape[1]) * atoms.cell[1,2] + \
        iz / float (shape[2]) * atoms.cell[2,2] 
    x = x.ravel()    
    y = y.ravel()    
    z = z.ravel()    
    #print  len (x), len(y), len(z), shape[0]*shape[1]*shape[2]
    for i in range ( len ( x ) ):
         fobji.write( "{} {} {}\n".fotmat( x[i],y[i],z[i] ) )
    fobj.write('POINT_DATA %10i \n' %( shape[0]*shape[1]*shape[2]))     
    fobj.write('SCALARS Density double \n')
    fobj.write('LOOKUP_TABLE default \n')
 
    volume = atoms.get_volume()
    if data is None:
        data = np.ones( ( 2, 2, 2 ) )
    
    #CHGCAR format uses DATA*volume
    data = np.asarray( data ) 

    if data.dtype == complex:
        data = np.abs( data )


    # Make a 1D copy of chg, must take transpose to get ordering right
    chgtmp = data.T.ravel()
    # Must be a tuple to pass to string conversion
    chgtmp = tuple(chgtmp)
    # Other formats - 5 columns
    # Write all but the last row
    for ii in range( ( len( chgtmp ) - 1 ) / 5 ):
         fobj.write(' %17.10E %17.10E %17.10E %17.10E %17.10E\n' % chgtmp[ ii*5: (ii+1) * 5])
    # If the last row contains 5 values then write them without a newline
    if len( chgtmp ) % 5 == 0:
         fobj.write(' %17.10E %17.10E %17.10E %17.10E %17.10E' % chgtmp[ len( chgtmp ) - 5 : len( chgtmp ) ] )
    # Otherwise write fewer columns without a newline
    else:
         for ii in range( len( chgtmp ) % 5 ):
              fobj.write((' %17.10E') % chgtmp[ len( chgtmp ) - len( chgtmp ) % 5 + ii ] )
    # Write a newline whatever format it is
    fobj.write('\n')
    # Clean up
    del chgtmp


def average(data):
#    system=read(CONTCAR)
    LCPOT=np.zeros([data.shape[2]])
    f=open('AVELCPOT','w')
    for k in range(data.shape[2]):
        LCPOT[k]=np.sum(data[:,:,k])/float(data.shape[0]*data.shape[1])
        f.write( "{} {}\n".format( k,LCPOT[k] ) )
 
