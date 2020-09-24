from ase import *
from ase.io import *
import numpy as np
from ..math import length

def rec2cart( recs, b ):
    '''
    recs: 3-elements array like. Reciprocal coordinate
    b: reciprocal vectors ( atoms.get_reciprocal_cell())
    output: carts in 2 pi / a
    '''

    carts = np.zeros( [ 3 ] )
    for i in range( 3 ):
        carts = recs[ 0 ] * b[ :, 0 ] + recs[ 1 ] * b[ :, 1 ] + recs[ 2 ] * b[ :, 2 ]

    return carts 
    
def vasp_band_gen_k(klist, npts0, ppf,b = None, hse=False):    

    klist_ = klist
    if b is None:
        coor_type = 'Reciprocal lattice'
    else:
        coor_type = 'Cartesian'
        klist = []
        for k in klist_:
            k_ = rec2cart( k, b )
            klist.append( k_ )
    
    d   = np.zeros( [ len( klist ) - 1    ] ) 
    x   = np.zeros( [ len( klist ) - 1, 3 ] )
    npd = np.zeros( [ len( klist ) - 1 ], dtype = int )

    for i in range(len(d)):
        x[ i ] = np.array( klist[ i + 1 ] ) - np.array( klist[ i ] )
        d[ i ] = length( x[ i ] )
    
    total_d=np.sum( d )

    for i in range( len( d ) ):
        npd[ i ] = int( d[ i ] / total_d * npts0 )
        print( 'Path ',i,' has ',npd[ i ] )
    
    npts = int( np.sum( npd ) + 1 )
    k  = np.zeros( [ npts, 3 ] )

    nstart=0
    for i in range( len( d ) ):
        a = np.array( klist[ i ] )
        b = np.array( klist[ i + 1] )
        k[ nstart : nstart + npd[ i ] ] = a + ( b - a ) * np.array( (range( npd[ i ] ), range( npd[ i ] ), range( npd[ i ] ) ) ).T / float( npd[ i ] )
        nstart += npd[ i ] 
    
    k[ npts - 1 ] = klist[ len( d ) ]
    
    nfile = int( npts / ppf )
    
    addition = False
    if nfile * ppf < npts:
        last_ppf = npts - nfile * ppf
        addition = True

    if  not hse:
        weight = 1.0
    else:
        weight = 0.0

    klist_file = open( 'KLIST', 'w' )
    for i in range( nfile ):
        file = 'KPOINTS.' + str( i )
        f = open( file, 'w' )
        if not hse:
            f.write( 'Automatically generated mesh\n' )
            f.write( str( ppf ) + '\n' )
            f.write( coor_type + '\n' )
        else:
            try:
                IBZKPT = open( 'IBZKPT', 'r' )
            except IOError:
                print ('ERROR: IBZKPT file is needed for hse = True')
                print ('SOLUTION: Copy IBZKPT file from SCF calculation to this folder')
                raise

            lines  = IBZKPT.readlines( )
            for i_line in range( len( lines ) ):
                if i_line == 1:
                    org_ppf = int( lines[ i_line ].strip( ) )
                    f.write( str( ppf + org_ppf ) + '\n' )
                else:
                    f.write( lines[ i_line ].strip( ) + '\n' )
            IBZKPT.close( )

        for ii in range( i * ppf, ( i + 1 ) * ppf ):
            f.write( "%15.12f %15.12f %15.12f %12.6f\n" %(k[ ii, 0 ],k[ ii, 1 ], k[ ii, 2 ], weight ) )
            klist_file.write( "%5i %15.12f %15.12f %15.12f \n" %( ii, k[ ii, 0 ],k[ ii, 1 ], k[ ii, 2 ] ) )
        f.close()
    
    if addition:
        file = 'KPOINTS.' + str( nfile )
        f=open( file, 'w' )
        if not hse:
            f.write( 'Automatically generated mesh\n' )
            f.write( str( last_ppf ) + '\n' )
            f.write( 'Reciprocal lattice\n' )
        else:
            IBZKPT = open( 'IBZKPT', 'r' ) 
            lines  = IBZKPT.readlines()
            for i_line in range( len( lines ) ):
                if i_line == 1:
                    org_ppf = int( lines[ i_line ].strip( ) )
                    f.write( str( last_ppf + org_ppf ) + '\n' )
                else:
                    f.write( lines[ i_line ].strip( ) + '\n' )
            IBZKPT.close()

        for ii in range( nfile * ppf, nfile * ppf + last_ppf):
            f.write( "%15.12f %15.12f %15.12f %15.6f\n" %(k[ ii, 0 ], k[ ii, 1 ], k[ ii, 2 ],weight) )
            klist_file.write( "%5i %15.12f %15.12f %15.12f \n" %(ii, k[ ii, 0 ], k[ ii, 1 ], k[ ii, 2 ] ) )
        f.close()
        klist_file.close()
