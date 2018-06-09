from ase import *
from ase.io import *
import numpy as np
from dlePy.math import length

def vasp_band_gen_k(klist, npts0, ppf,hse=False):    
    
    d   = np.zeros( [ len( klist ) - 1    ] ) 
    x   = np.zeros( [ len( klist ) - 1, 3 ] )
    npd = np.zeros( [ len( klist ) - 1 ], dtype = int )

    for i in range(len(d)):
        x[ i ] = np.array( klist[ i + 1 ] ) - np.array( klist[ i ] )
        d[ i ] = length( x[ i ] )
    
    total_d=np.sum( d )

    for i in range( len( d ) ):
        npd[ i ] = int( d[ i ] / total_d * npts0 )
        print 'Path ',i,' has ',npd[ i ]
    
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
            print >>f, 'Automatically generated mesh'
            print >>f, ppf
            print >>f, 'Reciprocal lattice'
        else:
            try:
                IBZKPT = open( 'IBZKPT', 'r' )
            except IOError:
                print 'ERROR: IBZKPT file is needed for hse = True'
                print 'SOLUTION: Copy IBZKPT file from SCF calculation to this folder'
                raise

            lines  = IBZKPT.readlines( )
            for i_line in range( len( lines ) ):
                if i_line == 1:
                    org_ppf = int( lines[ i_line ].strip( ) )
                    print >> f, ppf + org_ppf
                else:
                    print >> f, lines[ i_line ].strip( )
            IBZKPT.close( )

        for ii in range( i * ppf, ( i + 1 ) * ppf ):
            print >>f,"%15.12f %15.12f %15.12f %12.6f" %(k[ ii, 0 ],k[ ii, 1 ], k[ ii, 2 ], weight )
            print >>klist_file,"%5i %15.12f %15.12f %15.12f " %( ii, k[ ii, 0 ],k[ ii, 1 ], k[ ii, 2 ] )
        f.close()
    
    if addition:
        file = 'KPOINTS.' + str( nfile )
        f=open( file, 'w' )
        if not hse:
            print >>f, 'Automatically generated mesh'
            print >>f, last_ppf
            print >>f, 'Reciprocal lattice'
        else:
            IBZKPT = open( 'IBZKPT', 'r' ) 
            lines  = IBZKPT.readlines()
            for i_line in range( len( lines ) ):
                if i_line == 1:
                    org_ppf = int( lines[ i_line ].strip( ) )
                    print >>f, last_ppf + org_ppf
                else:
                    print >>f, lines[ i_line ].strip( )
            IBZKPT.close()

        for ii in range( nfile * ppf, nfile * ppf + last_ppf):
            print >>f,"%15.12f %15.12f %15.12f %15.6f" %(k[ ii, 0 ], k[ ii, 1 ], k[ ii, 2 ],weight)
            print >>klist_file,"%5i %15.12f %15.12f %15.12f " %(ii, k[ ii, 0 ], k[ ii, 1 ], k[ ii, 2 ] )
        f.close()
        klist_file.close()
