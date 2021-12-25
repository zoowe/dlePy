import os, json
import gzip  as gz
import numpy as np
import json
from sys import getsizeof
from ..small_tools import str_decode
from ..gaussian import gauss_smooth
from tqdm import tqdm

def parse_procar( procar, spin = False, soc = False ):
    print ( "Parsing %s" %( procar ) )
    try:
        if '.gz' in procar:
            decode = True
            with gz.open( procar, 'rb' ) as f:
                lines = f.readlines( )
        else:
            decode = False
            with open( procar, 'r' ) as f:
                lines = f.readlines( )
    except:
        print ( "Error with reading %s" %( procar ) )
        exit( ) 

    data = { }

    ( nkpts, nbands, nions ) = parse_info( lines[ 1 ] )
    ispin = 1
    if spin:
        ispin = 2
    data[ 'NKPTS'    ] = nkpts
    data[ 'NIONS'    ] = nions
    data[ 'NBANDS'   ] = nbands
    data[ 'ISPIN'    ] = ispin
    spins = { 0: 'SPINUP', 1: 'SPINDN' }
    i_start = 2
    for spin in range( ispin ):
        data[ spins[ spin ] ] = { }
        for ik in range( 1, nkpts + 1 ):
            for i in range( i_start, len( lines ) ):
                if 'k-point ' in str_decode( lines[ i ], decode ):
                    k_point = {}
                    k_point = parse_kline( str_decode( lines[ i ], decode), k_point )
                    if k_point[ 'IK' ] != ik:
                       raise( 'Index for k-point mismatched {} vs {}'.format( ik, k_point[ 'IK' ] ) )
                    k_point[ 'BANDS' ] = { }
                    for j in range( i, len( lines ) ):
                        if 'band' in str_decode( lines[ j ], decode ):
                            band = { }
                            band = parse_bline( str_decode( lines[ j ], decode ), band )
                            band[ 'PROJECTIONS' ] = { }
                            chars = [ str_decode( item, decode ) for item in lines[ j + 2 ].split( )[ 1: ] ]
                            band[ 'PROJECTIONS' ][ 'total' ] = {}
                            for k in range( j + 3, j + 3 + nions + 1 ):
                                tmp = lines[ k ].split( )
                                ion = str_decode( tmp[ 0 ], decode )    
                                band[ 'PROJECTIONS' ][ 'total' ][ ion ] = { }
                                for l in range( len( chars ) ):
                                    band[ 'PROJECTIONS' ][ 'total' ][ ion ][ chars[ l ] ] = tmp[ l + 1 ]

                            #For soc:
                            if soc:
                                ms = { 1 : 'mx', 2 : 'my', 3 : 'mz' }
                                for im in sorted( ms.keys( ) ):
                                    start_line = j + 3 + im * ( nions + 1 )
                                    end_line   = j + 3 + ( im + 1 ) * ( nions + 1  )
                                    band[ 'PROJECTIONS' ][ ms[ im ] ] = {}
                                    for k in range( start_line, end_line ):
                                        tmp = lines[ k ].split( )
                                        ion = str_decode( tmp[ 0 ], decode )
                                        band[ 'PROJECTIONS' ][ ms[ im ] ][ ion ] = { }
                                        for l in range( len( chars ) ):
                                            band[ 'PROJECTIONS' ][ ms[ im ] ][ ion ][ chars[ l ] ] = tmp[ l + 1 ]

                            k_point[ 'BANDS' ][ band[ 'IB' ] ] = band
                            if band[ 'IB' ] == nbands:
                                i_start = j  
                                break
                    data[ spins[ spin ] ][ ik ] = k_point
                    break
    del lines 
    return data 

def parse_info( line ):
    nkpts  = int( line.split( )[ 3  ] ) 
    nbands = int( line.split( )[ 7 ] ) 
    nions  = int( line.split( )[ -1  ] )
    return ( nkpts, nbands, nions )

def parse_kline( line, k_point ):
    tmp    = line.replace( '-', ' -').replace( 'k -','k-').split( )
    ik     = int( tmp[ 1 ] )
    coor   = ( float( tmp[ 3 ] ), float( tmp[ 4 ] ), float( tmp[ 5 ] ), )
    weight = float( tmp[ -1 ] )
    k_point[ 'IK' ]     = ik
    k_point[ 'COOR' ]   = coor
    k_point[ 'WEIGHT' ] = weight
    return k_point

def parse_bline( line, band ):
    tmp    = line.split( )
    ib     = int( tmp[ 1 ] )
    energy = float( tmp[ 4 ] )
    occ    = float( tmp[ -1 ] )
    band[ 'IB'     ] = ib
    band[ 'EIGENVAL' ] = energy
    band[ 'OCC'    ] = occ
    return band

def parse_procars( procars, spin = False, soc = False ):
    data = { }
    for i in sorted( procars.keys( ) ):
        data[ i ] = parse_procar( procars[ i ], spin = False, soc = False )
    return data

def get_bandstructure_( procar, spin = False, soc = False, atoms = [ ], characters = [ ], mag_components = [ 'total' ] ):
    '''
    mag_components = [ 'total', 'mx', 'my', 'mz' ]
    '''
    data = parse_procar( procar, spin = spin, soc = soc )
    nkpts  = data[ 'NKPTS' ]
    ispin  = data[ 'ISPIN' ]
    nbands = data[ 'NBANDS' ]
    spins = { 0: 'SPINUP', 1: 'SPINDN' }
    bandstructure_ = { }
    for i_spin in range( data[ 'ISPIN' ] ):
        kpoints  = [ ]
        eigenval = [ ]
        kweight  = [ ]
        pdos     = [ ]
        for ik in sorted( data[ spins[ i_spin ] ].keys( ) ):
            for ib in sorted( data[ spins[ i_spin ] ][ ik ][ 'BANDS' ].keys( ) ):
                kpoints.append( ik )
                eigenval.append( data[ spins[ i_spin ] ][ ik ][ 'BANDS' ][ ib ][ 'EIGENVAL' ] )
                kweight.append( data[ spins[ i_spin ] ][ ik ][ 'WEIGHT' ] )
                 
                if len( atoms ) == 0:
                    atoms = [ 'tot' ]
                if len( characters ) == 0:
                    characters = [ 'tot' ]
                tmp = 0
                for at in atoms:
                    for char in characters:
                        for mag in mag_components:
                            #print ( ik, ib, at, char, mag, data[ spins[ i_spin ] ][ ik ][ 'BANDS' ][ ib ][ 'PROJECTIONS' ][ mag ][ at ][ char ] )
                            print ( ik, ib, at, char, mag, i_spin )  
                            tmp += float( data[ spins[ i_spin ] ][ ik ][ 'BANDS' ][ ib ][ 'PROJECTIONS' ][ mag ][ at ][ char ] )
                #print ( tmp )
                pdos.append( tmp )
        bandstructure_[ spins[ i_spin ] ] = { 'NKPTS'    : nkpts, 
                                              'NBANDS'   : nbands, 
                                              'KPOINTS'  : kpoints, 
                                              'EIGENVAL' : eigenval, 
                                              'KWEIGHT'  : kweight,
                                              'PDOS'     : pdos }

    return bandstructure_

def get_bandstructure( procars, spin = False, soc = False, atoms = [ ], characters = [ ], mag_components = [ 'total' ] ):
    data = { }
    for i_procar in sorted( procars.keys( ) ):
        data[ i_procar ] = get_bandstructure_( procars[ i_procar ], spin, soc, atoms, characters, mag_components )
    bandstructure = { }
    for spin in data[ sorted( data.keys( ) )[ 0 ] ].keys( ):
        bandstructure[ spin ] = data[ sorted( data.keys( ) )[ 0 ] ][ spin ]
        for i_procar in sorted( data.keys( ) )[ 1: ]:
            starting_k = bandstructure[ spin ][ 'NKPTS' ]
            for key in bandstructure[ spin ].keys( ):
                if key == 'KPOINTS':
                    kpoints = ( np.array( data[ i_procar ][ spin ][ key ] ) + starting_k ).tolist( )
                    bandstructure[ spin ][ key ] += kpoints
                elif key == 'NBANDS':
                    pass
                else:
                    bandstructure[ spin ][ key ] += data[ i_procar ][ spin ][ key ]
    return bandstructure

def smooth_projection( nbands, k, e, pdos, emin, emax, de, sigma = None ):
    if sigma == None:
        sigma = 2 * de
    emin += -1
    emax +=  1
    erange = np.arange( emin, emax + de, de )
    ne  = len( erange )
    ne  = int( ( emax - emin ) / de ) + 1
    nk  = int( len( pdos ) / nbands )
    e   = e.reshape( nk, nbands )
    dos =  pdos.reshape( nk, nbands ) 
    dos_out = np.zeros( [ nk, ne ] )
    for i in tqdm( range( nk ) ):
        x0   =  e[ i, : ]
        amp  =  dos[ i, : ]
        e_, dos_ = gauss_smooth( x0, amp, sigma, 0, de, scale=1., xmin=emin, xmax=emax, verbose = False )
        dos_out[ i, : ] = dos_
    return dos_out, ( emin, emax, de ) 

def get_pdos( data, iions, spin = 'SPINUP' ):
    print ("Summing projected density of states for atoms: ", iions, "with " + spin )
    pdos = []
    characters  = data[ spin ][ 1 ][ 'BANDS'][ 1 ][ 'PROJECTIONS' ][ 'total' ][ '1' ].keys()
    for ik in [ i for i in sorted( data[ spin ].keys() ) if isinstance( i, int) ]:
        for ib in sorted( data[ spin ][ ik ][ 'BANDS'].keys( ) ):
            eigenval = data[ spin ][ ik ][ 'BANDS' ][ ib ][ 'EIGENVAL' ]
            weight   = data[ spin ][ ik ][ 'WEIGHT' ]
            proj     = [ ]
            for iion in iions:
                tmp = [ ]
                for char in characters:
                    tmp.append( float( data[ spin ][ ik ][ 'BANDS'][ ib ][ 'PROJECTIONS' ][ 'total' ][ str(iion) ][ char ] ) )
                proj.append( tmp )
                #proj.append( data[ spin ][ ik ][ 'BANDS'][ ib ][ 'PROJECTIONS' ][ 'total' ][ str(iion) ] )
            proj     = np.array( proj ).sum( axis = 0 )
            proj     = [ proj[ i ] for i in range( proj.shape[ 0 ] ) ]
            pdos.append( [ eigenval, weight ] + proj )
    pdos = np.array( pdos )
    character  = data[ spin ][ 1 ][ 'BANDS'][ 1 ][ 'PROJECTIONS' ][ 'total' ][ str(iion) ].keys()
    dataout = { }
    dataout[ 'EIGENVAL' ] = pdos[ :, 0 ]
    dataout[ 'WEIGHT' ]   = pdos[ :, 1 ]
    for i,char in enumerate( characters ):
        dataout[ char ] = pdos[ :, 2 + i ]

    print ( "Output is a dictionany with the following keys: " )
    print ( dataout.keys() )
    return dataout

def smooth_dos( dataout, fout, character = 'tot', sigma = 0.1, resolution = 0.01 ):
    d    = 20
    e    = dataout[ 'EIGENVAL' ]
    w    = dataout[ 'WEIGHT' ]
    pdos = dataout[ character ]
    ener, pdos = gauss_smooth( e, pdos * w, sigma, d, resolution)

    with open( fout, 'w' ) as f:
        for data in zip( ener, pdos ):
            f.write( "%12.8f %12.8f" %( data[ 0 ], data[ 1 ]) )
 
if __name__ == "__main__":
    import json
    data = parse_procar( 'D.0/PROCAR.gz' )
    '''
    with open( 'database.json', 'w' ) as f:
        json.dump( data, f )
    '''
    print ( data[ 'SPINDN' ].keys() )
