import numpy as np
import gzip as gz
import os

def read_file( filename ):
    if not os.path.isfile( filename ):   
        print ( "%s is not found" %( filename ) )
        exit

    if '.gz' in filename:
        with gz.open( filename, 'r' ) as f:
            lines = f.readlines( )
    else:
        with open( filename, 'r' ) as f:
            lines = f.readlines( )

    return lines

def add_key( data, name, line ):
    data[ name ] = int( line.split( )[ -1 ] )
    return data

def read_keys( data, lines, args ):
    for key in args.keys():
        print ( "Supplied {}:{}".format( key, args[ key ] ) )
        data = add_key( data, key, str( args[ key ] ) )
    list_key = [ 'natomwfc', 'nbnd', 'nkstot',  'npwx', 'nkb' ]
    for i in range( len( lines ) ):
       count_key = 0
       for name in list_key:
           if name in lines[ i ]:
               print ( "Add or update from output {}".format( name ) )
               data = add_key( data, name, lines[ i ] )
               count_key += 1
               if count_key == len( list_key ):
                   break
    nspin = 1
    for line in lines:
        if "spin up" in line:
            nspin  = 2
            data[ 'nkstot' ] = data[ 'nkstot' ] / 2
            break
    data[ 'nspin' ] = nspin
    
    return data

def read_atomic_states( data, lines ):
    # Read  Atomic states used for projection
    data[ "Atomic states" ] = {}
    count_state = 0
    for i in range( len( lines ) ):
        if "state #" in lines[ i ]:
            id = int( lines[ i ].split( ":" )[ 0 ].split( )[ - 1 ] ) - 1
            wfcname = lines[ i ].split( ":" )[ 1 ].replace( '\n', '' )
            data[ "Atomic states" ][ id ] = wfcname
            count_state += 1
            if count_state == data[ 'natomwfc' ]:
                break
    return data

def read_projection( data, lines ):
    ik = 0
    kdata = {}
    for i in range( len( lines ) ):
        if 'k =' in lines[ i ]:
            kdata[ ik ] = {}
            k_coor = [ float( x ) for x in lines[ i ].split()[-3:] ]
            kdata[ ik ][ 'coordinate' ] = k_coor
            ie = 0
            kdata[ ik ][ 'eigen' ]  = {}
            for j in range( i + 1, i + 1000 * data[ 'nbnd' ] ):
                if 'e' in lines[ j ]:
                    band = {}
                    e = float( lines[ j ].replace( '=','').split()[ -2 ] )
                    band[ 'e' ] = e
                    eline = lines[ j + 1 ] 
                    for k in range( j + 2, j + 20 ):
                        if '|psi|^2' in lines[ k ]:
                            break
                        else:
                            eline += lines[ k ].strip()
                    eline = eline.replace('\n','').replace( 'psi =', ' ' ) 
                    eline = eline.replace('*[#',' ').replace( ']+', ' ' )
                    eline = eline.split()
                    val = [ 0 ] * data[ 'natomwfc' ]
                    for ival in range( 0, len( eline ), 2 ):
                        val[ int(eline[ ival + 1 ] ) - 1] = float( eline[ ival ] )
                     
                    band[ 'psi' ] = val
                    kdata[ ik ][ 'eigen' ][ ie ] = band 
                    ie += 1
                    if ie == data[ 'nbnd']:
                        break
            ik += 1
            if ik == data[ 'nkstot' ] * data[ 'nspin' ]:
                break 
    data[ 'atomic projection' ] = kdata

    return data    

def read_lowdin( data, lines ):
    lowdin = { }
    for i, line in enumerate( lines ):
        if 'Lowdin Charges:' in line:
            start_line = i
            break
    for i in range( start_line, len( lines ) ):
        line = lines[ i ]
        if 'Atom #' in line:
            atom = {}
            tmp = line.replace( "#", '=').replace( ':',' ').replace( 'total charge', 'total' ).replace( ',',' ').replace('=', ' ').split()
            at  = int( tmp[ 1 ] ) - 1
            atom[ 'charge' ] = { }
            atom[ 'spin up' ] = { }
            atom[ 'spin down' ] = { }
            atom[ 'polarization' ] = { }
            for k in range( 1, len( tmp ) // 2 ):
                atom[ 'charge' ][ tmp[ 2*k ] ] =float ( tmp[ 2*k + 1 ] )
            for j in range( 100 ):
                line = lines[ i + j ]
                if 'polarization' in line:
                    tmp = line.replace( ',',' ').replace( '=',' ').replace( 'polarization', 'total' ).split()
                    for k in range( len( tmp ) // 2 ):
                        atom[ 'polarization' ][ tmp[ 2*k ] ] = float( tmp[ 2*k + 1 ] )
                    break
                if 'spin up' in line:
                    tmp = line.replace( 'spin up', 'total' ).replace( ',',' ').replace('=', ' ').split()
                    for k in range( len( tmp ) // 2 ):
                        atom[ 'spin up' ][ tmp[ 2*k ] ] = float( tmp[ 2*k + 1 ] ) 
                if 'spin down' in line:
                    tmp = line.replace( 'spin down', 'total' ).replace( ',',' ').replace('=', ' ').split()
                    for k in range( len( tmp ) // 2 ):
                        atom[ 'spin down' ][ tmp[ 2*k ] ] = float( tmp[ 2*k + 1 ] )
            lowdin[ at ] = atom
    data[ 'lowdin' ] = lowdin
    return data

def list_projections( data ):
    print ( "%-12s %-s" %( 'State Index', 'Projection' ) )
    for key in sorted( data[ "Atomic states" ].keys() ):
        print ( "State %4s: %-s" %( key, data[ "Atomic states" ][ key ]) )

def get_number_of_kpoinst( data ):
    return data[ 'nkstot' ]

def get_number_of_bands( data ):
    return data[ 'nbnd' ]

def get_nspin( data ):
    return data[ 'nspin' ]

def sum_states( data, spin = 'up', state_index = [ 0 ] ):
    print ( "Suming the following projections:" )
    print ( "%-8s %-s" %( 'Index', 'Projection' ) )
    for key in state_index:
        print ( "%-8s %-s" %( key, data[ "Atomic states" ][ key ]) )

    nkstotal = get_number_of_kpoinst( data )
    nbnd     = get_number_of_bands( data )
    nspin    = get_nspin( data )
    sum_data = np.zeros( [ nspin * nkstotal * nbnd, 3 ] ) # 3 columns for ik, e, sum
    ik_start = 0
    ik_end   = nkstotal

    if spin == 'down' and nspin == 2:
       ik_start = nkstotal
       ik_end   = 2 * nkstotal
    elif nspin == 1 and spin == 'down':
       print ( "Output file indicates non-spin-polarized calculation" )
       print ( "spin = \"down\" is ignored" )

    if spin == "up" and nspin == 2:
       print ( "Output file indicates spin-polarized calculation" )
       print ( "Use spin = \"down\" for spin-down data" )
       
    #for ik in sorted( data[ 'atomic projection' ].keys():
    for ik in range( ik_start, ik_end ):
        for ie in sorted( data[ 'atomic projection' ][ ik ]['eigen'].keys() ):
             proj = data[ 'atomic projection' ][ ik ]['eigen'][ ie ][ 'psi' ]
             pdos = 0
             for idx in state_index:
                 pdos += proj[ idx ]
             e = data[ 'atomic projection' ][ ik ]['eigen'][ ie ][ 'e' ]
             sum_data[ ik * nbnd + ie, : ] =  np.array( [ ik, e, pdos ] )
    return sum_data

def parse( output, args = {} ):
    lines = read_file( output )
    data = {}
    data = read_keys( data, lines, args )
    data = read_atomic_states( data, lines ) 
    data = read_projection( data, lines )
    data = read_lowdin( data, lines )
  
    return data


