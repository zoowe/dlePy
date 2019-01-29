import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib.ticker as ticker
import numpy as np
import gzip as gz
import os

def read_file( filename ):
    if not os.path.isfile( filename ):   
        print "%s is not found" %( filename )
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

def read_keys( data, lines ):
    list_key = [ 'natomwfc', 'nbnd', 'nkstot',  'npwx', 'nkb' ]

    for i in range( len( lines ) ):
       count_key = 0
       for name in list_key:
           if name in lines[ i ]:
               data = add_key( data, name, lines[ i ] )
               count_key += 1
               if count_key == len( list_key ):
                   break

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
            for j in range( i + 1, i + 10 * data[ 'nbnd' ] ):
                if 'e =' in lines[ j ]:
                    band = {}
                    e = float( lines[ j ].split()[ -2 ] )
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
            if ik == data[ 'nkstot' ]:
                break 
    data[ 'atomic projection' ] = kdata

    return data    

def parse( lines ):
    data = {}
    
    data = read_keys( data, lines )

    data = read_atomic_states( data, lines ) 
     
    data = read_projection( data, lines )
  
    return data


if __name__ == "__main__":
    loc = '/home/duy/PROJ_MetalOnSemiconductor/andre/TUTORIALS/DensityofStates/Dataaccumalation/k48/thesmear.020/u=2/'
    filename = 'outputpdos.dat'
    lines = read_file( loc + filename )
    data = parse( lines )
    nkstotal = data[ 'nkstotal' ]
    nbnd     = data[ 'nbnd' ]    
    band_data = np.zeros( [ nbnd * nbnd, 4 ] ) # 4 columns for ik, e, s, p
    for ik in sorted( data[ 'atomic projection' ].keys() ):
        for ie in sorted( data[ 'atomic projection' ][ ik ]['eigen'].keys() ):
             proj = data[ 'atomic projection' ][ ik ]['eigen'][ ie ][ 'psi' ]
             s = proj[ 0 ] + proj[ 4 ]
             p = proj[ 1 ] + proj[ 2 ] + proj[ 3 ] +  proj[ 5 ] + proj[ 6 ] + proj[ 7 ] 
             e = data[ 'atomic projection' ][ ik ]['eigen'][ ie ][ 'e' ]
             band_data[ ik * 16 + ie, : ] =  np.array( [ ik, e, s, p ] )
