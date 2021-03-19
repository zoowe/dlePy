from yaml import load
import numpy as np
import json

def yaml2json( f_bandyaml ):
    f = open ( f_bandyaml, 'r' )
    data = load( f )
    f.close ()
    
    fout = 'band.dict'
    with open( fout, 'w' ) as f:
        f.write( json.dumps( data ) )
    os.system( 'gzip ' + fout )
    

def parse_band( f_bandyaml ):
    f = open ( f_bandyaml, 'r' )
    data = load( f )
    f.close ()

    nqpoint  = data[ 'nqpoint' ]
    natom    = data[ 'natom' ]
    nband    = 3 * natom
    distance = np.zeros( [ nqpoint ] )
    label    = [ ]
    label_coor = [ ]
    freq     = np.zeros( [ nqpoint, nband ] )

    for iq in range( nqpoint ):
        phonon         = data[ 'phonon' ][ iq ]
        distance[ iq ] = phonon[  'distance' ]  
        if 'label' in phonon.keys( ):
            label.append( phonon[  'label' ] )
            label_coor.append( phonon[  'distance' ] )

        for iband in range( nband ):
            freq[ iq, iband ] = phonon[ 'band' ][ iband ][ 'frequency' ]

    return ( distance, freq, label, label_coor )

def parse_thermo( f_thermoyaml ):
    f = open ( f_bandyaml, 'r' )
    data = load( f )
    f.close ()

    nqpoint  = data[ 'nqpoint' ]
    natom    = data[ 'natom' ]
    nband    = 3 * natom
    distance = np.zeros( [ nqpoint ] )
    label    = [ ]
    label_coor = [ ]
    freq     = np.zeros( [ nqpoint, nband ] )

    for iq in range( nqpoint ):
        phonon         = data[ 'phonon' ][ iq ]
        distance[ iq ] = phonon[  'distance' ]
        if 'label' in phonon.keys( ):
            label.append( phonon[  'label' ] )
            label_coor.append( phonon[  'distance' ] )

        for iband in range( nband ):
            freq[ iq, iband ] = phonon[ 'band' ][ iband ][ 'frequency' ]

    return ( distance, freq, label, label_coor )


def FvibFunc( T, f ):
    # Equation 22, ab initio thermodynamics: a primer
    # T[K], omega[THz], output in eV
    hbar = 6.582119514 * 10 ** ( -16 ) # eV*s
    kB   = 8.6173303  *  10 ** ( -5 )  # eV K-1
    THz2eV = 1. / 241.79893
    omega = f * THz2eV
    ZPE = omega / 2.
    secondterm = 0.
    if T > 0:
        beta = kB * T
        secondterm = beta * np.log ( 1 - np.exp( - omega / beta ) )
    tmp = ZPE + secondterm
    return tmp
        
def get_Fvib( f, pdos, T ):
    df= f[ -1 ] - f[ -2 ]
    # Make every negative mode zero dos:
    #for i in range( len( f ) ):
    #    if f[ i ] < 10 ** (-5):
    #        pdos[ i ] = 0.
    #        f   [ i ] = 10 ** ( -5 )
    #total_pdos = np.sum( pdos * df )
    #print total_pdos
    #norm  = nmodes / total_pdos
    #print norm
    #pdos *= norm
    Fvib  = np.sum ( FvibFunc( T, f ) * pdos * df )

    return Fvib
