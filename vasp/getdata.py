import gzip as gz

def get_energy( outcar ):
    if '.gz' in outcar:
        with gz.open( outcar, 'r' ) as f:
            lines = f.readlines( )
    else:
        with open( outcar, 'r' ) as f:
            lines = f.readlines( )
    for iline in range( len( lines ) - 1, -1, -1 ):
        if "energy  w" in lines[ iline ]:
            ener = float( lines[ iline ].split( )[ 6 ] )
            break
    return ener
