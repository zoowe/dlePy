import gzip as gz

def get_lines_outcar( outcar ):
    if '.gz' in outcar:
        with gz.open( outcar, 'r' ) as f:
            lines = f.readlines( )
    else:
        with open( outcar, 'r' ) as f:
            lines = f.readlines( )
    #Check if calculation is done
    job_done = False
    for iline in range( len( lines ) - 1, -1, -1 ):
        if 'reached required accuracy - stopping structural energy minimisation' in lines[ iline ]:
            job_done = True
            break
    if not job_done:
        print "WARNING: perhaps OUTCAR is for an imcompleted calculations"
        print "Please check: ", outcar

    return lines

def get_energy( outcar ):

    lines = get_lines_outcar( outcar )

    for iline in range( len( lines ) - 1, -1, -1 ):
        if "energy  w" in lines[ iline ]:
            ener = float( lines[ iline ].split( )[ 6 ] )
            break

    return ener

def get_efermi( outcar ):

    lines = get_lines_outcar( outcar )

    for iline in range( len( lines ) - 1, -1, -1 ):
        if "E-fermi" in lines[ iline ]:
            efermi = float( lines[ iline ].split( )[ 2 ] )
            break

    return efermi

