import gzip as gz
from dlePy.math import str2bool

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
        print "WARNING: perhaps OUTCAR is for an imcompleted calculation"
        print "Please check: ", outcar

    return lines

def get_calculation_details( outcar ):
    
    lines = get_lines_outcar( outcar )

    for iline in range( len( lines ) - 1, -1, -1 ):
  
        if "PREC " in lines[ iline ]:
            PREC = lines[ iline ].split( )[ 2 ] 

        if "ENCUT " in lines[ iline ]:
            ENCUT = float( lines[ iline ].split( )[ 2 ] )
      
        if "LREAL " in lines[ iline ]:
            LREAL = str2bool( lines[ iline ].split( )[ 2 ] )

        if "EDIFF  =" in lines[ iline ]:
            EDIFF = float( lines[ iline ].split( )[ 2 ] )

        if "EDIFFG " in lines[ iline ]:
            EDIFFG = float( lines[ iline ].split( )[ 2 ] )

        if "ISMEAR =" in lines[ iline ]:
            ISMEAR = int( lines[ iline ].split( )[ 2 ].replace( ';', '' ) ) 
            SIGMA  = float( lines[ iline ].split( )[ 5 ] )

        if "NBANDS=" in lines[ iline ]:
            NKPTS = int( lines[ iline ].split( )[ 3 ].replace( ';', '' ) )
            NBANDS = int( lines[ iline ].split( )[ 14 ] )

        if "NELECT =" in lines[ iline ]:
            NELECT = float( lines[ iline ].split( )[ 2 ].replace( ';', '' ) )


    details = { 
                "PREC"       : PREC,
                "LREAL"      : LREAL,
                "ENCUT"      : ENCUT,
                "EDIFF"      : EDIFF,
                "EDIFFG"     : EDIFFG,
                "ISMEAR"     : ISMEAR,
                "SIGMA"      : SIGMA,
                "NKPTS"      : NKPTS,
                "NBANDS"     : NBANDS,
                "NELECT"     : NELECT
              }

    del lines

    return details
    

def get_energy( outcar ):

    lines = get_lines_outcar( outcar )

    for iline in range( len( lines ) - 1, -1, -1 ):
        if "energy  w" in lines[ iline ]:
            ener = float( lines[ iline ].split( )[ 6 ] )
            break

    del lines

    return ener

def get_efermi( outcar ):

    lines = get_lines_outcar( outcar )

    for iline in range( len( lines ) - 1, -1, -1 ):
        if "E-fermi" in lines[ iline ]:
            efermi = float( lines[ iline ].split( )[ 2 ] )
            break

    del lines

    return efermi

