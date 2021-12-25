import os, json
import gzip  as gz
import numpy as np
import pandas as pd 
from ..small_tools import str_decode

def parse_doscar( doscar, parse_pdos = False, assigned_characters = 's py pz px dxy dyz dz2 dxz dx2-y2', nline = 1 ):
    """
    Parsing a DOSCAR file.
  
    Input:
    ------
    doscar: str
    Path to a DOSCAR file. 

    parse_pdos: bool
    Parsing the atom projected dos

    character: str
    Characters of projected dos. This should be supplied to by user, 
    because number of columns and characters depend on setting in INCAR

    nline: int
    Number of lines for one energy point. Default is 1, good for many cases. 
    However, in some cases, the line is long and is breaked into many parts.
  
    Returns:
    --------
    data: dict
    Dictionary for dos and pdos data
 
    """
 
    print ( "Parsing %s" %( doscar ) )
    try:
        if '.gz' in doscar:
            with gz.open( doscar, 'r' ) as f:
                lines = f.readlines( )
            decode = True
        else:
            with open( doscar, 'r' ) as f:
                lines = f.readlines( )
            decode = False
    except:
        print ( "Error with reading %s" %( doscar ) )
        exit( ) 

    data = { }
    ( pdos, nions, efermi, npts, emin, emax, keyline ) = parse_info( lines, decode )
    print ( pdos, nions, efermi, npts, emin, emax, keyline )
    de = ( emax - emin ) / float( npts - 1 )
    data[ 'ENERGY' ] = np.arange( emin, emax, de ).tolist( ) + [ emax ] 
    data[ 'EFERMI' ] = efermi
    data[ 'NIONS'  ] = nions
    if not parse_pdos:
        energy, total = parse_total_dos( lines, npts, decode )
        if total.shape[ 1 ]   == 2:
            data[ 'tot' ]    = total[ :, 0].tolist( )
            data[ 'tot_up' ] = ( total[ :, 0] / 2.) .tolist( )
            data[ 'tot_dn' ] = ( total[ :, 0] / 2.) .tolist( )
        elif total.shape[ 1 ] == 4:
            data[ 'tot_up' ] = total[ :, 0].tolist( )
            data[ 'tot_dn' ] = total[ :, 1].tolist( )
            data[ 'tot' ]    = (total[ :, 0] + total[ :, 1 ]).tolist( )
        else:
            raise( 'DOSCAR format was not implemented' )

    if parse_pdos:
        if pdos == 0:
            raise( 'DOSCAR file does not contain partial DOS' )
        print ( "Assuming the following characters %s" %( assigned_characters ) )
        for iion in range( 1, nions + 1 ):
            pdos = parse_pdos_ion( lines, npts, iion, nline  = nline, decode = decode)
            pdos_dict = { }
            #print pdos.shape
            for i in range( pdos.shape[ 1 ] ):
                pdos_dict[ assigned_characters.split( )[ i ] ] = pdos[ :, i ].tolist()
            data[ iion ] = pdos_dict
     
    return data 

def parse_info( lines, decode = True ):
    line   = lines[ 0 ]
    pdos   = int( line.split( )[ 2 ] ) 
    nions  = int( line.split( )[ 0 ] )
    line   = str_decode( lines[ 5 ], decode )
    # Fix error with large NEDOS
    l1     = line[ :32 ].split()
    l2     = line[ 32: ].split()
    l      = l1 + l2
    line   = ' '.join( l )
    efermi = float( line.split( )[ 3 ] )
    npts   = int( line.split( )[ 2 ] )
    emax   = float( line.split( )[ 0 ] )
    emin   = float( line.split( )[ 1 ] )
    keyline = line.strip( )
    return ( pdos, nions, efermi, npts, emin, emax, keyline )

def parse_total_dos( lines, npts, decode = True ):
    ener  = [ ]
    total = [ ]
    for iline in range( 6, 6 + npts ):
        ener.append( float( lines[ iline ].split( )[ 0 ] ) )
        total.append( [ float( x ) for x in lines[ iline ].split( )[ 1: ] ] )
    return np.array( ener ), np.array( total )

def parse_pdos_ion( lines, npts, iion, nline = 1, decode = True ):
    pdos = [ ]
    startline = 6 + npts + 1 + ( iion - 1 ) * ( nline * npts + 1 )
    endline   = startline + nline * npts
    for iline in range( startline, endline, nline ):
        line = ''
        for i in range( nline ):
            line += ' ' + str_decode( lines[ iline + i ].strip(), decode )
        pdos.append( [ float( x ) for x in line.split( )[ 1: ] ] )
    return np.array( pdos )


def get_total_dos( doscars, fout ):
    """
    Get total DOS from list of DOSCAR files

    Input:
    ------
    doscars: dict
    is a dictionary containing path and weight to all DOSCAR
    doscars[ i ] = { 'doscar': 'PATH TO DOSCAR'; 'weight' : weight }

    fout: str
    Filename of output

    Return:
    -------
    energy: 1D-array
    dos: 1D-array
    """
    tot          = 0 
    tot_up       = 0 
    tot_dn       = 0 
    total_weight = 0
    for i in sorted( doscars.keys( ) ):
        if 'data' not in doscars[ i ].keys():
            data = parse_doscar( doscars[ i ][ 'doscar' ] )
        else:
            data = doscars[ i ][ 'data' ]
 
        energy        = data[ 'ENERGY' ]
        weight        = doscars[ i ][ 'weight' ]
        tot          += np.array( data[ 'tot' ] ) * weight
        tot_up       += np.array( data[ 'tot_up' ] ) * weight
        tot_dn       += np.array( data[ 'tot_dn' ] ) * weight
        total_weight += weight
            
    tot    /= float( total_weight )
    tot_up /= float( total_weight )
    tot_dn /= float( total_weight )
    df = pd.DataFrame()
    df[ 'Energy' ] = energy
    df[ 'TOTAL' ]  = tot
    df[ 'TOTAL_UP' ]  = tot_up
    df[ 'TOTAL_DN' ]  = tot_dn

    os.system( 'mkdir -p doscar_tmp' )
    #with open( 'doscar_tmp/' + fout, 'w' ) as f:
    #    for z in zip( energy, tot, tot_up, tot_dn ):
    #        f.write( "%12.8f %12.8f %12.8f %12.8f\n" %( z[ 0 ], z[ 1 ], z[ 2 ], z[ 3 ]) )
    df.to_csv( 'doscar_tmp/' + fout, index = False )
    return energy, tot, tot_up, tot_dn

def get_pdos( doscars, fout, iions = [ ], out_characters = [ ], save_data = False, assigned_characters = 's py pz px dxy dyz dz2 dxz dx2-y2', nline = 1 ):
    """
    Parsing and summing dos from multiple DOSCAR files.
  
    Input:
    ------
    doscars: dict
    is a dictionary containing path and weight to all DOSCAR
    doscars[ i ] = { 'doscar': 'PATH TO DOSCAR'; 'weight' : weight }

    fout: str
    Filename of output

    iions: list 
    A list of index of ions that will be summed together for pdos

    out_characters: list 
    A list of characters that will be summed together for pdos

    save_data: bool
    If true, parsed pdos data will be save in doscars[ i ]

    assigned_characters: str
    List of characters that will be assigned to data.
    example:
    Most cases: s py pz px dxy dyz dz2 dxz dx2-y2
    spin_pol:   s_up s_dn py_up py_dn pz_up pz_dn px_up px_dn ...
 
    nline: int
    Number of lines for one energy point. Default is 1, good for many cases. 
    However, in some cases, the line is long and is breaked into many parts.
  
    Returns:
    --------
    energy: 1D-array
    pdos: 1D-array
    doscars: dict
    Updated doscars with new data, if save_data = True
 
    """
    pdos = 0
    total_weight = 0

    for i in sorted( doscars.keys( ) ):
        if 'data' not in doscars[ i ].keys( ): 
            data   = parse_doscar( doscars[ i ][ 'doscar' ], parse_pdos = True, assigned_characters = assigned_characters, nline = nline)
            if save_data:
                doscars[ i ][ 'data' ] = data
        else:
            data = doscars[ i ][ 'data' ]

        energy = data[ 'ENERGY' ]
        weight = doscars[ i ][ 'weight' ]
        total_weight += weight

        if len( iions ) == 0:
            iions = range( 1, data[ 'NIONS' ] + 1 )

        if len( assigned_characters ) == 0:
            assigned_characters = data[ iions[ 0 ] ].keys( )

        pdos_tmp = [ ]
        for iion in iions:
            for character in out_characters:
                pdos_tmp.append( data[ iion ][ character ] )

        pdos_tmp  = np.abs( np.array( pdos_tmp ) ).sum( axis = 0 )
        pdos += pdos_tmp * weight
    pdos /= float( total_weight )
    os.system( 'mkdir -p doscar_tmp' )
    with open( 'doscar_tmp/' + fout, 'w' ) as f:
        for z in zip( energy, pdos ):
            f.write( "%12.8f %12.8f \n" %( z[ 0 ], z[ 1 ] ) )
   
    return energy, pdos, doscars 


def get_atom_pdos( doscars, fout, iion, save_data = False, assigned_characters = 's py pz px dxy dyz dz2 dxz dx2-y2', nline = 1 ):
    """
    Parsing pdos from multiple DOSCAR files for each iion, and save data for each atom in doscar_tmp/fout.
  
    Input:
    ------
    doscars: dict
    is a dictionary containing path and weight to all DOSCAR
    doscars[ i ] = { 'doscar': 'PATH TO DOSCAR'; 'weight' : weight }

    fout: str
    Filename of output

    iion: integer 
    Index of ion for pdos

    save_data: bool
    If true, parsed pdos data will be save in doscars[ i ]

    assigned_characters: str
    List of characters that will be assigned to data.
    example:
    Most cases: s py pz px dxy dyz dz2 dxz dx2-y2
    spin_pol:   s_up s_dn py_up py_dn pz_up pz_dn px_up px_dn ...
 
    nline: int
    Number of lines for one energy point. Default is 1, good for many cases. 
    However, in some cases, the line is long and is breaked into many parts.
  
    Returns:
    --------
    energy: 1D-array
    pdos: array contain pdos for atom
    doscars: dict
    Updated doscars with new data, if save_data = True
 
    """
    total_weight = 0

    pdos = pd.DataFrame( )

    for i in sorted( doscars.keys( ) ):
        if 'data' not in doscars[ i ].keys( ):
            data   = parse_doscar( doscars[ i ][ 'doscar' ], parse_pdos = True, assigned_characters = assigned_characters, nline = nline)
            if save_data:
                doscars[ i ][ 'data' ] = data
        else:
            data = doscars[ i ][ 'data' ]

        energy = data[ 'ENERGY' ]
        weight = doscars[ i ][ 'weight' ]
        total_weight += weight

        if len( assigned_characters ) == 0:
            assigned_characters = data[ iion ].keys( )

        #pdos_tmp = [ ]
        for character in assigned_characters.split():
            #pdos_tmp.append( data[ iion ][ character ] 
            if character in pdos.keys():
                pdos[ character ] += np.array( data[ iion ][ character ] ) * weight
            else:
                pdos[ character ]  = np.array( data[ iion ][ character ] ) * weight

    for character in assigned_characters.split():
        pdos[ character ] /= float( total_weight )
  
    pdos[ 'Energy' ] = np.array( energy )

    os.system( 'mkdir -p doscar_tmp' )

    pdos.to_csv( 'doscar_tmp/' + fout, index = False ) 

    return energy, pdos, doscars

def read_kpts( top_dir = './' ):
    if not os.path.isfile( top_dir + 'TMP/IBZKPT' ):
        raise( 'File ' + top_dir + 'TMP/IBZKPT' + ' not exits')
    with open( top_dir + 'TMP/IBZKPT', 'r' ) as f:
        lines = f.readlines( )

    nkpts = int(  lines[ 1 ].strip() )
    allkpts = [ ]
    for i in range( 3, 3 + nkpts ):
        allkpts.append( lines[ i ].split( ) )
    return allkpts

def write_kpoints( fout, kpts, top_dir = './' ):
    with open( top_dir + 'TMP/IBZKPT', 'r' ) as f:
        lines = f.readlines( )

    fout.write( lines[ 0 ].strip( ) + '\n')
    fout.write( str( len( kpts )  ) + '\n')
    fout.write( lines[ 2 ].strip( ) + '\n')
    for kpt in kpts:
        s = ' '
        fout.write( s.join( x for x in kpt ) + '\n')

def copy_file( ifile, top_dir = "./" ):
    file_list = [ 'IBZKPT', 'INCAR', 'POTCAR', 'POSCAR', 'job' ]
    for file in file_list:
        if not os.path.isfile( top_dir + 'TMP/' + file ):
            raise( 'File ' + top_dir + 'TMP/' + file + ' not exits')

    os.system( 'cp -a ' + top_dir + '/TMP K.' + str( ifile ) )
    os.system( 'cat ' + top_dir + '/TMP/job | sed \'s/$ifile/' + str( ifile ) + '/\' >  ' + top_dir + '/K.' + str( ifile ) + '/job' )

def gen_input( nkpfile, top_dir = './' ):
    allkpts = read_kpts( top_dir = './' )

    kid = range( 0, len( allkpts ), nkpfile )

    for ifile in range( len( kid ) ):
        begin = kid[ ifile ]
        if ifile < len( kid ) - 1:
            end   = kid[ ifile + 1]
        else:
            end   = len( allkpts )
        kpts = [ ]
        for i in range( begin, end ):
            kpts.append( allkpts[ i ] )
        copy_file( ifile )
        with open( top_dir + '/K.' + str( ifile ) + '/KPOINTS', 'w' ) as fout:
            write_kpoints( fout, kpts, top_dir )

def gen_doscars( nkpfile, top_dir = './' ):

    allkpts = read_kpts( top_dir )
    totalweight = np.sum( np.array( [ int( k[ -1 ] ) for k in allkpts ] ) )
    kid = range( 0, len( allkpts ), nkpfile )

    doscars = {}
    for ifile in range( len( kid ) ):
        begin = kid[ ifile ]
        if ifile < len( kid ) - 1:
            end   = kid[ ifile + 1]
        else:
            end   = len( allkpts )
        kpts = [ ]
        for i in range( begin, end ):
            kpts.append( allkpts[ i ] )
        weight = np.sum( np.array( [ int( k[ -1 ] ) for k in kpts ] ) )
        doscars[ ifile ] = { }
        if os.path.isfile( top_dir + 'K.' + str( ifile ) + '/DOSCAR' ):
            file_name = top_dir + 'K.' + str( ifile ) + '/DOSCAR'
        elif os.path.isfile( top_dir + 'K.' + str( ifile ) + '/DOSCAR.gz' ):
            file_name = top_dir + 'K.' + str( ifile ) + '/DOSCAR.gz'
        else:
            raise( 'NO DOSCAR FILE FOUND IN ' + 'K.' + str( ifile ) )

        doscars[ ifile ][ 'doscar' ] = file_name   
        doscars[ ifile ][ 'weight' ] = weight

    return doscars
     
