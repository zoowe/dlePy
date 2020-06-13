"""
****************************************
Author: Duy Le
Department of Physics
University of Central Florida

Email: duy.le@ucf.edu
Website: http://www.physics.ufc.edu/~dle
****************************************
"""

from ase    import Atoms
from ase.io import read, write
import os
import numpy as np

def rotate_group( atoms, group, axis, angle, center ):
    """ Rotate a group of atoms
   
    Parameters
    ----------
    atoms:  ase.Atoms object
         The system of which a group of atoms need to be rotated
    group:  list
         List of indices of atoms of a group
    axis: 3 element array
         Vector indicating the direction of rotation axis 
    angle: Scalar (degree)
         Rotate angle
    center: 3 element array
         Position of rotation center 
        
    Returns
    -------
    atoms: ase.Atoms object
         Modified ase.Atoms object with group rotated.

    Examples
    --------
 
    See also
    --------
    
    Update
    ------
    10/25/2016: duy.le@ucf.edu started this routine 
    """

    subsys = atoms[ group ]

    subsys.rotate( angle, axis, center, rotate_cell = False )

    atoms.positions[ group ] = subsys.positions

    return atoms

def translate_group( atoms, group, dir ):
    """ Translate a group of atoms
   
    Parameters
    ----------
    atoms:  ase.Atoms object
         The system of which a group of atoms need to be rotated
    group:  list
         List of indices of atoms of a group
    dir: 3 element array
         Vector indicating the direction of translation 
        
    Returns
    -------
    atoms: ase.Atoms object
         Modified ase.Atoms object with group rotated.

    Examples
    --------
 
    See also
    --------
    
    Update
    ------
    10/25/2016: duy.le@ucf.edu started this routine 
    """

    subsys = atoms[ group ]

    subsys.positions[ :, 0 ] += dir[ 0 ]
    subsys.positions[ :, 1 ] += dir[ 1 ]
    subsys.positions[ :, 2 ] += dir[ 2 ]

    atoms.positions[ group ] = subsys.positions

    return atoms


def sort_group( atoms, indices = [], symbols = [] ):
    """ Sort only the atoms with indices in `indices` or symbol in `symbols`

    Parameters
    ----------
    atoms:  ase.Atoms object
         The system of which a group of atoms need to be rotated
    indices:  list of indices
         List of indices of atoms for sorting  
    symbols:  list of symbols
         List of symbols of atoms for sorting

    Returns
    -------

    Update
    ------

    """

    if ( len( indices ) > 0 and len( symbols ) > 0 ) or ( len( indices ) == 0 and len( symbols ) == 0 ):
         raise RuntimeError( 'ERROR: indices or symbols, but not both,  must be specified' )
    
    atomsortlist = [ ]
    atomfixlist  = [ ] 

    if len( symbols ) > 0:

        for at in range( len( atoms ) -1, -1, -1 ):
            if atoms[ at ].symbol in symbols:
                atomsortlist.append( at )
            else:   
                atomfixlist.append( at )

    if len( indices ) > 0:

        for at in range( len( atoms ) ):
            if atoms[ at ].index in indices:
                atomsortlist.append( at )
            else:
                atomfixlist.append( at )


    sortatoms = atoms[ atomsortlist ]

    # Use sort=True to write it for now
    write('sort.POSCAR', sortatoms, format = 'vasp', direct = True, sort = True )

    # Read `sortatoms`
    sortatoms = read ( 'sort.POSCAR', format = 'vasp' )
    os.system( 'rm -f sort.POSCAR' )

    if len( atomfixlist ) > 0:
        fixatoms  = atoms[ atomfixlist  ]
        atoms_ = fixatoms + sortatoms
    else:
        atoms_ = sortatoms

    return atoms_ 

def sort_layers( system, layer = [ ] ):
    """
    This routine will create a proper POSCAR for doing partial phonopy run. 
    It will group atoms in different group, from bottom to top, so it will
    be convenient for seprate them later.

    system: ase.Atoms object
    layer: List of all list of atom.index in all layers.
    example: 
    layer = [ 
              [ 0, 1, 2 ],
              [ 3, 4, 5, 6, 7, 8],
              [ 9, 10, 11, 12]
            ]
    This will put [ 9, 10, 11, 12 ] atoms to the top of list, 
                  [ 3, 4, 5, 6, 7, 8] to the middle,
              and [ 0, 1, 2 ] to the top of the atoms list.

    """
    if system.get_number_of_atoms() != np.sum( len( x ) for x in layer ):
        raise RuntimeError("""
        Number of atoms are different.
        system.get_number_of_atoms( ) != np.sum( len( x ) for x in layer )
        %s VS %s
        """ %( system.get_number_of_atoms( ), np.sum( len( x ) for x in layer ) ) )

    if len( layer ) < 2:
        raise RuntimeError( """
        len( layer ) much be larger than 1
        """ )

    atoms = system[ layer[ -1 ] ]
    for ilayer in range( len( layer ) - 2, -1, -1):
        atoms += system[ layer[ ilayer ] ]

    return atoms


### For NEB
def move_atom( system, ats, dirs, ds ):
    '''
    Move a list of atoms along dirs, distances (ds) (in scaled coors
    ats  = [ 1,  2    3 ]
    dirs = [ 0,  1,   2 ]
    ds   = [ 1, -1, 0.5 ]
    system = move_atoms( system, ats, dirs, ds )
    '''
    scaled_coors =system.get_scaled_positions()

    for i in range( len( ats ) ):
        scaled_coors[ ats[ i ], dirs[ i ]] = scaled_coors[ ats[ i ], dirs[ i ]] + ds[ i ]

    system.set_scaled_positions( scaled_coors )

    return system

def scaled_position_positive( atoms ):
    scaled_coor  = atoms.get_scaled_positions( )
    print ( 'Inside:', scaled_coor[ 4 ] )
    for iat in range( len( atoms ) ):
        for i in range( 3 ):
            if scaled_coor[ iat ][ i ] < 0:
                scaled_coor[ iat ]  [ i ] += 1.
    atoms.set_scaled_positions( scaled_coor ) 
    return atoms
    

def regulate_pbc( init_sys, final_sys ):
    '''
    Compare position of atoms in two system. If any that move about 1 
    lattice cell, it is shifted so that interpolation between 
    init_sys and final_sys gives good path for NEB.

    '''
    scaled_coor_init  = init_sys.get_scaled_positions( )
    print ( scaled_coor_init[ 4 ] )
    init_sys          = scaled_position_positive( init_sys )
    scaled_coor_init  = init_sys.get_scaled_positions( )
    print ( scaled_coor_init[ 4 ] )

    final_sys         = scaled_position_positive( final_sys )
    scaled_coor_init  = init_sys.get_scaled_positions( )
    scaled_coor_final = final_sys.get_scaled_positions( )
    scaled_coor_diff  = scaled_coor_final - scaled_coor_init
    list = []
    for at in range( len( init_sys ) ):
        if scaled_coor_diff[ at, 0 ] > 0.9:
            list.append( [ at, 0, -1 ] )
        if scaled_coor_diff[ at, 0 ] < -0.9:
            list.append( [ at, 0, 1 ] )
        if scaled_coor_diff[ at, 1 ] > 0.9:
            list.append( [ at, 1, -1 ] )
        if scaled_coor_diff[ at, 1 ] < -0.9:
            list.append( [ at, 1, 1 ] )
        if scaled_coor_diff[ at, 2 ] > 0.9:
            list.append( [ at, 2, -1 ] )
        if scaled_coor_diff[ at, 2 ] < -0.9:
            list.append( [ at, 2, 1 ] )

    print ( list )
    ats   = [ item[ 0 ] for item in list ]
    dirs  = [ item[ 1 ] for item in list ]
    moves = [ item[ 2 ] for item in list ]

    final_sys = move_atom( final_sys, ats, dirs, moves )
    return init_sys, final_sys

