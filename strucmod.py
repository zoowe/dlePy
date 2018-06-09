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
    symbols:  list of symbols
         List of symbols of 

    Returns
    -------

    Update
    ------

    """

    if ( len( indices ) > 0 and len( symbols ) > 0 ) or ( len( indices ) == 0 and len( symbols ) == 0 ):
         exit( 'ERROR: indices or symbols, but not both,  must be specified' )
    
    atomsortlist = [ ]
    atomfixlist  = [ ] 

    if len( symbols ) > 0:

        for at in range( len( atoms ) -1, -1, -1 ):
            if atoms[ at ].symbol in symbols:
                atomsortlist.append( at )
            else:   
                atomfixlist.append( at )

    if len( indices ) > 0:

        for at in range( len( atoms ) -1, -1, -1 ):
            if atoms[ at ].index in indices:
                atomsortlist.append( at )
            else:
                atomfixlist.append( at )


    sortatoms = atoms[ atomsortlist ]
    fixatoms  = atoms[ atomfixlist  ]

    # Use sort=True to write it for now
    write('sort.POSCAR', sortatoms, format = 'vasp', direct = True, sort = True )

    # Read `sortatoms`
    sortatoms = read ( 'sort.POSCAR', format = 'vasp' )
    os.system( 'rm -f sort.POSCAR' )

    atoms_ = fixatoms + sortatoms

    return atoms_ 

