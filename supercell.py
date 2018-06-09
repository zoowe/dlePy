"""
****************************************
Author: Duy Le
Department of Physics
University of Central Florida

Email: duy.le@ucf.edu
Website: http://www.physics.ufc.edu/~dle
****************************************
"""

#!/usr/bin/env python
from   ase.io import read
import numpy as np
import os
from dlePy.math import dotproduct, length, area


def supercell( atoms, n0, n1, n2 ):
    '''
    Create supercell from primitive cell `atoms`
    This method is better than atoms * ( n0, n1, n2 ) because it keeps the order of atoms.
    '''

    org_atoms = atoms.copy( )
 
    new_atoms = sub_atoms = org_atoms[ [ 0 ] ] * ( n0, n1, n2 )
    for i in range( 1, org_atoms.get_number_of_atoms() ):
       new_atoms += org_atoms[ [ i ] ] * ( n0, n1, n2 ) 

    return new_atoms 


def remove_atoms_outside( atoms, shift = 0., threshold = 0.01 ):
    """
    Remove atoms outside supercell.
    
    atoms: ase atoms object

    shift: scalar. 
    It is used to calculate a shift applied to positions of atoms.
    It is helpful if atoms at the boder of supercell.

    threshold: scalar
    This routine calculate the diffrence between coordiate of the same atom
    using wrapped and unwrapped method. If the difference is larger than threshold, 
    atom is outside the supercell, then it is removed.

    """

    sys0 = atoms.copy( )

    # Shift whole system by delta. This is important if there are atoms
    # right on the cell boundary
    delta = shift * ( sys0.cell[ 0 ] + sys0.cell[ 1 ] + sys0.cell[ 2 ] )
    sys0.positions[ :, 0 ] += delta[ 0 ]
    sys0.positions[ :, 1 ] += delta[ 1 ]
    sys0.positions[ :, 2 ] += delta[ 2 ]

    # Here is the trick: get position of atom with wrap=True/False
    # If they are different, the atom is outside the cell.
    for at in range( len( sys0 ) -1, -1, -1 ):
        wrapped_pos   = sys0.get_positions( wrap = True )[ at ]
        unwrapped_pos = sys0.get_positions( wrap = False )[ at ]

        # Delete all atoms out side the cell
        if length(wrapped_pos - unwrapped_pos ) > threshold:
            del sys0[at]

    # Shift atoms back to original coordinates.   
    sys0.positions[ :, 0 ] -= delta[ 0 ]
    sys0.positions[ :, 1 ] -= delta[ 1 ]
    sys0.positions[ :, 2 ] -= delta[ 2 ]
    
    return sys0

def create_matrix_surface( atoms, matrix = ( 1, 0, 0, 1 ), pad=1, shift = 0.01, threshold = 0.01 ):
    """
    Create surface using matrix notation
    | m1  n1 |
    | m2  n2 |
    maxtrix = ( m1, n1, m2, n2 )

    pad: sometimes some atoms are missing. Increase pad to avoid this problem
    Shift and threshold: are for remove_atoms_outside
    """

    sys0= atoms.copy( )
    a1=sys0.cell[ 0, : ]
    a2=sys0.cell[ 1, : ]

    (m1, n1, m2, n2)  = matrix
   
    # New surface cell vector.
    A1_surf = a1 * m1 + a2 * n1
    A2_surf = a1 * m2 + a2 * n2

    # Calculate max length of A1, A2 or the diagonal
    max_length = np.max( ( length( A1_surf ), length( A2_surf ), length( A1_surf + A2_surf ) ) )

    # Estimate sample size from max_length. If result is not good, increase pad
    sample_size = ( 2 * int( max_length / length( a1 )  + pad ), 2 * int( max_length / length( a2 )  + pad ) )

    # Create large enough sample to cut the new supercell
    surface= supercell( sys0, sample_size[ 0 ], sample_size[ 1 ] , 1)

    #Make atom at center of cell become (0,0)
    # Find atom:
    indx  =  int( len( surface ) / 2 )
    for at in range( len( surface ) ):
        if     abs( surface.get_scaled_positions( )[ at ][ 0 ] - 0.5 ) < 0.1:
            if abs( surface.get_scaled_positions( )[ at ][ 1 ] - 0.5 ) < 0.1:
                # Find atom on the top surface
                if np.max( surface.positions[ :, 2 ] ) - surface[ at ].position[ 2 ] < 1.0:
                    indx = at
                    break

    center                     = surface[ indx ].position
    surface.positions[ :, 0 ] -= center[ 0 ]
    surface.positions[ :, 1 ] -= center[ 1 ]

    # Set new cell
    surface.cell[ 0 ] = A1_surf
    surface.cell[ 1 ] = A2_surf

    # Remove all atoms outside the new supercell
    surface = remove_atoms_outside( surface, shift = shift, threshold = threshold)

    # Calculate angle between A1_surface and x axis
    A1angle = np.arcsin( surface.cell[ 0, 1 ] / length( surface.cell[ 0 ] ) )
    if surface.cell[ 0, 0 ] / length( surface.cell[ 0 ] ) < 0:
        A1angle = np.pi - A1angle
 
    # Rotate the cell so that A1_surface is parallel to x axis.
    surface.rotate( -A1angle / np.pi * 180., 'z', center = ( 0., 0., 0. ), rotate_cell = True)

    # Verify: Calculate expecte number of atoms in new supercell based on area.
    # If it is different from number of atoms in surface, then we have problem.
    S_primitive  = area( a1, a2 )
    S_surface    = area( A1_surf, A2_surf )
    nat_expected = S_surface / S_primitive * len( sys0 )
    nat          = len( surface )
    if np.abs( nat_expected - nat ) > 0.1:
        print "WARNING: Expected number of atoms is %10i but real number of atoms is %10i" %(nat_expected, nat)
    
    return surface

if __name__ == "__main__":
    from ase.build import fcc111
    from ase.io import write
    primitive_cell = fcc111( 'Au', ( 1, 1, 5 ), vacuum=7.5)  
    matrix = ( 7, 0, 6, 6 )
    system = create_matrix_surface( primitive_cell, matrix  )
    write( 'POSCAR', system, format = 'vasp', direct = True, vasp5= True)
