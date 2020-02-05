import os

def get_reduce_atom_list ( atoms ):
    mol =  atoms.get_chemical_symbols( )
    if len ( mol ) > 1:
        for i in range (len( mol ) - 2 , -1 , -1):
            if mol [ i ] == mol [ i + 1 ]:
               del mol [ i ]
    return mol

def gen_POTCAR( system, potcar_loc, map = {}  ):
    atom_list= get_reduce_atom_list( system )
    map_in = map
    map  = { }
    for at in atom_list:
        map[ at ] = at
    for at in map_in.keys():
        map[ at ] = map_in[ at ]
    print ( map )
    for at in range( len( atom_list ) ):
        if at == 0:
            os.system( 'cat ' + potcar_loc + 'POTCAR.' +  map[ atom_list[ at ] ] + ' > POTCAR' )
        else:
            os.system( 'cat ' + potcar_loc + 'POTCAR.' +  map[ atom_list[ at ] ] + ' >> POTCAR' )

