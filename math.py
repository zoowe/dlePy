from   ase.io import read
import numpy as np
import os


def dotproduct( v1, v2 ):
    return sum( ( a * b ) for a, b in zip( v1, v2 ) )

def length(v):
    return np.sqrt(dotproduct(v, v))

def area( v1, v2):
    '''
    calculate area created by two vectors
    '''
    return length( np.cross( v1, v2 ) )

