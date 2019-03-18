from   ase.io import read
import numpy as np
import os

dotproduct    = lambda v1, v2    : np.dot (v1, v2 ) 
length        = lambda v         : np.linalg.norm( v )
area          = lambda v1, v2    : length( np.cross(v1, v2 ) )
str2bool      = lambda s         : s.lower() in ( "yes", "true", "t", "1" )
tripleproduct = lambda v1, v2, v3: dotproduct( v3, np.cross(v1, v2 ) )
