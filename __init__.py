import sys

str='''
dlePy V 0.2         
Updated: 01/02/2021 

Author: Duy Le
Department of Physics
University of Central Florida

Email: duy.le@ucf.edu
Website: http://www.physics.ufc.edu/~dle
----------------------------------------
'''

print ( str )


if sys.version_info[ 0 ] == 2:
    raise ImportError( 'dlePy requires Python3. This is Python2.' )

