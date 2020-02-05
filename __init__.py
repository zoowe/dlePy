import sys

str='''
dlePy V 0.2         
Updated: 08/23/2018 

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

"""
Cleaning log:
06/06/2018: strucmod.py and supercell.py 

Update log:
Disregard all log before 06/06/2018.
"""
