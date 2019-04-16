import matplotlib.pyplot as plt
import matplotlib.cm as cmx
import numpy as np
from  dlePy.qe.projwfc import *

def get_eG( band_data, nbnd ):
    ik  = 21
    bands = band_data[ :, 1 ]
    G4 = bands[ ik * nbnd + 3 ]
    G5 = bands[ ik * nbnd + 4 ]

    return G4, G5

def get_eL( band_data, nbnd ):
    ik  = 0
    bands = band_data[ :, 1 ]
    L4 = bands[ ik * nbnd + 3 ]
    L5 = bands[ ik * nbnd + 4 ]

    return L4, L5

def plot( band_data, Ef, nkstotal, legend, minlabel = 'Min', maxlabel = 'Max' ):
    nKPTS = nkstotal 
    nL = 0
    nG = nL + 21  #21 is G point
    nX  = nKPTS - 1
    xmin= 0.
    xmax= 1.
    ymin= -2.
    ymax= 2.

    plt.figure( figsize = ( 4.5, 3.5 ) )

    map = band_data[ :, 2]
    plt.scatter( band_data[ :, 0 ]/ float(nKPTS-1) , band_data[ :, 1] - Ef, c =map, cmap = cmx.jet, s = 10, edgecolors = 'none', label = legend )

    plt.axhline(y=0, linewidth=0.5, color='gray')
    plt.axvline(x=(nG )/float(nKPTS-1), linewidth=0.5, color='gray')

    xticks = [      (nL )/float(nKPTS-1),
                    (nG )/float(nKPTS-1),
                    (nX )/float(nKPTS-1),
             ]

    xticklabels = ['L', r'$\Gamma$', 'X']

    plt.xticks( xticks, xticklabels )
    plt.ylabel('Binding Energy [eV]' )

    plt.xlim(xmin,xmax)
    plt.ylim(ymin,ymax)
    plt.legend(loc=2,fontsize=10)
    plt.colorbar( ticks = [ ] )
    plt.text( 1.07, -2.2, minlabel ) 
    plt.text( 1.07, 2.1, maxlabel ) 
    plt.savefig( legend + '.png', dpi=300, bbox_inches='tight')

if __name__ == "__main__":
    loc = 'STEP_03_PROJWFC/'
    filename = 'output.dat'
    output = loc + filename
    
    # Reading output file of a projwfc run
    data = parse( output )

    # number of kpoints
    nkstotal = get_number_of_kpoinst( data )
 
    # number of band
    nbnd  = get_number_of_bands( data )

    # List all projection   
    list_projections( data ) 

    # Indices of s projection
    sdos_index = [ 0, 4 ]    

    # Sum all projection with indices in sdos_index
    sdos       = sum_states( data, state_index = sdos_index ) 

    # Define Fermi energy
    G3, G4 = get_eG( sdos, nbnd )
    L3, L4 = get_eL( sdos, nbnd )
    Ef = ( L4 + G3 ) / 2.  #Ef can also be supplied manually

    plot( sdos, Ef, nkstotal, legend = 's-projection' )

    # Indices of p projection
    pdos_index = [ 1, 2, 3, 5, 6, 7 ]

    # Sum all projection with indices in pdos_index
    pdos       = sum_states( data, state_index = pdos_index )

    plot( pdos, Ef, nkstotal, legend = 'p-projection' )

    # Calculate sdos /( sdos + pdos )
    ratio = sdos[ :, -1 ] / ( sdos[ :, -1 ] + pdos[ :, -1 ] + 1.e-10 )  #1.e-10 to avoid divide by 0

    # Construct data for plotting
    ratio_data = sdos.copy( )
    ratio_data[ :, -1 ] = ratio

    plot( ratio_data, Ef, nkstotal, legend = 's_p-decomposition', minlabel = 'p', maxlabel = 's' )
