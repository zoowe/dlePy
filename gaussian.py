import numpy as np

def f_gauss(x,muy,sigma):
    result= 1. / ( sigma * np.sqrt ( 2 * np.pi ) ) 
    result*= np.exp (-1./2. * ( (x-muy) / sigma ) **2 )
    return result

def gauss_smooth(muydata, amp, sigma, d, resolution,scale=1., xmin=-100, xmax=100, verbose = True):
    width= np.max( muydata ) - np.min( muydata )
    if xmin == -100:
        xmin = np.min( muydata ) - d / 100. * width 
    if xmax == 100:
        xmax = np.max( muydata ) + d / 100. * width 

    if verbose:
        print ( 'Data width=', width, 'Xmin =',xmin,'Xmax =',xmax )

    npts = int( ( xmax - xmin ) / resolution ) + 1

    xgrid = np.zeros( [ npts ] )

    xgrid = np.array( range( npts ) ) * resolution
    xgrid[ : ] += xmin
    ampgrid  = np.zeros( [ npts ] )

    for i in range( npts ):
        ampgrid[ i ] = np.sum( amp[ : ] * f_gauss( xgrid[ i ], muydata[ : ], sigma ) * scale )

    return xgrid, ampgrid 
