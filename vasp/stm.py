'''
Has not been cleaned up
'''
"""
***********************************
Tersoff-Hamann approximation of STM

Usage:
from dlescripts.getstm import tersoff_hamann

tersoff_hamann('PARCHG','OUTPUT_PREFIX',TOP,LDOS_IN,ref)

************************************
"""
#!/usr/bin/env python
from ase import *
from ase.io import *
import numpy as np
from ..vasp.chgcar import *
from ..small_tools import file_exist
from scipy.interpolate import interp1d
#from datetime import datetime

def tersoff_hamann(CONTCAR,INDATA,OUTPRE,TOP,LDOS_IN,ref):    
    print ('***********************************' )
    print ('Tersoff-Hamann approximation of STM' )
    print ('' )
    print ('The constant current STM image is  ' )
    print ('simulated for the following:' )
    print (''    )
    print ('Input data: ',INDATA )
    print ('Value of LDOS: ',LDOS_IN,' e/A^3' )
    print ('The program will add ',TOP,' point' )
    print ('in z directions' )
    print ('STM image in CHGCAR format is written' )
    print ('in ',OUTPRE+'.STM.vasp' )
    print ('and in 2D crystal coordinates in: ',OUTPRE+'.STM.vasp' )
    print ('' )
    print ('************************************' )
    print ('' )
    print ('' )

    #Read CONTCAR to define the system.
    print ("checking file: "+CONTCAR )
    file_exists, txt=file_exist(CONTCAR)
    if not file_exists:
        print ( txt )
    else:
        print ( "checking file: "+INDATA )
        file_exists, txt=file_exist(INDATA)
        if not file_exists:
            print ( txt )
    if file_exists:
        system=read(CONTCAR)
        total0 = read_chgcar(INDATA,CONTCAR=CONTCAR)
        total0 = total0.chg[ 0 ]
        ng = total0.shape
        min=np.min(total0)
        max=np.max(total0)
        print ('Min ldos ',min,' e/A^3' )
        print ('Max ldos ',max,' e/A^3' )
        total = np.zeros([ng[0],ng[1],ng[2] + TOP])
        total[:,:,0:ng[2]]=total0[:,:,0:ng[2]]
        total[:,:,ng[2]:ng[2]+TOP]=total0[:,:,0:TOP]

        # deallocate total0 for saving memory
        del total0   
 
        ldos = LDOS_IN  
        print ('Searching for ISO ',LDOS_IN, ' e/A^3' )
        if (ldos < min) or (ldos > max):
            print ('Error: iso value must be in range [',\
                   min,\
                   max,\
                   ']' )
        else:
            fileout=open(OUTPRE+'.STM','w')
            data=np.empty([total.shape[0],total.shape[1],1])
            for xx in range(total.shape[0]):
                for yy in range(total.shape[1]):
                    print ('Scanning :', xx, yy )
                    col_xy = total [xx,yy,:]
                    zz = col_xy.shape[0] -  np.argmax ( col_xy[::-1] >= ldos ) -1
                    if zz == 0:
                        fileout.write( str(xx) + ' ' + str( yy ) +  ' 0. \n' )
                        data[xx,yy,0]=0.
                    else:
                        zarray = [ zz-2, zz-1, zz, zz+1, zz+2, zz+3] 
                        darray = np.take ( col_xy , zarray )
                        f3=interp1d(zarray,darray,kind='cubic')
                        zz1=np.arange(zz,zz+1,0.0001)
                        dd1=f3(zz1)
                        dist = ( dd1 - ldos ) ** 2 
                        min_indx = dist.argmin ()
                        z1 = zz1 [min_indx]
                        fileout.write( str( xx/float(ng[0]) ) + ' ' + \
                                       str( yy/float(ng[1]) ) + ' ' +\
                                       str( z1/float(ng[2])-ref ) + '\n' )
                        data[xx,yy,0]=z1/float(ng[2]) - ref
 
            write_chgcar(OUTPRE+'.STM.vasp',system,data=data)
            fileout.close()


def convert_2D_to_3D(DATA2D,system,DATA3D):
    """
    Coded on 02/22/2013
    """
    file_exists, txt=file_exist(DATA2D)
    if not file_exists:
        exit(txt)
    file_exists, txt=file_exist(DATA3D)
    if file_exists:
        exit('File '+DATA3D+' exist. Use different filename or remove '+DATA3D+'.')
    #del system.constraints
    print ('***********************************' )
    print ('Converting 2D (three column data) STM to' )
    print ('3D (in CHGCAR format) STM.' )
    print ('' )
    print ('2D Data: ',DATA2D )
    print ('3D output: ',DATA3D )
    print ('************************************' )
    print ('')
    print ('')
   
    f=np.loadtxt(DATA2D)
    x=f[:,0]
    y=f[:,1]
    z=f[:,2]
    
    ndata=len(x)
    #Counting NY
    ny=0
    for i in range(ndata):
        if x[i]==0.:
            ny+=1
    nx = int (ndata / ny )

    if ndata !=  (ny * nx):
        exit('Someting went wrong. The estimated values of nx, ny are not correct')    

    data = np.zeros([nx,ny,1])
    item = 0
    for i in range(nx):
        for j in range(ny):
            data[i,j,0] = z[item]
            item += 1
    write_chgcar(DATA3D,system,data=data)


def expand_2D(DATA2D,system,OUT,n, m):
    """
    Coded on 04/10/2015
    """
    file_exists, txt=file_exist(DATA2D)
    if not file_exists:
        file_exists_gz, txt_gz=file_exist(DATA2D+'.gz')
        if not file_exists_gz:
            exit(txt +'\n' + txt_gz)
        else:
            DATA2D +='.gz'
            
    file_exists, txt=file_exist(OUT)
    if file_exists:
        exit('File '+OUT+' exist. Use different filename or remove '+OUT+'.')
    #del system.constraints
    print ('***********************************')
    print ('Expand 2D (three column data) STM ')
    print ('')
    print ('2D Data: ',DATA2D)
    print ('Output: ',OUT)
    print ('Expand: ', n, 'x ', m)
    print ('************************************')
    print ('')
    print ('')
    
    f=np.loadtxt(DATA2D)
    ix=f[:,0]
    iy=f[:,1]
    iz=f[:,2]
    
    ndata=len(ix)
    fout=open(OUT, 'w')
    for i in range(ndata):
        xn = np.array([ix[i],iy[i],iz[i]])
        x  = np.dot (xn, system.cell)
        for j in range(-1,n):
            for k in range(-1, m):
                x_out = x + j*system.cell[0,:] + k*system.cell[1,:]
                fout.write( "%15.10f %15.10f %15.10f \n" %(x_out[0], x_out[1], x_out[2]) )
    fout.close()
        
