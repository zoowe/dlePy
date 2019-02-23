'''
Has not been cleaned up
'''
from numpy import *
import sys
from ..small_tools import file_exist 
import gzip as gz

def get_band_atom(ATOM,NATOMS,NKPTS,NBANDS,PROCAR,NKPTS0,NSPIN,hse=False,nskip=0,LSORBIT=False):
   """
   To read data for a ATOM and put them in files
   ATOM.SPIN.XXX.
   
   ATOM: Atom index, start from 1 to NATOMS
   NATOMS: number of atoms
   NKPTS: number of k-points in the concerned PROCAR file
   NBANDS: number of bands
   PROCAR: procar file name, PROCAR.X
   NKPTS0: number of kpoints from PROCAR.0 to PROCAR.X-1
   NSPIN: 0 or 1, up or down
   """
   import linecache

   file_exists, txt = file_exist(PROCAR)
   if not file_exists:
        PROCAR+='.gz'
        file_exists, txt = file_exist(PROCAR)
        #print file_exists, txt
        f = gz.open( PROCAR )
        if not file_exists:
            print 'STOP'    
            exit(txt)
   else:
      f = open( PROCAR )

   lines = f.readlines()
   f.close()

   if ATOM  < 10:
      subfix='00'+str(ATOM )
   elif ATOM < 100:
      subfix='0'+str(ATOM )
   else:
      subfix=str(ATOM )
   #print subfix
   dos=[0.]*10
   dataline=0
   for ns in range(NSPIN):
      if ATOM > NATOMS:
          file=open('TOTAL','a')
      else:
          file=open('ATOM.'+str(ns)+'.'+subfix,'a')
      for ik in range(nskip+1,NKPTS+1):
          if NATOMS == 1:
              if ns == 0:
                 kline=4+(ik-1)*((NATOMS+4)*(NBANDS)+3)
              if ns == 1:
                 kline=4+(NKPTS-1)*((NATOMS+4)*(NBANDS)+3)+NBANDS*(NATOMS+4)
                 kline+=4+(ik-1)*((NATOMS+4)*(NBANDS)+3)
              if LSORBIT:
                 kline=4+(ik-1)*((NATOMS+5+(NATOMS+1)*3)*(NBANDS)+3) 
          else:
              if ns == 0:
                 kline=4+(ik-1)*((NATOMS+5)*(NBANDS)+3)
              if ns == 1:
                 kline=4+(NKPTS-1)*((NATOMS+5)*(NBANDS)+3)+NBANDS*(NATOMS+5)
                 kline+=4+(ik-1)*((NATOMS+5)*(NBANDS)+3)
              if LSORBIT: 
                 kline=4+(ik-1)*((NATOMS+5+(NATOMS+1)*3)*(NBANDS)+3) 
#          print kline
          #line= linecache.getline(PROCAR,kline)
          line= lines[ kline ]
          print line
          lenline=len(line.split())
          weight=float(line.split()[lenline-1])
          for iband in range(1,NBANDS+1):
              if NATOMS == 1:
                  eline= kline + (iband-1)*(NATOMS+3) + iband + 1
                  if LSORBIT: 
                      eline= kline + (iband-1)*(NATOMS+3+(NATOMS+1)*3) + iband + 1
              else:
                  eline= kline + (iband-1)*(NATOMS+4) + iband + 1
                  if LSORBIT:
                      eline= kline + (iband-1)*(NATOMS+4+(NATOMS+1)*3) + iband + 1 
#              print eline
              #line= linecache.getline(PROCAR,eline)
              line= lines[ eline - 1]
              ener=float(line.split()[4])
              dataline=eline+2+ATOM
              #line= linecache.getline(PROCAR,dataline)
              line= lines[ dataline - 1]
              for i in range(1,11):
                  dos[i-1]=float(line.split()[i]) 
              print >> file, '%5i %6i %10.6f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f' %(ik-1+NKPTS0-nskip,iband, ener,dos[0],dos[1],dos[2],dos[3],dos[4],dos[5],dos[6],dos[7],dos[8],dos[9],weight)
      file.close()
      del lines

def get_band(ATLIST,NATOMS,NKPTS,NBANDS,c,FILENAME,efermi,SPIN,REF,THRESHOLD,hse=False,nskip=0):
    """    
    Sum projected dos for atoms
    
    ATLIST: List of atoms you want to sum.
    NATOMS: number of atoms
    NKPTS: total number of k-points
    c: column you want to sum. c=1 >> sum column 3, s states.
                               c=10, get the total
    FILENAME: The filename of output 
    efermi: Fermi energy, for reference it will be set to Zero.
    SPIN: 0 or 1 == spin up or down
    REF: File name of the reference data. 
         Usually, it is the sum of projected dos for all atoms.
         The program will print out % contribution with respect to REF
         'NA' to ignore this option.
    THRESHOLD: If the sum is smaller than THRESHOLD, 
               the data will not be printed out.
    """
    from numpy import loadtxt
    
    ATOM=ATLIST[0]
    if ATOM < 10:
       subfix='00'+str(ATOM)
    elif ATOM < 100:
       subfix='0'+str(ATOM)
    else: 
       subfix=str(ATOM)    

    # Read a template
    file_exists, txt = file_exist('ATOM.'+str(SPIN)+'.'+subfix)
    if not file_exists:
        print 'STOP'
        exit(txt)

    f1=loadtxt('ATOM.'+str(SPIN)+'.'+subfix)
    k=f1[:,0]
    b=f1[:,1]
    e=f1[:,2]
    w=f1[:,13] 
    dos=f1[:,c+2]
    dos=0.
    
    for ATOM in ATLIST:
       if ATOM < 10: 
          subfix='00'+str(ATOM) 
       elif ATOM < 100: 
          subfix='0'+str(ATOM) 
       else:  
          subfix=str(ATOM)    

       file_exists, txt = file_exist('ATOM.'+str(SPIN)+'.'+subfix)
       if not file_exists:
           print 'STOP'
           exit(txt)

       print 'ATOM.'+str(SPIN)+'.'+subfix 
       f1=loadtxt('ATOM.'+str(SPIN)+'.'+subfix) 
       dostmp=f1[:,c+2]
       dos+=dostmp
    file=open(FILENAME,'w')
    print 'The columns in ',FILENAME,' are:'
    if ( REF == 'NA' ):
       print '%10s %10s %7s %6s %7s' %('k#/NPTS','e-ef','pdos','band#','weigh')
    if ( REF != 'NA' ): 
       print '%10s %10s %7s %14s %6s %7s' %('k#/NPTS','e-ef','pdos','pdos/pdosref','band#','weigh')
 
    if ( REF != 'NA' ):

       file_exists, txt = file_exist(REF)
       if not file_exists:
            print 'STOP'
            exit(txt)

       reff=loadtxt(REF)
       ref=reff[:,2]
    for i in range(len(k)):
       if ( REF == 'NA' ):
          if NKPTS > 1:
              print >>file,'%10.6f %10.6f %7.4f %6i %7.4f' %(k[i]/float(NKPTS-1),e[i]-efermi,dos[i],b[i],w[i])
          else:
              print >>file,'%10s %10.6f %7.4f %6i %7.4f' %('0.',e[i]-efermi,dos[i],b[i],w[i])
#       if ( REF != 'NA' ) and (dos[i]/ref[i] > THRESHOLD) :
       if ( REF != 'NA' ) and ( dos[i] > THRESHOLD ):
          if NKPTS > 1:
              print >>file,'%10.6f %10.6f %7.4f %7.4f %6i %7.4f' %(k[i]/float(NKPTS-1),e[i]-efermi,dos[i],dos[i]/ref[i],b[i],w[i])
          else:
              print >>file,'%10s %10.6f %7.4f %7.4f %6i %7.4f' %('0.',e[i]-efermi,dos[i],dos[i]/ref[i],b[i],w[i])
  
############################# 

def split_bands(NPROCAR,NSPIN,hse=False,nskip=0,LSORBIT=False):
    """
    Split PROCAR
    
    NPROCAR: number of PROCAR files, they have to be named 
             PROCAR.0, PROCAR.1, .. PROCAR.NPROCAR-1
    NSPIN: spin, 1 for non-spinpolarized, 2 for spin-polarized
    """
    import linecache
    from sys import exit
    
    for i in range(1000):
       if i < 10:
          subfix='00'+str(i)
       elif i  < 100:
          subfix='0'+str(i)
       else:
          subfix=str(i)

       file_exists, txt = file_exist('ATOM.0.'+subfix)
       if file_exists:
           print 'STOP'
           exit('File ATOM.0.'+subfix+' is already exist. You should delete all ATOM.X.XXX files before continue.')
       file_exists, txt = file_exist('ATOM.1.'+subfix)
       if file_exists:
           print 'STOP'
           exit('File ATOM.1.'+subfix+' is already exist. You should delete all ATOM.X.XXX files before continue.')

       file_exists, txt = file_exist('TOTAL')
       if file_exists:
           print 'STOP'
           exit('File TOTAL is already exist. You should delete TOTAL file before continue.')

    print 'The columns in ATOM.X.XXX files are:'
    print '%5s %6s %10s %7s %7s %7s %7s %7s %7s %7s %7s %7s %7s %7s' %('k#','band#','e','s','py','pz','px','dxy','dyz','dz2','dxz','dx2','tot','w')

    NKPTS0=0
    for pro in range(NPROCAR):
        PROCAR='PROCAR.'+str(pro)
        print 'Parsing file ',PROCAR
        file_exists, txt = file_exist(PROCAR)
        if not file_exists:
            PROCAR='PROCAR.'+str(pro)+'.gz'
            file_exists, txt = file_exist(PROCAR)
            #print file_exists, txt
            if not file_exists:
                print 'STOP'     
                exit(txt)
            else:
                f = gz.open( PROCAR )
        else: 
            f = open( PROCAR )
 
        lines = f.readlines( )
        f.close()
#        line=linecache.getline(PROCAR,2)
        line=lines[ 1 ]
#        print 'x',line
#        NATOMS=int(line.split()[11])
        NATOMS=get_NATOMS()
#        line=linecache.getline(PROCAR,2)
        NKPTS=int(line.split()[3])
#        line=linecache.getline(PROCAR,2)
#        NBANDS=int(line.split()[7])
        NBANDS=get_NBANDS()
    
        print NATOMS, NKPTS, NBANDS,NSPIN
    
        for ATOM in range(1,NATOMS+2):
            print 'COLLECTING DATA FOR ATOM ', ATOM
            get_band_atom(ATOM,NATOMS,NKPTS,NBANDS,PROCAR,NKPTS0,NSPIN,hse=hse,nskip=nskip,LSORBIT=LSORBIT)
        NKPTS0+=NKPTS-nskip
        #linecache.clearcache()

def get_number_of_kpoints(NPROCAR,hse=False, nskip=0):
    """
    To count number of k-points from all provided PROCAR files
   
    NPROCAR: number of PROCAR files, they have to be named 
             PROCAR.0, PROCAR.1, .. PROCAR.NPROCAR-1
    """
#    import linecache
    import os
    from numpy import loadtxt, sum

    NKPTS0=0
    os.system('rm -f .tmp')
    for pro in range(NPROCAR):
        PROCAR='PROCAR.'+str(pro)
        print PROCAR
        file_exists, txt = file_exist(PROCAR)
        if not file_exists:
            PROCAR='PROCAR.'+str(pro)+'.gz'
            file_exists, txt = file_exist(PROCAR)
            #print file_exists, txt
            if not file_exists:
                print 'STOP'    
                exit(txt)
            else:
                f = gz.open( PROCAR )
        else:
            f = gz.open( PROCAR )

        lines = f.readlines( )
        f.close()
        #os.system('head -2 '+PROCAR+' | tail -1 | awk \'{print $4}\' >> .tmp') i
        os.system( 'echo ' + lines[ 1 ] + ' >> .tmp' )
        del lines
    os.system('echo 0 >> .tmp') 
    f=loadtxt('.tmp')
    nk=f[:]
    NKPTS0=int(sum(nk))-NPROCAR*nskip
    return NKPTS0

def get_NATOMS():
    """
    To get number of ATOMS
    
    Need only PROCAR.0
    """
    import linecache
 
    PROCAR='PROCAR.0'
    file_exists, txt = file_exist(PROCAR)
    if not file_exists:
        PROCAR+='.gz'
        file_exists, txt = file_exist(PROCAR)
        #print file_exists, txt
        f = gz.open( PROCAR, 'r' )
        if not file_exists:
            print 'STOP'    
            exit(txt)
    else:
        f = open( PROCAR, 'r' )

    line = f.readline()
    line = f.readline()
    #line=linecache.getline(PROCAR,2)
    
    LASTITEM=line.split()[len(line.split())-1]
    NATOMS=int(line.split()[len(line.split())-1].replace('ions:','').strip())
 
    f.close()
    return NATOMS
 
def get_NBANDS():
    """
    To get number of NBANDS 
    
    Need only PROCAR.0
    """
    import linecache

    PROCAR='PROCAR.0'
    file_exists, txt = file_exist(PROCAR)
    if not file_exists:
        PROCAR+='.gz'
        file_exists, txt = file_exist(PROCAR)
        #print file_exists, txt
        if not file_exists:
            print 'STOP'    
            exit(txt)
        else:
            f = gz.open( PROCAR )
    else:
        f = gz.open( PROCAR )
    lines = f.readlines()
    f.close()
    #line=linecache.getline(PROCAR,2)
    line= lines[ 1 ]
    for i in range(len(line.split())):
        if 'bands' in line.split()[i]:
            item=line.split()[i].replace('bands:','').strip()
            if item == '':
                NBANDS=int(line.split()[i+1].strip())
            else:
                NBANDS=int(item)
            break
#    NBANDS=int(line.split()[7])

    #linecache.clearcache()
    del lines 
    return NBANDS

def get_bands_structure(NPROCAR,NSPIN,DATA,LSORBIT=False,hse=False, nskip=0):
    """
    Get band structure for  PROCAR files
    
    NPROCAR: number of PROCAR files, they have to be named 
             PROCAR.0, PROCAR.1, .. PROCAR.NPROCAR-1
    NSPIN: spin, 1 for non-spinpolarized, 2 for spin-polarized
    LSORBIT: True or false. If true, spin orbit coupling
    """
    import linecache
    from sys import exit
    import os
    from numpy import loadtxt
 
    TOTALKNPTS= get_number_of_kpoints(NPROCAR,hse=hse,nskip=nskip)

    fileout=open(DATA,'w')
    NKPTS0=0 
    for pro in range(NPROCAR):
        PROCAR='PROCAR.'+str(pro)
        print 'Parsing file ',PROCAR
        file_exists, txt = file_exist(PROCAR)
        if not file_exists:
            PROCAR='PROCAR.'+str(pro)+'.gz'
            file_exists, txt = file_exist(PROCAR)
            #print file_exists, txt
            if not file_exists:
                print 'STOP'    
                exit(txt)


        line=linecache.getline(PROCAR,2)
        NATOMS=get_NATOMS()
        NKPTS=int(line.split()[3])
        NBANDS=get_NBANDS()

        print NATOMS, NKPTS, NBANDS,NSPIN

        os.system('grep "band " '+PROCAR+' | awk \'{print $2,$5}\'  > .DATA')
        f=loadtxt('.DATA')
        e=f[:,1]
        for i in range(1,NKPTS+1):
            knum=i
            for j in range(NBANDS):
                ener=e[(knum-1)*NBANDS + j]
                print >>fileout, float(knum-1+NKPTS0)/float(TOTALKNPTS-1),ener
        NKPTS0+=NKPTS

