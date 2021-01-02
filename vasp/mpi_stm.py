"""
Has not been cleaned up 
"""
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
from mpi4py import MPI
import sys
import numpy as np
from datetime import datetime

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
name = MPI.Get_processor_name()

def print_mpi(rank,txt,endline=True):
    if rank == 0:
        if endline==True:
            print ( txt )
        else:
            print ( txt ),

def print_time(rank, time):
    print_mpi (rank,'Total time: '+str(time))

def tersoff_hamann(CONTCAR,INDATA,OUTPRE,TOP,LDOS_IN,ref):    
    from scipy.interpolate import interp1d
   
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()
    name = MPI.Get_processor_name()
    print_mpi(rank,'***********************************')
    print_mpi(rank,'Tersoff-Hamann approximation of STM')
    print_mpi(rank,'Run with '+str(size)+' cpus on'+str( name))
    print_mpi(rank,'The constant current STM image is  ')
    print_mpi(rank,'simulated for the following:')
    print_mpi(rank,''    )
    print_mpi(rank,'Input data: '+INDATA)
    print_mpi(rank,'Value of LDOS: '+str(LDOS_IN)+' e/A^3')
    print_mpi(rank,'The program will add '+str(TOP)+' point')
    print_mpi(rank,'in z directions')
    print_mpi(rank,'STM image in CHGCAR format is written')
    print_mpi(rank,'in '+OUTPRE+'.STM.vasp')
    print_mpi(rank,'and in 2D crystal coordinates in: '+OUTPRE+'.STM')
    print_mpi(rank,'')
    print_mpi(rank,'************************************')
    print_mpi(rank,'')
    print_mpi(rank,'')
    #Read CONTCAR to define the system.
    print_mpi(rank,"checking file: "+CONTCAR)
    file_exists, txt=file_exist(CONTCAR)
    if not file_exists:
        print_mpi(rank,txt)
    else:
        print_mpi(rank,"checking file: "+INDATA)
        file_exists, txt=file_exist(INDATA)
        if not file_exists:
            print_mpi(rank,txt)
#    file_exists = comm.bcast(file_exists,root=0)
    if file_exists:
        system=read(CONTCAR)
        #del system.constraints
    
    
    #We don't need this, thus we need to count how many lines
    #we need to ignore
        nline=9
        startline=nline+system.get_number_of_atoms()
   
    # Reading INDATA
        print_mpi ( rank, 'Reading '+INDATA+' ...')
        t1=datetime.now() 
    #Open file LOCPOT
        file=open(INDATA,'r')
    
    #Read all ignored lines
        for i in range(startline):
            dump=file.readline()
    
    
        ngr=file.readline().split()
    #Make it becomes integers 
        ng = (int(ngr[0]), int(ngr[1]), int(ngr[2]))
        print_mpi(rank,'Grid info: '+str(ng[0]) + ' '+str(ng[1])+' '+str(ng[2]))
        total0 = np.empty(ng)
        for zz in range(total0.shape[2]):
            for yy in range(total0.shape[1]):
                total0[:, yy, zz] = np.fromfile(file, count = total0.shape[0],
                                           sep=' ')
    
        file.close()
        min=np.min(total0)
        max=np.max(total0)
        print_mpi(rank, 'Min ldos '+str(min/system.get_volume())+' e/A^3')
        print_mpi(rank, 'Max ldos '+str(max/system.get_volume())+' e/A^3')
        total = np.zeros([ng[0],ng[1],2*ng[2]])
        total[:,:,0:ng[2]]=total0[:,:,0:ng[2]]
        total[:,:,ng[2]:2*ng[2]]=total0[:,:,0:ng[2]]
        del total0
        #Done reading data
        t2=datetime.now()
        print_time (rank, t2-t1)   

        VALUE=LDOS_IN
        ldos=VALUE*system.get_volume()
        print_mpi(rank, 'Searching for ISO '+str(VALUE)+' e/A^3')
        if (ldos < min) or (ldos > max):
            print_mpi(rank, 'Error: iso value must be in range ['+str(min/system.get_volume())+','+str(max/system.get_volume())+']')
        else:
            fileout=open('data.'+str(rank),'w')

            ### Parallel this part
            #Arranging point for each cpu
            print_mpi(rank, 'Arranging data ...') 
            t1=datetime.now()

            npts=total.shape[0]*total.shape[1]  #Total point
            ppc=int(npts/size)                  #point per cpu
            pre=npts-ppc*size                   #remaining point
            rppc=np.zeros(size,dtype=int)       #real ppc 
            rppc[:]=ppc
            if pre > 0:
                for i in range(pre):
                    rppc[i]+=1
            for i in range(size):
                print_mpi (rank, 'CPU '+str(i)+' calculates '+str(rppc[i])+' points')
            data_point=np.zeros([npts,2],dtype=int)
            data=np.zeros([npts])
            for xx in range(total.shape[0]):
                for yy in range(total.shape[1]):
                    i=xx*total.shape[1]+yy
                    data_point[i,0]=xx
                    data_point[i,1]=yy
            #Assign index for each cpu
            start_index=np.zeros(size,dtype=int)
            end_index=np.zeros(size,dtype=int)
            start_index[0]=0
            end_index[0]=rppc[0]-1

            for i in range(1,size):
                for j in range(i):
                    start_index[i]+=rppc[j]
                end_index[i]=start_index[i]+rppc[i]-1   
            for i in range(size):
                print_mpi(rank, str(i)+' '+str(start_index[i])+' '+str(end_index[i])) 
            t2=datetime.now()
            print_time (rank, t2-t1)

            print_mpi (rank, 'SCANNING...')
            t1 = datetime.now()

            for i in range(start_index[rank],end_index[rank]+1):
                xx=data_point[i,0]
                yy=data_point[i,1]
                for zz in range(ng[2]-1+TOP,-1,-1):
                    if zz == 0:
                        print >>fileout,xx,yy,'0.'
                        data[i]=0.-ref
                    else:
                        if (total[xx,yy,zz]<=ldos) and (total[xx,yy,zz-1]>=ldos):
                            zarray=[zz-3,zz-2,zz-1,zz,zz+1,zz+2]
                            darray=np.take( total [ xx,yy, :] , zarray ) 
                            f3=interp1d(zarray,darray,kind='cubic')
                            zz1=np.arange(zz-1,zz,0.0001)
                            dd1=f3(zz1)
                            dist = ( dd1 - ldos ) ** 2 
                            min_indx = dist.argmin ()
                            z1 = zz1 [min_indx]
  
                            print>>fileout,xx,yy,z1
                            data[i]=z1/float(total.shape[2]/2.)-ref
                            break

            t2 = datetime.now()
            fileout.close()
            print 'CPU '+str(rank)+' TIME '+str(t2-t1)
            comm.Barrier()
            # Broadscating data
            print_mpi(rank, 'Broadscating data...',endline=False)
            t1=datetime.now()
            tmp_data = np.zeros([npts])
            if rank == 0:
                for i in range(1, size):
                    comm.Recv(tmp_data, source= MPI.ANY_SOURCE)
                    data += tmp_data
            else:
                comm.Send(data,dest=0)
            comm.Barrier()
            t2=datetime.now()
            print_time(rank, str(t2-t1))
            
            if rank == 0:
                stm=np.zeros([total.shape[0],total.shape[1],1])
                for xx in range(total.shape[0]):
                    for yy in range(total.shape[1]):
                        i = xx*total.shape[1] + yy
                        stm[xx,yy,0]=data[i]

                write_chgcar(OUTPRE+'.STM.vasp',system,data=stm)

            comm.Barrier()

