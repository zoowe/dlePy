import numpy as np
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms
from phonopy.interface.vasp import write_supercells_with_displacements 
from ase.io import read, write
from ase.visualize import view
import os
import gzip as gz

def gen_POSCAR( system, top_dir = "VIB", last_run_top_dir = "", disp_indices = [ ], disp   = 0.01, dim = ( 1, 1, 1 ), script_only = False, more="" ):
    """
    atoms: ase Atoms object of relaxed system.
    top_dir: path to directory where POSCARs will be saved to.
    last_run_top_dir: top_dir of a previous run. Empty string if this is the first run.
    disp_indices: list of indices of atoms allowed to be displaced.
    disp: In angstrom, displacement
    dim: dimention of supercell. (1, 1, 1) for vib freq.
    """

    dim_str = str( dim[ 0 ] ) + ' ' + str( dim[ 1 ] ) + ' '+ str( dim[ 2 ] ) 

    if len( disp_indices ) == 0 or len( disp_indices ) > len( system ):
        print ("""
        This module is for doing vibritional calculations with fixing a part of the system.
        The provided disp_indices indicates that no atom is allowed to displace.
        Condition: 0 < len( disp_indices ) <= len( system ) ( 0 < %i <= %i) does not sastify.
        
        """ %( len( disp_indices ) ,  len( system ) ) )
        exit( 'EXITTING' )

    if os.path.isdir( top_dir ):
        if not script_only:
            print ("""
            %s exists. Change `top_dir` or remove %s.
            """ %( top_dir, top_dir ) )
            exit( 'EXITTING' )

    vasp_input = [ 'INCAR', 'POTCAR', 'KPOINTS', 'job' ]
    for file in vasp_input:
        if not os.path.isfile( file ):
            print ("""
            %s for vasp run is needed at current directory.
            """ % ( file ) )
            if not script_only: 
                exit( 'EXITTING' )

    #Generate displacemnt for displaceble atoms
    disp_system = system[ disp_indices ]
    os.system( 'mkdir -p ' + top_dir )
    os.chdir ( top_dir )
    if not script_only:
        write( 'POSCAR', disp_system, format = 'vasp', direct = True, vasp5 = True )
        os.system( 'phonopy -d --dim="' + dim_str  + '" --amplitude=' + str( disp ) + " " + more)

    # Create POSCARs with all atoms for VASPRUN
    list_POSCAR = get_list_POSCAR( './' )

    # Get fixed atoms, the ones that were not displaced in disp_system
    fixed_indices = [ at.index for at in system if at.index not in disp_indices ]
    if len( fixed_indices ) > 0:
        fixed_system = system[ fixed_indices ]

    # Create POSCARs for vasp run
    os.system( 'mkdir -p VASPRUNS' )
    if not script_only:
        for item in list_POSCAR:
            disp_system = read( item, format = 'vasp' )
            vasprun_system = system.copy( ) # Original system
            vasprun_system.positions[ disp_indices ] = disp_system.positions
            write( 'VASPRUNS/' + item, vasprun_system, format = 'vasp', direct = True, vasp5 = True )

    os.chdir ( '../' )
   
    if not script_only: 
        os.system( 'cp INCAR KPOINTS POTCAR job ' + top_dir + '/VASPRUNS/' )
    
    gen_analyze_script( top_dir, last_run_top_dir )
    gen_check_script( top_dir )
    gen_collect_vasprun_script( top_dir, fixed_indices, disp_indices  )
    gen_get_BORN_script( top_dir, fixed_indices, disp_indices )

    print ("""
    ###########################################
    All POSCAR files for vasp runs are in folder %s
    cd %s 
    Run analyze.py to for further instruction.
    """ %( top_dir + '/VASPRUNS/', top_dir + '/VASPRUNS/' ) )
    check_atomic_ordering( dir = top_dir )

def check_atomic_ordering( dir = './' ):
    """
    To verify the consistency of POSCAR and SPOSCAR
    """
    sys_POSCAR = read( dir + '/POSCAR' )
    sys_SPOSCAR = read( dir + '/SPOSCAR' )
    dx = sys_SPOSCAR.positions - sys_POSCAR.positions
    if ( np.max( np.abs( dx ) ) > 1.e-5 ):
        print ( 'WARNING: ATOMS MAY HAVE BEEN SORTED, MAKE SURE IT CAUSES NO ERROR' )
 
def get_list_POSCAR( dir ):
    list_POSCAR = []
    for item in os.listdir( dir ):
        if ( 'POSCAR' in item ) and ( len( item ) > 6  ):
            list_POSCAR.append( item )
    return list_POSCAR

def analyze( last_run_top_dir = 'VIB' ): 
    current_run_POSCAR = get_list_POSCAR( './' )
    print ( "There are %i POSCARs that need to be analyze" %( len( current_run_POSCAR ) ) )
    if last_run_top_dir.strip() != "": 
        last_run_dir = '../../' + last_run_top_dir + '/VASPRUNS/'
        last_run_POSCAR = get_list_POSCAR( last_run_dir ) 
        last_run_POSCAR_ = get_list_POSCAR( last_run_dir ) 
        print ("There are %i POSCARs in last run" %( len( last_run_POSCAR ) ) )

        print ('Read last run structure...', end = "" )
        last_frame = [ ]
        last_name  = [ ]
        for last_item in sorted( last_run_POSCAR ):
            last = read( last_run_dir + last_item, format ='vasp' )
            last_frame.append( last )
            last_name.append( last_item )
        print ( 'DONE' )

        print ( 'Read current run structure...', end = "" )
        current_frame = [ ]
        current_name = [ ]
        for current_item in sorted( current_run_POSCAR):
            current = read( './' + current_item, format ='vasp' )
            current_frame.append( current )
            current_name.append( current_item )
        print ( 'DONE' )

        
        map_done = { }
        for i, current in enumerate( current_frame ):
            print ( "%i/%i: Analyze %s" %( i + 1, len( current_run_POSCAR  ), current_name[ i ] ), end = "" )
            diff_list = [ ]
            for j, last in enumerate( last_frame ):
                d = current.positions - last.positions
                diff_list.append( np.sum( np.abs( d ) ) )
                if np.sum( np.abs( d ) ) < 0.000001:
                    map_done[ current_name[ i ] ] = last_name[ j ] 
                    print ("=> already DONE, mapped to %s of the last run" %( last_name[ j ] ), end = "" )
                    #del last_run_POSCAR_[ iitem_last ]
                    break
                else:
                    if j == len( last_frame ) - 1:
                        print ("=> not done", end = "")
                        print (" Min diff ", np.min( np.array( diff_list ) ), end = "" )
            print ("!")

        print ("================================")
        print ("The following POSCARs have been done previously in:")
        print (last_run_dir)
        print ("They will be mapped to this folder")

    if last_run_top_dir.strip() == "":
        print ("This is the first run, so all POSCARs in this folder needs to run")
        map_done = { }

    done_list = [ ]
    for item_current in map_done:
        print ( item_current, map_done[ item_current ] )
        done_list.append( item_current )

    print ("================================")
    print ("The following POSCARs have not been calculated")
    not_done_list = [ ]
    for item in current_run_POSCAR:
        if item not in done_list:
            not_done_list.append( item )

    return map_done, not_done_list

def link_last_run( map_done, last_run_top_dir = 'VIB' ):
    last_run_dir = '../../' + last_run_top_dir + '/VASPRUNS/'
    print ("Linking already-done calculations to this folder")
    for item_current in map_done:
        id_current = item_current.replace( 'POSCAR-', '' ).replace( 'SPOSCAR', 'S' )
        id_last    = map_done[ item_current ].replace( 'POSCAR-', '' ).replace( 'SPOSCAR', 'S' )
        folder_current = 'disp-' + id_current
        folder_last    = 'disp-' + id_last
        os.system( 'ln -sv ' + last_run_dir +  folder_last + ' ' + folder_current )

    return
 
def submit( list ):
    for item in list:
        folder = item.replace( 'POSCAR-', 'disp-' ).replace( 'SPOSCAR', 'disp-S' )
        os.system( 'mkdir -p ' + folder )
        os.chdir ( folder )
        os.system( 'ln -s ../INCAR' )
        os.system( 'ln -s ../POTCAR' )
        os.system( 'ln -s ../KPOINTS' )
        os.system( 'ln -s ../' + item + ' POSCAR' )
        os.system( 'cat ../job | sed \'s/$sys/'+ item.replace( 'POSCAR','' ) +'/\' > job' )
        os.system( 'qsub job' )
        os.chdir ( '../' )

def gen_analyze_script( top_dir, last_run_top_dir = '' ):

    s ="""from dlePy.phonopy.generate import analyze, link_last_run, gen_submit_script

last_run_top_dir = '%s'  # if this is the first run, use empty string 

map_done, not_done_list = analyze( last_run_top_dir )

#This is for mapping the one that has already run to this folder"
if len( map_done ) > 0:
    link_last_run( map_done, last_run_top_dir = last_run_top_dir )

#This is for generating submit.py
if len( not_done_list ) > 0 :
    gen_submit_script( not_done_list ) 
    
    print ("###########################################")
    print ("Run python submit.py to submit nessesary job.")
    print ("Make sure to have job file before running submit.py")
    print ("Once all jobs are done, run collect_vasprun.py (in parent folder of VASPRUNS) to collect and convert vasprun.xml")
 
""" %( last_run_top_dir )

    f = open( top_dir + '/VASPRUNS/analyze.py', 'w' )
    f.write(s)
    f.close()
    return 

def gen_check_script( top_dir ):
    s = """
#Check if any disp- folder has been completed 
import os
from dlePy.jobmon.monitor import if_vasp_done
from dlePy.phonopy.generate import gen_submit_script
from submit import list

fout = open( 'CHECK.log', 'w' )

fout.write( "%-20s %6s %5s %20s \\n" %( 'FOLDER', 'DONE', 'NITER', 'TIME_WRITE') )

not_done_list = []
for item in sorted( list ):
    folder = item.replace( 'POSCAR-', 'disp-' ).replace( 'SPOSCAR', 'disp-S' )
    data = if_vasp_done( folder + '/OUTCAR' )
    fout.write( "%-20s %6s %5s %20s \\n" %( folder, data[ 0 ], data[ 1 ], data[ 3 ]  ) )
    if ( not data[ 0 ] ):
        not_done_list.append( item )
fout.close()

if len( not_done_list ) > 0:
    print ( "There are %i jobs that have not been done" % len( not_done_list ) )
    print ( "Run resubmit.py" )
    gen_submit_script( not_done_list, fout = 'resubmit.py' )
"""
    f = open( top_dir + '/VASPRUNS/check.py', 'w' )
    f.write(s)
    f.close()

def gen_submit_script( not_done_list, fout = 'submit.py' ):
    s = """import os
submit_cmd = 'sbatch' # replace this with an appropriate one

list = %s
if __name__ == "__main__":
    for item in sorted( list ):
        folder = item.replace( 'POSCAR-', 'disp-' ).replace( 'SPOSCAR', 'disp-S' )
        os.system( 'rm -fv ' + folder )
        os.system( 'mkdir -p ' + folder )
        os.chdir ( folder )
        os.system( 'ln -s ../INCAR' )
        os.system( 'ln -s ../POTCAR' )
        os.system( 'ln -s ../KPOINTS' )
        os.system( 'ln -s ../' + item + ' POSCAR' )
        os.system( 'cat ../job | sed \\'s/$sys/'+ item.replace( 'POSCAR','' ) +'/\\' > job' )
        os.system( submit_cmd +  ' job' )
        os.chdir ( '../' )
    """ %( not_done_list )
    
    f = open( fout, 'w' )
    f.write( s )
    f.close()
    return

def gen_collect_vasprun_script( top_dir, fixed_indices, disp_indices ):
    str = '''import os 
from dlePy.phonopy.generate import convert_vasprun
for item in os.listdir( './' ):
    if ( 'POSCAR' in item ) and ( len( item ) > 6 and ('APOSCAR' not in item) ):
        id = item.replace( 'POSCAR-', '' ).replace( 'SPOSCAR', 'S' )
        path = 'VASPRUNS/disp-' + id + '/'
        file= path + 'vasprun.xml'
        fileout= 'vasprun_' + id + '.xml'
        if id != 'S':
            fileout= 'vasprun_' + str( int( id ) ) + '.xml'
        # list of indices of fixed atoms
        fixed_indices = %s 
        disp_indices = %s 
        convert_vasprun( file, fileout, fixed_indices, disp_indices )

print ("""###########################################
After this, run the following:
phonopy -f data/vasprun_{XXX..YYY}.xml  # Replace XXX and YYY by appropriate numbers
phonopy -p mesh.conf

In addition, if BORN file is needed, convert_BORN_file.py script is available. 

    """)
    ''' %( fixed_indices, disp_indices )
    f = open( top_dir + '/collect_vasprun.py', 'w' )
    f.write( str )
    f.close()
    return

def convert_vasprun( file, fileout, fixed_indices, disp_indices ):
    natoms = len( fixed_indices ) + len( disp_indices )
    outdir = 'data/'
    os.system( 'mkdir -p ' + outdir )
    if len( fixed_indices )  == 0:
        os.system( 'cp -v ' + file + ' ' + outdir + fileout )
        return
    print ( 'Converting: ', file + ' -> ' + outdir + fileout )
    if os.path.isfile( file ):
        f_in = open( file, 'r' )
        f_out = open( outdir + fileout, 'w' )
        for i in range( 1000000 ):
            line = f_in.readline( )
            f_out.write( line.strip() + '\n' )
            if line.strip() == '<varray name="positions" >':
                for j in range( natoms ):
                    line = f_in.readline( )
                    if j in disp_indices:
                        f_out.write( line.strip() + '\n' )

            if line.strip() == '<varray name="forces" >':
                for j in range( natoms ):
                    line = f_in.readline( )
                    if j in disp_indices:
                        f_out.write( line.strip() + '\n' )

            if line.strip() == "":
                break

        f_in.close()
        f_out.close()
    else:
        print ( file, 'NOT FOUND' )

    return

'''
def convert_vasprun( file, fileout, fixed_indices, disp_indices ):
    natoms = len( fixed_indices ) + len( disp_indices )
    outdir = 'data/'
    os.system( 'mkdir -p ' + outdir )
    if len( fixed_indices )  == 0:
        os.system( 'cp -v ' + file + ' ' + outdir + fileout )
        return
    print ( 'Converting: ', file + ' -> ' + outdir + fileout )
    if os.path.isfile( file ):
        f_in = open( file, 'r' )
        f_out = open( outdir + fileout, 'w' )
        for i in range( 1000000 ):
            line = f_in.readline( )
            f_out.write( line.strip() + '\n' )
            if line.strip() == '<field type="int">atomtype</field>':
                line = f_in.readline( )
                f_out.write( line.strip() + '\n' )
                for j in range( natoms ):
                    line = f_in.readline( )
                    if j in disp_indices:
                        f_out.write( line.strip() + '\n' )

                    if line.strip() == '<varray name="positions" >':
                        for j in range( natoms ):
                            line = f_in.readline( )
                            if j in disp_indices:
                                f_out.write( line.strip() + '\n' )

                    if line.strip() == '<varray name="forces" >':
                        for j in range( natoms ):
                            line = f_in.readline( )
                            if j in disp_indices:
                                f_out.write( line.strip() + '\n' )

                    if line.strip() == "":
                        break

                f_in.close()
                f_out.close()
            else:
                print ( file, 'NOT FOUND' )

    return
'''
def gen_get_BORN_script( top_dir, fixed_indices, disp_indices ):
    s = '''
from dlePy.phonopy.generate import get_BORN

original_BORN = 'PATH_TO_ORIGINAL_BORN' # with all atoms, path to calculation. Thr program will read OUTCAR or OUTCAR.gz

fixed_indices = %s 
disp_indices = %s 
get_BORN( original_BORN, fixed_indices, disp_indices )

    ''' %( fixed_indices, disp_indices )
    f = open( top_dir + '/get_BORN_file.py', 'w' )
    f.write( s )
    f.close()
    return

def convert_BORN( original_BORN, fixed_indices, disp_indices ):
    print ( """ This is wrotking, however, will be replace by `get_BORN` """ )
    f_in = open( original_BORN, 'r' )
    f_out = open( 'BORN', 'w' )
    natoms = len( fixed_indices ) + len ( disp_indices )
    line = f_in.readline ()
    s = '# epsilon and Z* of atoms ',
    for at in disp_indices:
        s += ' ' + str( at + 1 )
    s += '\n'
    f_out.write( s )

    line = f_in.readline ()
    f_out.write( line.strip() + '\n' )
    for i in range( natoms ):
        line = f_in.readline()
        if i in disp_indices:
            f_out.write( line.strip() + '\n' )
    f_in.close()
    f_out.close()
    return

def get_BORN( original_BORN, fixed_indices, disp_indices, read_from_vasprun = False ):
    if not read_from_vasprun:
        if os.path.isfile( original_BORN + '/OUTCAR' ):
            f = open( original_BORN + '/OUTCAR', 'r' )
        elif os.path.isfile( original_BORN + '/OUTCAR.gz' ):
            f = gz.open( original_BORN + '/OUTCAR.gz', 'r' )
        else:
            print ( """ NO OUTCAR or OUTCAR.gz found in %s """ %( original_BORN ) )
            print ( """ EXITTING """ )
            exit( )
    else:
        if os.path.isfile( original_BORN + '/vasprun.xml' ):
            f = open( original_BORN + '/vasprun.xml', 'r' )
        elif os.path.isfile( original_BORN + '/vasprun.xml.gz' ):
            f = gz.open( original_BORN + '/vasprun.xml.gz', 'r' )
        else:
            print (""" NO vasprun.xml or vasprun.xml.gz found in %s """ %( original_BORN ))
            print (""" EXITTING """)
            exit( )

    lines = f.readlines( )
    f_out = open( 'BORN', 'w' )
    natoms = len( fixed_indices ) + len ( disp_indices )
    s = '# epsilon and Z* of atoms '
    for at in disp_indices:
        s += ' ' + str( at + 1 )
    s += '\n'
    f_out.write( s )

    if not read_from_vasprun:
        for iline in range( len( lines ) - 1, -1, -1 ):
            if 'MACROSCOPIC STATIC DIELECTRIC TENSOR ' in lines[ iline ]:
                s = '   '
                dielectrictensor = s.join( lines[ iline + 2 ].split( ) + lines[ iline + 3 ].split( ) +  lines[ iline + 4 ].split( ) )
                f_out.write( dielectrictensor.strip() + '\n' )
                break
    else:
        for iline in range( len( lines ) - 1, -1, -1 ):
            if '<varray name="epsilon" >' in lines[ iline ]:
                s = '   '
                dielectrictensor = s.join( lines[ iline + 1 ].replace('<v>','').replace('</v>','').split( ) + lines[ iline + 2 ].replace('<v>','').replace('</v>','').split( ) +  lines[ iline + 3 ].replace('<v>','').replace('</v>','').split( ) )
                f_out.write( dielectrictensor.strip() + '\n' )
                break


    if not read_from_vasprun:
        for iline in range( len( lines ) - 1, -1, -1 ):
            if 'BORN EFFECTIVE CHARGES ' in lines[ iline ]:
                start_line = iline
                break

        for iline in range( start_line, len( lines ) ):
            s = '   '
            if 'ion' in lines[ iline ].split()[0] and int ( lines[ iline ].split()[-1] ) - 1 in disp_indices:
                f_out.write( ' ' )
                borncharge = s.join( lines[ iline + 1 ].split()[ 1:: ] + lines[ iline + 2 ].split()[ 1:: ]  + lines[ iline + 3 ].split()[ 1:: ]  )
                f_out.write( borncharge.strip() + '\n' )
                if int( lines[ iline ].split()[ -1 ] ) - 1 == disp_indices[ -1 ]:
                    break
    else:
        for iline in range( len( lines ) - 1, -1, -1 ):
            if '<array name="born_charges" >' in lines[ iline ]:
                start_line = iline
                break
        for ion in disp_indices:
            l1 = iline + 1 + 5 * ion + 2
            l2 = iline + 1 + 5 * ion + 3
            l3 = iline + 1 + 5 * ion + 4
            dataline = lines[ l1 ].strip() + ' ' + lines[ l2 ].strip() + ' ' + lines[ l3 ].strip()
            dataline = dataline.replace('<v>','').replace('</v>','')
            s = '   '
            borncharge = s.join( dataline.split() )
            f_out.write( borncharge.strip() + '\n' )

    f_out.close()
    return

if __name__ == "__main__":
    print ( "Analizing POSCARs to see if any of them have been calculated in previous directory" )
    last_run_dir = '../../FIRST_RUN/VASPRUNS/'
    map_done, not_done_list = analize( last_run_dir )
    
    #This is for mapping the one that has already run to this folder"
    #link_last_run( map_done, last_run_dir = '../../FIRST_RUN/VASPRUNS/' ) 
    
    #To submit the remaining"
    submit( not_done_list )
