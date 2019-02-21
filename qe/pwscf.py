"""
PWscfInput class was design to speed up  
input-generation processes for QUANTUM ESPRESSO 
runs. 

***Required*** 
Atomic Simulation Environment (ASE), numpy

***How to use*** 
If you know ASE, you know how to use this class.
See the provided example (make_pwscf_input.py) for details.

***Note***
This class is distributed in the hope
that it will benefit QE and/or ASE users. 
There is NO WARRANTY of any form. Users are supposed
to carefully check the class and its obtained results.

***History***
+++ The python class was started back in 2010. 
+++ Its author did not keep track of what have been changed.
+++ The code was last modified on 06/19/2013 for cleaning
up before being distributed. It was renamed to PWscfInput.

 
****************************************
Author: Duy Le
Department of Physics
University of Central Florida

Email: duy.le@ucf.edu
Website: http://www.physics.ufc.edu/~dle
****************************************
"""

import numpy as np
from ase.constraints import FixAtoms, FixScaled
from ..small_tools import file_exist

class PWscfInput:
    def __init__(self, atoms):
        """
        Used multiple sub-classes

        """
        self.atoms          = atoms
        self.control        = Control ()
        self.system         = System ( atoms )
        self.electrons      = Electrons ()
        self.ions           = Ions ()
        self.cell           = Cell ()
        self.atomic_species = AtomicSpecies ( self.atoms )
        self.kpoints        = Kpoints ( )

    def write_input ( self , filename):
        write_pwscf_input ( self , filename ) 

class Control:
    def __init__(self):
        """
        &CONTROL section
        """
        self.settings       = ControlSettings ()
        self.io             = ControlIO ()
        self.ion_relax      = ControlIonRelax ()

class ControlSettings:
    def __init__(self):
        self.calculation    = 'scf'
        self.restart_mode   = 'from_scratch'
        self.prefix         = 'pwscf'
        self.pseudo_dir     = 'PATH_TO_PSEUDO_DIR'

class ControlIO:
    def __init__(self):
        self.outdir         = './tmp/'
        self.verbosity      = 'low'
        self.disk_io        = 'none'
        self.wf_collect     = True 

class ControlIonRelax:
    def __init__(self):
        self.tprnfor        = True
        self.tstress        = True 
        self.forc_conv_thr  = 1.0e-3
        self.etot_conv_thr  = 1.0e-4
        self.nstep          = 100

class System:
    def __init__(self, atoms ): 
        self.structure      = SystemStructure ( atoms )
        self.ecut           = SystemEcut ( )
        self.occupations    = SystemOccupations ( )
        self.spin_pol       = SystemSpinPol ( )
        self.starting_magnetization  \
                            = StartingMagnetization \
                              (  self.structure.ntyp )

 
class SystemStructure:
    def __init__(self , atoms ):
        self.ibrav          = 0
        self.a              = 1.0 
        self.nat            = len(atoms)
        mol                 = get_reduce_atom_list ( atoms )
        self.ntyp           = len ( mol )
        # incorrect if atoms are not grouped.
        # self.nbnd= 100

class SystemEcut:
    def __init__(self):
        self.ecutwfc        = 50 
        self.ecutrho        = 200 

class SystemOccupations:
    def __init__(self):
        self.occupations    = 'smearing'
        self.smearing       = 'gaussian'
        self.degauss        = 0.0073498618

class SystemSpinPol:
    def __init__(self):
        self.nspin          = 1

class StartingMagnetization:
    def __init__(self , ntyp ):
        self.starting_magnetization \
                            = [ 0. ] * ntyp

class Electrons:
    def __init__(self ):
        self.diagonalization \
                            = 'david'
        self.conv_thr       = 1.0e-6
        self.electron_maxstep \
                            = 100
        self.mixing_mode    = 'plain'
        self.mixing_beta    = 0.7

class Ions:
    def __init__(self ):
        self.ion_dynamics   = 'bfgs'
        self.wfc_extrapolation \
                            = 'none'
        self.pot_extrapolation \
                            = 'atomic'

class Cell:
    def __init__(self ):
        self.press_conv_thr  = 0.5

class AtomicSpecies:
    def __init__(self , atoms ):
        mol                 =  get_reduce_atom_list ( atoms ) 
        ntyp                = len ( mol )
        self.ntyp           = ntyp
        self.symbol         = mol
        mass                = np.zeros ( [ ntyp ] , dtype = np.float )
        self.mass           = mass
        pseudo_potential    = np.zeros ( [ ntyp ] , dtype='S30')
        pseudo_potential [ : ] = 'Please_set_pseudo_file'
        self.pseudo_potential = pseudo_potential 

class Kpoints:
    def __init__( self ):
        self.type           = 'automatic' 
        self.mesh           = [ 1 , 1 , 1 ]
        self.smesh          = [ 0 , 0 , 0 ]
           
def update_keyword( input_block, keyword, value ):
    setattr ( input_block, keyword, value )

def get_reduce_atom_list ( atoms ):
    """
    Get list of atomic symbol then reduce it.
    New version of ASE should have this option 
    already.
    """
    mol =  atoms.get_chemical_symbols( )
    if len ( mol ) > 1:
        for i in range (len( mol )-1, 0 , -1):
            for j in range ( i ):
                if mol [ j ] == mol [ i ]:
                    del mol [ i ]
                    break
    return mol

def write_k_points ( kpoints , f ):
    if kpoints.type.lower() == 'gamma':
        kpoints.mesh[:]=1
        kpoints.smesh[:]=0
        
    print >>f, '! .kpoints' 
    print >>f, 'K_POINTS '+ kpoints.type  + '   ! .kpoints.type '
    print >>f, "%4i %4i %4i %3i %3i %3i %s" % ( 
               kpoints.mesh[0] ,  kpoints.mesh[1] ,  kpoints.mesh[2], 
               kpoints.smesh[0] ,  kpoints.smesh[1] ,  kpoints.smesh[2], '! .kpoints.mesh and .kpoints.smesh' 
               )

def write_key ( item , dict ):
    value =   vars( dict )[item]
    if type ( value  ) is str:
        str_value   = '\''+value+'\','
    if type ( value )  is float:
        str_value   = str ( value )+','
    if type ( value ) is int:
        str_value   = str ( value )+','
    if type ( value ) is bool:
        if ( value ):
            str_value  = '.true.,'
        else:
            str_value  = '.false.,'
    item_len = item.__len__()
    default_len = 25 
    add_str = ''
    for i in range ( item_len , default_len):
        add_str +=' ' 
    string  = item + add_str + ' = ' + str_value
    return string

def write_array_key ( item , dict , f ):
    array_value =   vars( dict )[item]
    for i in range (len (array_value) ):
        value = array_value [ i ]
        item_len = item.__len__()
        default_len = 25
        add_str = ''
        for j in range ( item_len + 3, default_len):
            add_str +=' '
 
        string  = item + '('+str(i+1)+')'+add_str + ' = ' 
        print >> f, string, value, ','

def write_atomic_species ( atomic_species , f ):
    if len( atomic_species.mass ) != atomic_species.ntyp:
        print 'ERROR: len( mass ) != ntyp'
        exit( 'EXITTING' )
    if len( atomic_species.pseudo_potential ) != atomic_species.ntyp:
        print 'ERROR: len( pseudo_potential ) != ntyp'
        exit( 'EXITTING' )

    for i in range ( atomic_species.ntyp ):
        print >> f, "%5s %8.4f %s" % ( \
                    atomic_species.symbol [ i ] , \
                    atomic_species.mass [ i ] , \
                    atomic_species.pseudo_potential [ i ]
                    )

def write_structure ( atoms, f, ibrav):

    print >>f,'ATOMIC_POSITIONS  crystal '
    sflags = np.zeros((len(atoms), 3), dtype=bool)
    newsflags = np.ones((len(atoms), 3), dtype=np.int)
    if atoms.constraints:
        for constr in atoms.constraints:
            if isinstance(constr, FixScaled):
                sflags[constr.a] = constr.mask
            elif isinstance(constr, FixAtoms):
                sflags[constr.index] = [True, True, True]
            for i in range(len(atoms)):
                for j in range(3):
                    if  sflags[i,j]:
                        newsflags[i,j]=0
                
    for i in range(len(atoms)):
        x = atoms.get_scaled_positions()[i,:]
        for j in range(3):
            if x[j] > 0.5:
                x[j]-=1.
        print >>f, '%3s %20.14f %20.14f %20.14f %3i %3i %3i' %( \
                    atoms.get_chemical_symbols()[i], \
                    x[0] , x[1] , x[2],  \
                    newsflags[i,0] , newsflags[i,1], newsflags[i,2] )
#                    atoms.get_scaled_positions()[i,0], \
#                    atoms.get_scaled_positions()[i,1], \
#                    atoms.get_scaled_positions()[i,2], \
#                    newsflags[i,0] , newsflags[i,1], newsflags[i,2] )

    if ibrav == 0: 
        print >>f,'CELL_PARAMETERS'
        for i in range (3):
            print >>f,'%20.14f %20.14f %20.14f' %( \
                 atoms.cell[i,0], \
                 atoms.cell[i,1], \
                 atoms.cell[i,2]) 

def verify_potential( object ):
    pseudo_dir = object.control.settings.pseudo_dir
    ecutwfc = []
    ecutrho = [] 
    for i in range ( object.atomic_species.ntyp ):
        pot = object.atomic_species.pseudo_potential [ i ]
        exist, txt = file_exist( pseudo_dir + '/' + pot )
        if not exist:
            print 'WARNING :' + txt
        else:
            with open( pseudo_dir + '/' + pot, 'r' ) as fpot:
                lines = fpot.readlines( )
            for line in lines:
                if 'Suggested minimum cutoff for wavefunctions' in line:
                    ecutw = float( line.split()[ - 2] )    
                    ecutwfc.append( ecutw )
                    break
            for line in lines:
                if 'Suggested minimum cutoff for charge density' in line:
                    ecutc = float( line.split()[ - 2] )    
                    ecutrho.append( ecutc )
                    break
    if len( ecutwfc ) * len( ecutrho ) > 0:
        print "*************************"
        print "Suggested minimum cutoff for wavefunctions: ecutwfc = %4.1f Ry" %( np.max( np.array( ecutwfc ) ) )
        print "Suggested minimum cutoff for charge density: ecutrho = %4.1f Ry" %( np.max( np.array( ecutrho ) ) )
        print "You are using ecutwfc = %4.1f Ry and ecutrho = %4.1f Ry" % ( object.system.ecut.ecutwfc, object.system.ecut.ecutrho )
    else:
        print "*************************"
        print "No suggestion for cutoff values as potential files cannot be read"

def write_pwscf_input ( object , filename, verify_pot = False):
    f = open ( filename, 'w' )
    """ Write CONTROL section """
    print >>f, '&CONTROL'
    print >>f, '! .control.settings'
    dict = object.control.settings
    for item in vars( dict ):
        print >>f, write_key ( item , dict )

    print >>f, ''
    print >>f, '! .control.io'
    dict = object.control.io 
    for item in vars( dict ):
        print >>f, write_key ( item , dict )
 
    if vars( object.control.settings )[ "calculation" ] in [ "relax", "vc-relax" ]:
        print >>f, ''
        print >>f, '! .control.ion_relax'
        dict = object.control.ion_relax
        for item in vars( dict ):
            print >>f, write_key ( item , dict )

    print >>f, '/'

    print >>f, ''

    """ &SYSTEM section """
    print >>f, '&SYSTEM'
    print >>f, '! .system.structure'
    dict = object.system.structure 
    for item in vars( dict ):
        print >>f, write_key ( item , dict )

    print >>f, ''
    print >>f, '! .system.ecut'
    dict = object.system.ecut
    for item in vars( dict ):
        print >>f, write_key ( item , dict )

    print >>f, ''
    print >>f, '! .system.occupations'
    dict = object.system.occupations
    for item in vars( dict ):
        print >>f, write_key ( item , dict )

    if object.system.spin_pol.nspin is 2:
        print >>f, ''
        print >>f, '! .system.spin_pol'
        dict = object.system.spin_pol
        for item in vars( dict ):
            print >>f, write_key ( item , dict )

        print >>f, '! .system.starting_magnetization'
        dict = object.system.starting_magnetization
        item = 'starting_magnetization'
        write_array_key ( item ,  dict , f)

    print >>f, '/'
    print >>f, ''

    """ &ELECTRONS section """
    print >>f, '&ELECTRONS' 
    print >>f, '! .electrons'
    dict = object.electrons
    for item in vars( dict ): 
        print >>f, write_key ( item , dict )

    print >>f, '/'
    print >>f, ''

    """ &IONS section """
    if vars( object.control.settings )[ "calculation" ] in [ "relax", "vc-relax", "md", "vc-md" ]:
        print >>f, '&IONS'
        print >>f, '! .ions'
        dict = object.ions
        for item in vars( dict ):
            print >>f, write_key ( item , dict )

        print >>f, '/'
        print >>f, ''

    """ &CELL section """
    if vars( object.control.settings )[ "calculation" ] in [ "vc-relax", "vc-md" ]:
        print >>f, '&CELL'
        print >>f, '! .cell'
        dict = object.cell
        for item in vars( dict ):
            print >>f, write_key ( item , dict )

        print >>f, '/'
        print >>f, ''

    """ ATOMIC_SPECIES section """
    print >>f, '! .atomic_species'
    print >>f, 'ATOMIC_SPECIES'
    write_atomic_species ( object.atomic_species , f )
    print >>f, ''
    write_structure      ( object.atoms, f, object.system.structure.ibrav)

    print >>f, ''
    write_k_points       ( object.kpoints, f)
    
    if verify_pot:
        verify_potential( object )
