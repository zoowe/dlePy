from ase import Atoms, Atom
from ase.io.vasp import read_vasp_out
from ase.io import Trajectory
from ase.calculators.singlepoint import SinglePointCalculator
from ase.visualize import view
import numpy as np
import gzip as gz
import os
import random as rd
import json as js
from tqdm import tqdm
from ..small_tools import str_decode

def read_OUTCAR( outcar, data = {}, maxiter = 1000 ):
    print ( "Parsing %s" %( outcar ) )
    try:
        if '.gz' in outcar:
            decode = True
            with gz.open( outcar, 'rb' ) as f:
                lines = f.readlines( )
        else:
            decode = False
            with open( outcar, 'r' ) as f:
                lines = f.readlines( )
    except:
        print ( "Error with reading %s" %( outcar ) )
        exit( )

    try:  # try to read constraints, first from CONTCAR, then from POSCAR
        constr = read_vasp(path + 'CONTCAR').constraints
    except Exception:
        try:
            constr = read_vasp(path + 'POSCAR').constraints
        except Exception:
            constr = None

    natoms = 0
    images = []
    atoms = Atoms(pbc=True, constraint=constr)
    energy = 0
    species = []
    species_num = []
    stress = None
    symbols = []
    ecount = 0
    poscount = 0
    magnetization = []
    imagecount = 0
    chargesMod  = []

    iter = len( data.keys() )
    ZVAL = {}
    
    FERMI_SHIFT = 0
    VACPOT_PSP  = 0
    E_fermi     = 0
 
    for n in tqdm( range( len( lines ) ) ):
        line = str_decode( lines[ n ], decode ) #.decode( 'utf-8' )
        if 'TITEL  =' in line:
            temp = line.split()[3].split( '_' )[0]
            species += [temp]
        
        if 'mass and valenz' in line:
            tmp = line.split( )[5]
            ZVAL[ species[ -1 ] ] = float( tmp )


        if 'ions per type' in line:
            #species = species[:len(species) // 2]
            temp = line.split()
            ntypes = min(len(temp)-4, len(species))
            for ispecies in range(ntypes):
                species_num += [int(temp[ispecies + 4])]
                natoms += species_num[-1]
                for iatom in range(species_num[-1]):
                    symbols += [species[ispecies]]

        if 'NELECT =' in line:
            temp = line.split()
            nElect = float(temp[2])

        if 'direct lattice vectors' in line:
            cell = []
            for i in range(3):
                temp = lines[n + 1 + i].split()
                cell += [[float(temp[0]), float(temp[1]), float(temp[2])]]
       
        if 'total charge    ' in line:
            charges = []
            charge = 0
            for iatom in range(natoms):
                temp = lines[n + 4 + iatom].split()
                charges += [ float(temp[4]) ]

            charges = np.array(charges)            

            # For future reference this total charge value below
            # from the file must be rounded becuase it does not give an accurate value
            # as compared to individually summing up the charge values so will NOT be using it
            # as the total charge value
            charge = float( lines[n + 4 + natoms + 1].split()[4] )

            # Modifier that will charge the charge values so they all equal the 
            # total number of electrons in the system
            electMod = nElect/np.sum(charges)
            chargesMod = charges * electMod
            
            for iatom in range(natoms):
                chargesMod[iatom] = chargesMod[iatom] - ZVAL[ symbols[iatom] ] 

        if len( chargesMod ) == 0:
            chargesMod  = [0]* natoms

        if 'FERMI_SHIFT' in line:
            tmp = float( line.split( )[2] )
            if np.abs( tmp ) > 0.000001:
                FERMI_SHIFT = tmp
        if 'VACPOT_PSP' in line: 
            tmp = float( line.split( )[2] )  
            if np.abs( tmp ) > 0.000001:  
                VACPOT_PSP = tmp  
        if 'E-fermi' in line:
            tmp = float( line.split( )[2] )
            E_fermi = tmp

        if 'POSITION          ' in line:
            imagecount += 1
            niter = get_niter( lines, n, decode )
            if niter < maxiter:
                forces = []
                positions = []
                atoms = Atoms(pbc=True, constraint=constr)
                atoms.set_cell(cell)
                for iatom in range(natoms):
                    temp = lines[n + 2 + iatom].split()
                    atoms += Atom(symbols[iatom],
                                  [float(temp[0]), float(temp[1]), float(temp[2])])
                    forces += [[float(temp[3]), float(temp[4]), float(temp[5])]]
                    positions += [[float(temp[0]), float(temp[1]), float(temp[2])]]

                energy = float( lines[n + 2 + natoms + 12].split()[6] )

                atoms.set_initial_charges(charges=chargesMod)

                atoms.set_calculator(SinglePointCalculator(atoms,
                                                           energy=energy,
                                                           forces=forces,
                                                           stress=stress,
                                                           charges=chargesMod))
                data_ = {}
                data_[ 'E-fermi' ] = E_fermi
                data_[ 'FERMI_SHIFT' ] = FERMI_SHIFT
                data_[ 'VACPOT_PSP' ] = VACPOT_PSP
                for x in range( n + 2 + natoms + 13, n + 2 + natoms + 12 + 50, 1 ):
                    line = lines[ x ].decode( 'utf-8' )
                    if 'ion-electron   TOTEN' in line:
                        data_ [ 'TOTEN' ] = float( line.split( "=")[-1].split()[0] )
                    if 'kinetic energy EKIN' in line:
                        data_[ 'EKIN' ]   = float( line.split( "=")[-1].split()[0] )
                    if 'kin. lattice  EKIN_LAT' in line:
                        data_[ 'EKIN_LAT' ] = float( line.split( "=")[-1].split()[0] )
                        data_[ 'TEMP' ]     = float( line.split( "temperature")[-1].split()[0] )
                    if 'nose potential ES' in line:
                        data_[ 'ES'  ]       = float( line.split( "=")[-1].split()[0] )
                    if 'nose kinetic   EPS' in line:
                        data_[ 'EPS' ]       = float( line.split( "=")[-1].split()[0] )
                    if 'total energy   ETOTAL' in line:
                        data_[ 'ETOTAL' ]  = float( line.split( "=")[-1].split()[0] )
                        break
                
                data[ int( iter ) ] =  data_
                iter += 1
 
                #print(len( atoms), atoms.get_number_of_atoms(), energy)
                images += [atoms]
            else:
                print ( "WARNING: Image %i has %i iterations and is removed from database" %( imagecount, niter ) )

    
    return images, data

def add_to_traj( traj, images, first, last = 500):

    if last > len( images ):
        last = len( images )
    for i in range( first, last ):
        traj.write( images[ i ] )

def get_niter( lines, endline, decode ):
    finish = False
    i = endline
    while not finish:
        i -= 1
        line = str_decode( lines[ i ], decode ) #.decode( 'utf-8 ' )
        if 'Iteration' in line:
            finish = True
            niter = line.replace( '-', '' ).replace( '(', ' ' ).replace( ')', ' ' ).split()[ -1 ]
            niter = int( niter )
    return niter

if __name__ == "__main__":
    images, data = read_OUTCAR( 'OUTCAR.gz'  )
    with open( 'DATA.js', 'w' ) as f:
        js.dump( data, f )
    traj = Trajectory('MOVIE.traj', 'w')
    add_to_traj( traj=traj, images=images, first = 0, last = 1000)
