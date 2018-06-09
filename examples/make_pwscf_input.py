import numpy as np
from ase import Atoms 
from ase.build.surface import fcc111
from ase.constraints import FixAtoms
from dlePy.qe.pwscf import PWscfInput, write_pwscf_input

# Generate 5 layer of (1x1) Cu(111) in a 20 angstrom vacuum
system = fcc111('Cu',size=(1,1,5),a=3.16,vacuum=20./2.)

# Fix bottom 3 layer
indices = [ a.index for a in system \
            if a.position[2] < 0.51 * system.cell[2,2] ]
system.set_constraint ( FixAtoms ( indices = indices ) )
pwscf = PWscfInput ( system )

# Set pseudo_dir
pwscf.control.settings.pseudo_dir = '/home/PS_LIBRARY'

# Set high disk_io
pwscf.control.io.disk_io = 'high'

# set non spinpolarize
pwscf.system.spin_pol.nspin = 1

# Do not need it if nspin = 1. This is how to set starting_magnetization
pwscf.starting_magnetization.starting_magnetization[0] = 1.

# Set mass and pseudo potential file for each type
pwscf.atomic_species.mass[0] = 63.546 
pwscf.atomic_species.pseudo_potential[0] = 'Cu.pbe-n-nc.UPF' 

# Set k-point mesh. Only automatic 
pwscf.kpoints.nk = [ 3 , 3 , 1]
pwscf.kpoints.sk = [ 0 , 0 , 0]


write_pwscf_input ( pwscf , 'input.inp' )

# set more input
setattr ( pwscf.ions , 'upscale', 100.0 )
write_pwscf_input ( pwscf , 'input.1.inp' )

# remove some setting
del  pwscf.ions.upscale
del  pwscf.control.ion_relax.tstress
write_pwscf_input ( pwscf , 'input.2.inp' )

# change default setting
setattr ( pwscf.electrons , 'mixing_beta', 0.2 )
write_pwscf_input ( pwscf , 'input.3.inp' )


