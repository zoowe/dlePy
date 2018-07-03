from ase.build import bulk
from dlePy.qe.pwscf import PWscfInput, write_pwscf_input, update_keyword

# Create your structure
latt = 5.65
Gebulk = bulk( 'Ge', 'diamond', a = latt )

# Create `pwscf` object for `Gebulk`
pwscf = PWscfInput ( Gebulk )

# The following are for changing default values, keyword

# Change calculation to scf
update_keyword( pwscf.control.settings, 'calculation', 'scf' )

# Update pseudo_dir
update_keyword( pwscf.control.settings, 'pseudo_dir', '/shared/ESPRESSO/PSLIBRARY/1.0.0/pbe/PSEUDOPOTENTIALS/' )

# Set mass and pseudo potential file for each type
mass = [ 72.6300 ]
pseudo_potential = [ 'Ge.pbe-n-kjpaw_psl.1.0.0.UPF' ]
update_keyword( pwscf.atomic_species, 'mass' , mass )
update_keyword( pwscf.atomic_species, 'pseudo_potential', pseudo_potential )

# Set k-point mesh. Only automatic 
update_keyword( pwscf.kpoints, 'mesh',  [ 15, 15, 15] )
update_keyword( pwscf.kpoints, 'smesh', [  0,  0, 0 ]  )

# write input file
write_pwscf_input ( pwscf , 'input_01.inp' )

