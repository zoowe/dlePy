from ase.build import bulk
from dlePy.qe.pwscf import PWscfInput, write_pwscf_input, update_keyword

# Create your structure
latt = 5.65
Gebulk = bulk( 'Ge', 'diamond', a = latt )

# Create `pwscf` object for `Gebulk`
pwscf = PWscfInput ( Gebulk )

# Write input file
write_pwscf_input ( pwscf , 'input_default.inp' )

