from ase.build import molecule, fcc111, add_adsorbate
from ase.io import write
from dlePy.strucmod import sort_group

slab = fcc111('Al', size=(2,2,3), vacuum=7.5)
mol  = molecule( 'H2O' )

add_adsorbate(slab, mol , 1.5, 'ontop')
add_adsorbate(slab, mol , 2.0, 'bridge')

indices = [ atom.index for atom in slab if atom.symbol in ['H', 'O'] ]
slab_sort_group_indices = sort_group( slab, indices = indices  )
slab_sort_group_symbols = sort_group( slab, symbols = [ 'H', 'O' ]  )

# Write diffent POSCAR files for the same system with different sort options.
write('POSCAR.no_sort', slab, format='vasp', vasp5 = True, direct = True, sort = False )
write('POSCAR.sorted',  slab, format='vasp', vasp5 = True, direct = True, sort = True )
write('POSCAR.sort_group_indices', slab_sort_group_indices, format='vasp', 
       vasp5 = True, direct = True, sort = False )
write('POSCAR.sort_group_symbols', slab_sort_group_symbols, format='vasp', 
       vasp5 = True, direct = True, sort = False )
