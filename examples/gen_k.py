"""
This example is to generate npts=38 kpoints from
Gamma -> K -> M -> Gamma, each KPOINTS file has 
ppf=5 k-points. 

This code will adjust the npts.
"""

from dlePy.vasp.kpoints import vasp_band_gen_k

klist=(
       ( 0.0,    0.0,   0.0 ),   #Gamma
       ( 2/3.,   1/3.,  0.0 ),      #K
       ( 0.5,    0.0,   0.0 ),     #M
       ( 0.0,    0.0,   0.0 )    #Gamma
      )

npts=38

ppf=5

vasp_band_gen_k(klist,npts,ppf, hse = False )


