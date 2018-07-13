"""
Replace XX with appropriate values, which are coordinate of L, G, X points
"""
from dlePy.vasp.kpoints import vasp_band_gen_k

klist=(
   ( XX,     XX,     XX), #     L
   ( XX,     XX,     XX), #     G
   ( XX,     XX,     XX), #     X
      )

npts = ppf = 31

vasp_band_gen_k(klist,npts,ppf)


