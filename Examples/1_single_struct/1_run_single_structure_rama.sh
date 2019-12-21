
## generate dihedral angle distribution of all residues in PDB
## according to their type: General AA, Glycine, Proline, 
##   one residue before Proline

../../1_rama_single_structure.py  \
  -in 3anr.pdb.bz2                \
  -img 3anr.rama_plot.png         \
  -ref ../../pyrama_data

## ../../pyrama_data
#  directory to datafiles of background dihedral angle distribution
#
## 3anr.pdb.bz2
#  input pdb file, zipped format is okay
#
## 3anr.rama_plot.png
#  result plots
#
#
#  19.12.20

