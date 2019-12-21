#!/bin/csh

# generate the distribution of torsion angles of a residue
# thoughout the MD trajectory in a heatmap format

../../2_rama_md_heatmap.csh   \
  2 11                        \
  fgf21-wt-s2.fasta           \
  200                         \
  fgf21-wt                    \
  png                         \
  templ_metrics.wt.traj       \
  1 1                         \
  ../..                       \
  /home/software/ctraj/bin/cpptraj


## fgf21-wt.A208.rama.out.bz2
#  output of the phi/psi angles of residue A208
#
## fgf21-wt.A208.rama_histo.png
#  heatmap of the distribution
#
## fgf21-wt-s2.fasta
#  Fasta sequence corresponding to the simulated peptide; no header
#
## fgf21-wt-s2.prot.prmtop
#  topology of peptide simulation
#
## fgf21-wt-s2.prot.gamd.nc
#  trajectory of peptide simulation
#
## templ_metrics.wt.traj
#  template CPPTRAJ input file, contains topology info and trajectory
#  info
#
## tmp.fgf21-wt.traj  
#  generated from the template 'templ_metrics.wt.traj' to extract
#  dihedral angles data from trajectory
#
#  19.12.20
