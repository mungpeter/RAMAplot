# RAMAplot
_**Generate Ramachandran plots for single PDB structure and as heatmap for MD trajectory**_

There are 2 parts to this particular Ramachandran (amino acid backbone dihedral angle) plot scripts:
1) Generation of standarad Ramachandran plots for one input PDB structure.
2) Generation amino acid backbone dihedral angle distribution/population (from MD simulations) as contoured heatmap.

*_Reference 1_: [PyRAMA](https://github.com/gerdos/PyRAMA)

*_Reference 2_: [Lovell et al. Structure validation by Calpha geometry: phi, psi and Cbeta deviation. Proteins, (2003) 50(3): 437-450](https://doi.org/10.1002/prot.10286)

#######################################################################################
```
> rama_single_structure.py
      -in  [ PDB file for Ramachandran density plot ]
      -png [ Output PNG name ]

     optional:
      -int [ Resolution (def: 2-deg interval) ]
      -ref [ Directory to Density data for reference Ramachandran distribution ]
      -dpi [ PNG figure DPI resolution (def: 300) ]


e.g.> ./rama_single_structure_gen.py \
        -in 3anr.pdb.bz2             \
        -png 3anr.rama_plot.png      \
        -ref ../pyrama_data
```

This script manages the generation of Ramachandran plots for one input PDB structure, where amino acids are segregrated into **"Proline"**, **"Pre-Proline"** (1 residue before Proline), **"Glycine"**, and all other **"General"** amino acids. Each plot has a typical distribution of amino acid dihedral angle distribution in the background as a reference.
![Ramachandran plot of a PDB structure](https://github.com/mungpeter/RAMAplot/blob/master/1_example/3anr.rama_plot.png)

#######################################################################################

```
> rama_md_heatmap.csh
      [ Starting residue for phi/psi in PTRAJ (prmtop) ]
      [ Ending residue for phi/psi in PTRAJ (prmtop)   ]
      [ fasta_file (one line of fasta seq) ]
      [ Actual starting residue number ]
      [ Output prefix ]
      [ Template cpptraj input with "parm" and all "trajin" ]
      [ Run cpptraj? : 0|1 ]
      [ Generate figure with AA Dihedral density reference? : 0|1 ]
      [ Path to this script directory ]
      [ Path to GitHub cpptraj ]

*** Require the main processing script: rama_md_heatmap_gen.py ***

e.g.> ./rama_md_heatmap.csh    \
        2 11 fgf21-wt-s2.fasta \
        200                    \
        fgf21-wt               \
        templ_metrics.wt.traj  \
        1 1                    \
        .                      \
        /home/software/cpptraj/bin/cpptraj
```
This shell script calls up **rama_md_heatmap_gen.py** to generate the distribution/population of amino acid backbone dihedral angle as contoured heatmap in a Ramachandran plot. Data is usually from MD trajectories and can only do 1 residue at a time. Most useful for comparison of amino acid behaviour change in a protein with mutant residues, or to examine the difference in structural or conformational changes.

![Ramachandran plot of a PDB structure](https://github.com/mungpeter/RAMAplot/blob/master/1_example/fgf21-wt.A208.rama_histo.png)

#######################################################################################
Stable packages:
```
python      # 3.6.9+
numpy       # 1.17.4
pandas      # 0.25.3
matplotlib  # 3.0.2+
scipy       # 1.3.2
biopython   # 1.72+
gzip        #
bzip2       #
argparse    # 1.1
```

**Disclaimer:** These scripts are inspired by the PyRAMA repository by @gerdos so you will see resemblance of the figures. Redid it most to add a slightly difference artistic touch of my own, and to enable the heatmap generation for MD trajectory.
