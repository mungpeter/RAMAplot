# RAMAplot
_**Generate Ramachandran plots for single PDB structure and as heatmap for MD trajectory**_

```
  Author: Peter M.U. Ung @ MSSM / Yale

  vers:   1.0
```
There are 2 parts to this particular Ramachandran (amino acid backbone dihedral angle) plot scripts:
1) Generation of standarad Ramachandran plots for one input PDB structure.
2) Generation amino acid backbone dihedral angle distribution/population (from MD simulations) as contoured heatmap.

*_Reference 1_: [PyRAMA](https://github.com/gerdos/PyRAMA)

*_Reference 2_: [Lovell et al. Structure validation by Calpha geometry: phi, psi and Cbeta deviation. Proteins, (2003) 50(3): 437-450](https://doi.org/10.1002/prot.10286)

#######################################################################################
- **Convert PDB into its corresponding FASTA sequence**
```
> 0_pdb2fasta.py
      [PDB filename: .pdb]
      [Name of PDB]
      [FASTA Output Name]

  Note: Capping Groups (ACE/NME) and nonstandard amino acids are not recognized
        Need to add 'X' to the result FASTA manually

e.g.> ./0_pdb2fasta.py \
        test.pdb       \
        test           \
        test.fasta
```

This script convert PDB file into FASTA sequence. However, nonstandard amino acids such as capping groups (ACE/NME) and modified residues (PTO/MSE) and alternative names (HIE/HIP/HID for HIS) are not read in properly by BioPython. User need to add ACE/NME (as 'X') to the result fasta file; alternative residue name _(sed 's/HIE|/HIP|/HID/HIS/g')_ and non-standard AA fixings have to be done before hand.

#######################################################################################
- **Generate Ramachandran Plot for _all_ residues in _ONE_ PDB structure**
```
> 1_rama_single_structure.py
      -in  [ PDB file for Ramachandran density plot ]
      -img [ Output image name, extension as format (e.g. .png,svg,eps,ps,pdf) ]

     optional:
      -ref [ Directory to Density data for reference Ramachandran distribution ]
      -dpi [ Figure DPI resolution (def: 300) ]


e.g.> ./1_rama_single_structure_gen.py \
        -in 3anr.pdb.bz2               \
        -img 3anr.rama_plot.png        \
        -ref ../pyrama_data
```

This script manages the generation of Ramachandran plots for one input PDB structure, where residues are groups into **"Proline"**, **"Pre-Proline"** (1 residue before Proline), **"Glycine"**, and all other **"General"** amino acids. Each plot has a typical distribution of amino acid dihedral angle distribution in the background as a reference.

- Ramachandran plot of a PDB structure
![Ramachandran plot of a PDB structure](https://github.com/mungpeter/RAMAplot/blob/master/Examples/1_single_struct/3anr.rama_plot.png)

#######################################################################################
- **Generate Ramachandran Plot for _ONE_ residues in a MD trajectory**
```
> 2_rama_md_heatmap.csh
      [ Starting residue for phi/psi in PTRAJ (prmtop) ]
      [ Ending residue for phi/psi in PTRAJ (prmtop)   ]
      [ fasta_file (one line of fasta seq) ]
      [ Actual starting residue number in protein ]
      [ Output prefix ]
      [ Output image extension, (e.g. png,svg,eps,ps,pdf) ]
      [ Template cpptraj input with "parm" and all "trajin" ]
      [ Run cpptraj? : 0|1 ]
      [ Generate figure with AA Dihedral density reference? : 0|1 ]
      [ Path to this script directory ]
      [ Path to GitHub cpptraj ]

*** Require the main processing script: rama_md_heatmap_gen.py ***

e.g.> ./2_rama_md_heatmap.csh    \
          2 11                   \
          fgf21-wt-s2.fasta      \
          200                    \
          fgf21-wt               \
          png                    \
          templ_metrics.wt.traj  \
          1 1                    \
          .                      \
          /home/software/cpptraj/bin/cpptraj
```
This shell script calls up **2x_rama_md_heatmap_gen.py** to generate the distribution/population of amino acid backbone dihedral angle as contoured heatmap in a Ramachandran plot. Data is usually from MD trajectories and can only do 1 residue at a time. Most useful for comparison of amino acid behaviour change in a protein with mutant residues, or to examine the difference in structural or conformational changes.

```
> 2x_rama_md_heatmap_gen.py
      -in  [ phi-psi file for Ramachandran density plot ]
      -img [ output image name, extension as format (e.g. .png,svg,eps,ps,pdf) ]
      
    optional:
      -res [ 1-char AminoAcid code for reference selection (def: Pro, PreP, Gly, Gen) ]
      -int [ Resolution (def: 2-deg interval) ]
      -ref [ Density data for reference Ramachandran distribution ]
      -fraction [ Cutoff fraction of the maximum Histogram value (def: 33) ]
      -smooth   [ Histogram data smoothening (def: 1.15) ]
      -t_step   [ Colorbar tick spacing per histogram digits unit (def: 4) ]
      -c_step   [ Histo Contour spacing per histogram digits unit (def: 4) ]
      -dpi      [ Figure DPI resolution (def: 300) ]\n
```
This is the **actual** script that generates the heatmap figure, but it can only do **one** residue at a time so it should be coupled with the above _rama_md_heatmap.csh_ script to run through a list of residues.

- Distribution of AA backbone dihedral angle of a residue throughout MD trajectories as heatmap
![Distribution of AA backbone dihedral angle of a residue throughout MD trajectories as heatmap](https://github.com/mungpeter/RAMAplot/blob/master/Examples/2_md_heatmap/fgf21-wt.P205.rama_histo.png)

#######################################################################################
- **Extract percentage of population within a range of x- and y-axes**
```
> 3_rama_extract_popul.py
      -in      [ phi-psi file for Ramachandran density plot ]
      -list    [ list of ranges in degree to extraction population ]
      
    Optional:
      -pwd     [ path to the data file ]
      -xrange  [ Phi-range, if not using -list ]
      -yrange  [ Psi-range, if not using -list ]

#####
# example input for -list file: 
  # -180 -120 , 120 180     # population_1: Phi-range [-180,-120] ; Psi-range [120,180]
  # -120 -30  ; 120 180     # population_2: Phi-range [-120,-30] ;  Psi-range [120,180]
  # 30 90     : -30 60      # population_2: Phi-range [30,90] ;     Psi-range [-30,60]
#####  

e.g.> ./3_rama_extract_popul.py                  \
          -in   fgf21-wt.A208.rama.out.bz2       \
          -list fgf21-wt.A208.rama.cluster.txt
          
---------------------------------
# X-range [ -180.0 , -120.0 ] | Y-range [ 120.0 , 180.0 ]
 % Population: 11.62

# X-range [ -120.0 , -30.0 ] | Y-range [ 120.0 , 180.0 ]
 % Population: 31.39

# X-range [ -120.0 , -30.0 ] | Y-range [ -60.0 , 30.0 ]
 % Population: 21.10

# X-range [ 30.0 , 90.0 ] | Y-range [ -30.0 , 60.0 ]
 % Population: 4.36

# % Total Population: 68.48

```
This script calculate the population of dihedral angles that fall within the defined ranges. 

- Population of different clusters of dihedral angles in MD trajectory
![Population of different clusters of dihedral angles in MD trajectory](https://github.com/mungpeter/RAMAplot/blob/master/Examples/3_extract_popul/fgf21-wt.A208.rama_histo.popul.png)


#######################################################################################
- **Stable packages:**
```
csh/tcsh      # shell

python        # 3.6.9+
  numpy       # 1.17.4
  pandas      # 0.24.1+
  matplotlib  # 3.0.2+
  scipy       # 1.3.2
  biopython   # 1.72+
  gzip        #
  bzip2       #
  argparse    # 1.1
```

**Disclaimer:** These scripts are inspired by the PyRAMA repository by @gerdos so you will see resemblance of the figures. Redid it mostly to add a slightly different artistic touch of my own, and to enable the heatmap generation for MD trajectory.
