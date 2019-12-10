#!/bin/csh

##########################################################################
## 
##  Peter MU Ung @ UMich COPharm
##  v1.0 -- 09.01.22
##
##  Peter MU Ung @ MSSM / Yale
##  v2.0 -- 19.12.03	generalize process up to generating figures
##
##  To generate PTRAJ input files and run PTRAJ to generate separate  
##  Phi and Psi angle files for generating Ramachandran plot.
##
##  Use "phi_psi_dih_combine.csh" to combine the separate phi/psi files
##  Phi+Psi file will be used for "rama_png_gen.pl" to generate Rama plot
##  in PNG format.
##
##  Use with:
##    -) cpptraj (use GitHub V4.14.12+ for extra functions)
##    -) combine_column.py
##    -) rama_md_heatmap_gen.py
##    -) ecoDnaK.fasta (only 1 line of the sequence)
##
##  Note: Torsion angle calculation will always start with i+1 and end i-1
##        because the first residue will not have Phi and last residue will
##        not have Psi angle
##
##
##########################################################################

if ($#argv != 10) then
  echo ''
  echo '  > x.csh '
  echo '      [ Starting residue for phi/psi in PTRAJ (prmtop) ]'
  echo '      [ Ending residue for phi/psi in PTRAJ (prmtop)   ]'
  echo '      [ fasta_file (one line of fasta seq) ]'
  echo '      [ Actual starting residue number ]'
  echo '      [ Output prefix ]'
  echo '      [ Template cpptraj input with "parm" and all "trajin" ]'
  echo '      [ Run cpptraj? : 0|1 ]'
  echo '      [ Generate figure with AA Dihedral density reference? : 0|1 ]'
  echo ''
  echo '      [ Path to this script directory ]'
  echo '      [ Path to GitHub cpptraj ]'
  echo ''
  exit
endif

set startResid = $argv[1]	# Starting residue for phi/psi in PTRAJ
set endResid   = $argv[2]	# Ending residue for phi/psi in PTRAJ
set fasta_file = $argv[3]
set fastaStart = $argv[4]

set out_pref   = $argv[5]	# Folder name, and file prefix
set img_ext    = $argv[6] # Image format. png,svg,jpg,eps,ps,pdf
set templ_traj = $argv[7]	# Template PTRAJ input file
set run_traj   = $argv[8]	# run dihedral data generation
set ref_rama   = $argv[9]	# use reference AA dihedral density map

#set scptdir = '/home/pmung/Dropbox/9_scripts/3_program/plotting/RAMAplot'
#set ctraj   = '/home/software/ctraj/bin/cpptraj'
set scptdir = $argv[10]
set ctraj   = $argv[11]

##########################################################################


## split fasta into array of characters
set fasta = (`cat ./$fasta_file | grep -o .`)

cp ./$templ_traj tmp.$out_pref.traj

## Write out all individual temp traj input files, then run CPPTRAJ
@ i = $startResid
@ a = $fastaStart
while ($i <= $endResid)	# Generate tmp input for Dihedral angle calc.
  @ j = $i - 1
  @ k = $i + 1

  ## Phi = C1-N2-CA2-C2
  echo dihedral $i\_phi :$j@C :$i@N :$i@CA :$i@C \
    out $out_pref.${fasta[$i]}$a.phi.out >> tmp.$out_pref.traj

  ## Psi = N2-CA2-C2-N3
  echo dihedral $i\_psi :$i@N :$i@CA :$i@C :$k@N \
    out $out_pref.${fasta[$i]}$a.psi.out >> tmp.$out_pref.traj

  @ i++
  @ a++
end

## generate phi-psi data 
if ($run_traj == 1) then
  $ctraj tmp.$out_pref.traj
endif


########################################################################


## combine phi/psi files into one, then convert to figures with scripts
@ i = $startResid
@ a = $fastaStart
while ($i <= $endResid)

  ## Generate phi-psi data
  echo " ## Generate combined phi-psi files: $out_pref.${fasta[$i]}$a"
  if (! -e $out_pref.${fasta[$i]}$a.rama.out.bz2) then
    $scptdir/combine_column.py \
      -in1 $out_pref.${fasta[$i]}$a.phi.out \
      -in2 $out_pref.${fasta[$i]}$a.psi.out \
      -cols 2                               \
      > $out_pref.${fasta[$i]}$a.rama.out
    echo $out_pref.${fasta[$i]}$a.rama.out

    bzip2 -f $out_pref.${fasta[$i]}$a.rama.out
    rm $out_pref.${fasta[$i]}$a.phi.out $out_pref.${fasta[$i]}$a.psi.out
  endif


  ## Set up parameters for Ramachandran density heatmap,
  ## use different reference Ramachandran densities (Pro, Gly, all others)
  ## according the the residue to be plotted
  echo ''
  echo " ## Generate Ramachandran plot: $out_pref.${fasta[$i]}$a"
  set t_step = '-t_step 1'
  set c_step = '-c_step 1'

  @ j = $i + 1

  if ($ref_rama == 0) then
    set ref = ''
    set res = ''
  else if (${fasta[$i]} == 'P') then
    set ref = "-ref $scptdir/pyrama_data/pref_proline.data.bz2"
    set res = '-res Pro'
  else if (${fasta[$i]} == 'G') then
    set ref = "-ref $scptdir/pyrama_data/pref_glycine.data.bz2"
    set res = '-res Gly'
  else if (${fasta[$j]} == 'P') then
    echo  "  * Pre-proline: ${fasta[$i]}$a"
    set ref = "-ref $scptdir/pyrama_data/pref_preproline.data.bz2"
    set res = '-res PreP'
  else
    set ref = "-ref $scptdir/pyrama_data/pref_general.data.bz2"
    set res = '-res Gen'
    set t_step = ''
    set c_step = ''
  endif

  ## Generate dihedral density heatmap
  $scptdir/2x_rama_md_heatmap_gen.py                  \
    -in  $out_pref.${fasta[$i]}$a.rama.out.bz2        \
    -img $out_pref.${fasta[$i]}$a.rama_histo.$img_ext \
    $ref $res $t_step $c_step

  @ i++
  @ a++
end
  

