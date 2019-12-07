#!/usr/bin/env python3

#############################################################################
##
##  Peter MU Ung @ MSSM / Yale
##
##  v2.0 19.12.03 - generate figure with Seaborn/Matplotlib directly
##
##  use to convert dihedral angle of residues of a PDB structure 2D 
##  Ramachandran plot
##
##
############################################################################

import sys,os,gzip,bz2,re
import aa_residue
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
p = PDBParser(PERMISSIVE=1,QUIET=True)

from argparse import ArgumentParser
from aa_residue import AA, UnnaturalAA, CheckUnnaturalAA

font = {'family': 'Arial'}
mpl.rc('font', **font)

msg = '''  > {0}
      -in  [ PDB file for Ramachandran density plot ]
      -png [ Output PNG name ]\n
    optional:
      -int [ Resolution (def: 2-deg interval) ]
      -ref [ Directory to Density data for reference Ramachandran distribution ]
      -dpi [ PNG figure DPI resolution (def: 300) ]\n
'''.format(sys.argv[0])
#if len(sys.argv) < 3 or len(sys.argv) > 4: sys.exit(msg)

ref_df = pd.DataFrame(
  { 'Pro':{'file': 'pref_proline.data.bz2', 'bounds': [0,0.002, 0.02,1],
           'cmap': mpl.colors.ListedColormap(['#FFFFFF','#D0FFC5','#7EE589'])},
    'Gly':{'file': 'pref_glycine.data.bz2', 'bounds': [0,0.002, 0.02,1],
           'cmap': mpl.colors.ListedColormap(['#FFFFFF','#FFE8C5','#FED479'])},
    'Gen':{'file': 'pref_general.data.bz2', 'bounds': [0,0.0005,0.02,1],
           'cmap': mpl.colors.ListedColormap(['#FFFFFF','#B3E8FF','#26C3FF'])},
    'PreP':{'file': 'pref_preproline.data.bz2', 'bounds': [0,0.002, 0.02,1],
            'cmap': mpl.colors.ListedColormap(['#FFFFFF','#FFD7DA','#F3ABB0'])}
  } )

############################################################################

def main( ):

  args = UserInput()
  if args.actual_start is None:
    args.actual_start = 1
  else: args.actual_start = int(args.actual_start)
  if args.interval is None:   # default degree interval
    args.interval = 2
  else: args.interval = float(args.interval)
  if args.dpi is None:        # figure DPI
    args.dpi      = 300
  else: args.dpi = int(args.dpi)

  # if reference rama density is available, generate reference figure settings
  ref_dict = {}
  for key in ref_df:
    ref_dict[key] = RefRamaData( args.ref_dir, ref_df[key] )

  # extract input residue dihedral angles and generate figure settings
  res_dict = InputRamaData( args.in_file, args.actual_start, args.interval )

  # generate figure
  GeneratePNG( res_dict, ref_dict, args.png_name, args.dpi )


############################################################################
############################################################################
## extract input residue dihedral angles
def InputRamaData( pdb_file, actual_start, interval ):

  Res_Dih = []
  ## extract sequence and dihedral angle data from PDB structure
  pdb_id = pdb_file.split('.pdb')[0].split('/')[-1]

  # going into each chain in model(s), then each of broken pieces in each
  # chain, get peptide and info
  for model in p.get_structure(pdb_id, file_handle(pdb_file) ):
    for chain in model:
      chain_id = chain.get_id()
      polypeptides = PPBuilder().build_peptides(chain, aa_only=False)

      for poly in polypeptides:
        phi_psi = poly.get_phi_psi_list()

        for res_idx, residue in enumerate(poly):
          resname = residue.get_resname()
          resid   = residue.get_id()[1]
          phi, psi = phi_psi[res_idx]

          # convert 3-char AA to 1-char, check if AA is non-standard AA        
          if CheckUnnaturalAA(resname):
            resname = UnnaturalAA(resname)
          one_resname = AA(resname)
        
          Res_Dih.append([one_resname, resid, chain_id, phi, psi])


  # build dataframe of dihedral angle data
  pdb_df = pd.DataFrame(Res_Dih, 
              columns=['resname','resid','chain','PHI','PSI']).dropna()

  # convert angle from radian to degree
  pdb_df['phi'] = radian2deg(pdb_df.PHI.to_numpy())
  pdb_df['psi'] = radian2deg(pdb_df.PSI.to_numpy())

  # select subsets of residues: Pro, Gly, PrePro, General
  pro_df = pdb_df[pdb_df.resname == 'P']
  gly_df = pdb_df[pdb_df.resname == 'G']
  pp_df  = pdb_df.loc[(pro_df.index-1)]  # pre-Proline residues
  pp_df  = pp_df[pp_df.resname != 'P']
  not_x  = list(pro_df.index) + list(gly_df.index) + list(pp_df.index)
  gen_df = pdb_df.drop(not_x)             # general residues

  res_dict = {'Pro': pro_df, 'Gly': gly_df, 'PreP': pp_df, 'Gen': gen_df}
  
  return res_dict


############################################################################
# if reference rama density is available, generate settings
def RefRamaData( ref_dir, ref_inf ):

  rama_data = pd.read_csv(ref_dir+'/'+ref_inf.file, delimiter=' ', comment='#')
  rama_ref  = rama_data.pivot(index='phi',columns='psi',values='density')

  # data is transpose to get correct orientation
  ref_obj = PNGData( histo2d=rama_ref.transpose() )

  # unique setting for density data to generate correct plot axis order
  # extent is different from res_obj.extent
  # color is (white, light cyan, cyan) at specific contour level
  ref_obj.extent = (-181,181,-181,181)  # different from res_obj.extent
  ref_obj.colors = ref_inf.cmap
  ref_obj.norm   = mpl.colors.BoundaryNorm(ref_inf.bounds, ref_obj.colors.N)

  return ref_obj


############################################################################
## Generate Ramachandran heat map
def GeneratePNG( res_dict, ref_dict, png_name, dpi ):

  for idx, (key, ref_obj) in enumerate(sorted(ref_dict.items(), key=lambda x:x[0].lower())):
    
    plt.figure(2, figsize=(8,8))
    plt.subplot(2,2, idx+1)
    plt.title(key, fontsize=16)

    # reference dihedral angle density as background
    plt.imshow( ref_obj.histo2d[::-1], cmap=ref_obj.colors, 
                norm=ref_obj.norm,     extent=ref_obj.extent )

#    plt.grid(linestyle='--')

    # PDB AA backbone dihedral angle distribution
    plt.scatter( res_dict[key].phi, res_dict[key].psi,
                 alpha=0.67, s=8 )

    ## add additional items
    plt.xlim([-180,180])
    plt.ylim([-180,180])
    plt.plot([-180, 180], [0, 0], color="black", linewidth=1)
    plt.plot([0, 0], [-180, 180], color="black", linewidth=1)
    plt.xticks(np.arange(-180,210, step=60), fontsize=14)
    plt.yticks(np.arange(-180,210, step=60), fontsize=14)

    plt.xlabel(r'Phi $\phi$', fontsize=16)
    plt.ylabel(r'Psi $\psi$', fontsize=16)

  plt.tight_layout()
  plt.savefig(png_name, bbox_inches=0, dpi=dpi)


############################################################################

class PNGData(object):
  def __init__( self, histo2d='', **kwargs):
    self.histo2d = histo2d


# convert radian into degree
def radian2deg(rad):
  return np.rad2deg(rad)


def file_handle(file_name):
  if re.search(r'.gz$', file_name):
    handle = gzip.open(file_name, 'rt')
  elif re.search(r'.bz2$', file_name):
    handle = bz2.open(file_name, 'rt')
  else:
    handle = open(file_name)

  return handle

##########################################################################
def UserInput():
  p = ArgumentParser(description='Command Line Arguments')

  p.add_argument('-in', dest='in_file', required=True,
                 help='PDB structure')
  p.add_argument('-png', dest='png_name', required=True,
                 help='Output PNG filename')

  p.add_argument('-stt', dest='actual_start', required=False)

  p.add_argument('-int', dest='interval', required=False,
                 help='Ramachandran plot resolution (def: 2-degree)')
  p.add_argument('-ref', dest='ref_dir', required=False,
                 help='Directory of Density data for reference Ramachandran distribution')
  p.add_argument('-dpi', dest='dpi', required=False,
                 help='Figure DPI (def: 300)')
  
  return p.parse_args()


##########################################################################
if __name__ == '__main__':
  main()
