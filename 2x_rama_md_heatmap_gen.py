#!/usr/bin/env python3   

#############################################################################
##
##  Peter MU Ung @ MSSM / Yale
##
##  v2.0 19.12.03 - generate figure with Seaborn/Matplotlib directly
##
##  use to convert dihedral angle of a residue over time (snapshots) to 
##  HISTOGRAM format and generate figure from here directly
##
##  Required input file format (for each frame): <time> <phi> <psi>
##
##  Use with:
##    -) combine_column.pl
##    -) ecoDnaK.fasta
##    1) rama_md_heatmap.csh
##
############################################################################

import sys, os
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

from scipy import ndimage
from argparse import ArgumentParser

font = {'family': 'Arial'}
mpl.rc('font', **font)

msg = '''  > {0}
      -in  < >     [ phi-psi file for Ramachandran density plot ]
      -img < >     [ output image name (format: png,svg,eps,ps,pdf ]\n
    optional:
      -res < >     [ 1-char AminoAcid code for ref. selection (def: Pro, PreP, Gly, Gen) ]
      -int < >     [ Ramachandran plot resolution, in angular degree (def: 2l) ]
      -ref < >     [ Density data for reference Ramachandran distribution ]
      -fraction < >[ Cutoff fraction of the maximum Histogram value (def: 33) ]
      -smooth < >  [ Histogram data smoothening coefficient (def: 1.15) ]
      -t_step < >  [ Colorbar tick spacing per histogram digits unit (def: 4) ]
      -c_step < >  [ Histo Contour spacing per histogram digits unit (def: 4) ]
      -dpi    < >  [ Figure DPI resolution (def: 300) ]\n
'''.format(sys.argv[0])
if len(sys.argv) < 3: sys.exit(msg)


ref_df = pd.DataFrame(
  { 'Pro':{ 'file': 'pref_proline.data.bz2', 'bounds': [0,0.002, 0.02,1],
            'cmap': mpl.colors.ListedColormap(['#FFFFFF','#D0FFC5','#7EE589'])},
    'Gly':{ 'file': 'pref_glycine.data.bz2', 'bounds': [0,0.002, 0.02,1],
            'cmap': mpl.colors.ListedColormap(['#FFFFFF','#FFE8C5','#FED479'])},
    'Gen':{ 'file': 'pref_general.data.bz2', 'bounds': [0,0.0005,0.02,1],
            'cmap': mpl.colors.ListedColormap(['#FFFFFF','#B3E8FF','#26C3FF'])},
    'PreP':{'file':'pref_preproline.data.bz2', 'bounds': [0,0.002, 0.02,1],
            'cmap':mpl.colors.ListedColormap(['#FFFFFF','#FFD7DA','#F3ABB0'])}
  } )


############################################################################

def main( ):

  args = UserInput()
  interval = float(args.interval)
  fraction = float(args.fraction)
  smooth = float(args.smooth)
  t_step = float(args.t_step)
  c_step = float(args.c_step)

  dpi = int(args.dpi)

  # if reference rama density is available, generate reference figure settings
  ref_obj = RefRamaData( args.rama_ref, ref_df, args.residue )

  # extract input residue dihedral angles and generate figure settings
  res_obj = InputRamaData(args.in_file, interval, fraction,
                          smooth, t_step, c_step    )

  # generate figure
  GenerateImage( res_obj, ref_obj, args.img_name, dpi )


############################################################################
############################################################################
## extract input residue dihedral angles and generate figure settings
def InputRamaData( in_file, interval, fraction, smooth, t_step, c_step ):

  rama_inp = pd.read_csv(in_file, delimiter='\s+').drop(columns=['#Frame'])

  # Generate Ramachandran data in X-axis and Y-axis
  # Edges is array of [-180, 180] at a certain interval
  rama_x, rama_y = zip(*rama_inp.to_numpy())
  edges  = range(-180,180+interval,interval)

  # Generate normalized 2D histogram (array of array)
  Histo, xedges, yedges = np.histogram2d(rama_x, rama_y, 
                                bins=(edges,edges), density=True)

  # Smoothening the 2D histogram data
  Sigma = [ (max(xedges)-min(xedges)) * smooth/len(xedges),
            (max(yedges)-min(yedges)) * smooth/len(yedges) ]
  smooth_hist = ndimage.filters.gaussian_filter(Histo, sigma=Sigma)

  # get scientif notation, then set the cutoff to a fraction of maximum
  max_nm = float(np.max(smooth_hist))
  powers = int('{:e}'.format(max_nm).split('e')[1])
  digits = np.ceil(float(('{:e}'.format(max_nm)).split('e')[0]))
  h_max  = np.float('{0}e{1}'.format(digits, powers))

  # introduce a cutoff to histogram data
  histo2d = smooth_hist - h_max/fraction

  # side bar tick, maximum = histogram value
  cbar_ticks = np.linspace( 0, h_max, num=(digits*t_step)+1 )

  # Contour levels, set to be 'c_step' of the histo value, default is 4x
  levels = np.linspace( 0, h_max, num=(digits*c_step)+1 )

  # X- and Y-axes min and max, will be stretch to be equal
  extent = [  xedges[0] -1, xedges[-1]+1,
              yedges[-1]+1, yedges[0] -1  ]

  # data is transpose to get correct orientation
  # im_extent and im_colors are dummy value to generate holder blank plot
  res_obj = ImageData( histo2d=histo2d.transpose() )
  res_obj.cbar_ticks = cbar_ticks
  res_obj.levels = levels
  res_obj.extent = extent
  res_obj.edges  = len(edges)
  res_obj.im_extent = (-181,181,-181,181)
  res_obj.im_colors = mpl.colors.ListedColormap(['#FFFFFF'])

  return res_obj


############################################################################
# if reference rama density is available, generate settings
def RefRamaData( rama_ref, ref_df, residue ):

  if rama_ref is None:
    return None

  rama_data = pd.read_csv(rama_ref, delimiter='\s+', comment='#')
  rama_ref  = rama_data.pivot(index='phi',columns='psi',values='density')

  max_nm  = float(np.max(rama_data[['density']]))
  powers  = int('{:e}'.format(max_nm).split('e')[1])
  digits  = np.ceil(float(('{:e}'.format(max_nm)).split('e')[0]))
  ref_max = np.float('{0}e{1}'.format(digits, powers))
  ref_levels = np.linspace( 0, ref_max, num=(digits)+1 )

  # data is transpose to get correct orientation
  ref_obj = ImageData( histo2d=rama_ref.transpose() )

  # unique setting for density data to generate correct plot axis order
  # extent is different from res_obj.extent
  # color is (white, light cyan, cyan) at specific contour level
  ref_obj.levels = ref_levels
  ref_obj.extent = (-181,181,-181,181)  # different from res_obj.extent
  ref_obj.colors = ref_df[residue].cmap
  ref_obj.norm   = mpl.colors.BoundaryNorm(ref_df[residue].bounds, ref_obj.colors.N)

  return ref_obj


############################################################################
## Generate Ramachandran heat map
def GenerateImage( res_obj, ref_obj, img_name, dpi ):

  plt.figure(2, figsize=(7,5.5))
  colors = mpl.cm.jet

  bar_extend  = 'neither'
  plot_extend = 'neither'

  ## if data is available, generate contour map for Reference Ramachandran 
  ##density map of general AA.
  ## imshow forces figure to have axis ratio 1:1
  ## issue with y-axis data ordering, deal with it by inversing y-axis
  if ref_obj:
    print('  ## INFO: Reference Ramachandran density map is used ##')
    plt.imshow( ref_obj.histo2d[::-1], 
                extent=ref_obj.extent, cmap=ref_obj.colors, norm=ref_obj.norm )
  else:
    ## generate a fake imshow with dummy matrix to get axis ratio 1:1
    plt.imshow( np.zeros(shape=(res_obj.edges-1, res_obj.edges-1)),
      extent=res_obj.im_extent, cmap=res_obj.im_colors)

  ## overlay input AA histogram heat map on top of reference map, if exists
  plt.contourf( res_obj.histo2d, 
                origin='upper', extend=plot_extend, alpha=0.6,
                extent=res_obj.extent, levels=res_obj.levels,
                cmap=mpl.cm.get_cmap(colors, len(res_obj.levels)) )

  ## create colorbar instance on side based on last data input
  cbar = plt.colorbar(ticks=res_obj.cbar_ticks, format=('%.1e'),
                      extend=bar_extend, aspect=20 )
  bar_label = '% Population'
  cbar.ax.set_ylabel(bar_label, rotation=270, fontsize=18, labelpad=20)

  ## then overlay contour lines on top of heat map
  plt.contour( res_obj.histo2d, 
              extent=res_obj.extent, levels=res_obj.levels,
              origin='upper', colors='black', linewidths=0.67, alpha=0.4 )

  ## add additional items
  plt.xlim([-180,180])
  plt.ylim([-180,180])
  plt.plot([-180, 180], [0, 0], color="black", linewidth=1)
  plt.plot([0, 0], [-180, 180], color="black", linewidth=1)
  plt.xticks(np.arange(-180,210, step=60), fontsize=14)
  plt.yticks(np.arange(-180,210, step=60), fontsize=14)

  plt.xlabel(r'Phi $\phi$', fontsize=14)
  plt.ylabel(r'Psi $\psi$', fontsize=14)
#  plt.grid(linestyle='--')

  plt.savefig(img_name, bbox_inches=0, dpi=dpi)


############################################################################

class ImageData(object):
  def __init__( self, histo2d='', **kwargs):
    self.histo2d = histo2d


##########################################################################
def UserInput():
  p = ArgumentParser(description='Command Line Arguments')

  p.add_argument('-in', dest='in_file', required=True,
                  help='Phi-Psi time-series file of one residue (accept zipped file)')
  p.add_argument('-img', dest='img_name', required=True,
                  help='Output Image filename, extension as format (e.g. .png,svg,eps,ps,pdf)')

  p.add_argument('-int', dest='interval', required=False, default=2,
                  help='Ramachandran plot resolution, in angular degree (def: 2)')
  p.add_argument('-res', dest='residue', required=False, default='Gen',
                  help='AA type (Gen,Gly,Pro,PreP) for background dihedral density (def: Gen)')
  p.add_argument('-ref', dest='rama_ref', required=False,
                  help='Density data for reference Ramachandran distribution')
  p.add_argument('-smooth', dest='smoothen', required=False, default=1.15,
                  help='Histogram data smoothening (def: 1.15)')
  p.add_argument('-fraction', dest='fraction', required=False, default=33,
                  help='Cutoff fraction of the maximum Histogram value (def: 33)')
  p.add_argument('-t_step', dest='t_step', required=False, default=4,
                  help='Colorbar tick spacing per Histogram digits value (def: 4)')
  p.add_argument('-c_step', dest='c_step', required=False, default=4,
                  help='Histogram Contour spacing per Histogram digits value (def: 4)')
  p.add_argument('-dpi', dest='dpi', required=False, default=300,
                  help='Figure DPI (def: 300)')

  return p.parse_args()


##########################################################################
if __name__ == '__main__':
  main()
