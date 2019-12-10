#!/usr/bin/env python3

import sys,os,re
import pandas as pd

if int(pd.__version__.split('.')[1]) < 24:
  sys.exit('\n  ERROR: Require Pandas 0.24+ \n')

from argparse import ArgumentParser

##################################################################
def main( ):
  args = UserInput()

  delim = '\s+'
  if args.delim:
    delim = args.delim
  if args.sep:
    sep = args.sep
  else:
    sep = '\t'

  f1_df = pd.read_csv(args.file_1, delimiter=delim)
  f2_df = pd.read_csv(args.file_2, delimiter=delim)

  if len(f1_df.index) != len(f2_df.index):
    sys.exit('\n  ERROR: file_1 {} and file_2 {} do not match\n'.format(
              len(f1_df.index), len(f2_df.index)) )

  final_df = f1_df
  for col in args.cols:
    final_df = pd.concat( [final_df, f2_df.iloc[:,int(col)-1]] , axis=1)

  print(final_df.to_csv(sep=sep, index=False))

#################################################################
def UserInput():
  p = ArgumentParser(description='Command Line Arguments')

  p.add_argument('-in1', dest='file_1', required=True,
                 help='file_1')
  p.add_argument('-in2', dest='file_2', required=True,
                 help='file_2')
  p.add_argument('-cols', dest='cols', required=True, nargs='+',
                 help='Column(s) from file_2 to graft onto file_1 (multiple ok)')
  p.add_argument('-delim', dest='delim', required=False,
                 help='Any Delimiter (def: "\s+"')
  p.add_argument('-sep', dest='sep', required=False,
                 help='Delimiter for output file (def: "\t")')
  return p.parse_args()

##################################################################
if __name__ == '__main__':
  main( )

