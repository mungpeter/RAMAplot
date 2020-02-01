#!/usr/bin/env python3

##########################################################################
##
## Peter M.U. Ung @ MSSM
##
##  18.06.19  v2.0  add chain information
##
##  Take in a PDB structure and generate FASTA sequence from it
##  Can only handle standard amino acids
##
##########################################################################

import sys,glob
from Bio.PDB.PDBParser import PDBParser 
from Bio.PDB.Polypeptide import PPBuilder

msg = """
    > {0}\n      [PDB filename: .pdb] [Name of PDB] [FASTA Output Name]

  Note: Capping Groups (ACE/NME) and nonstandard amino acids are not recognized
        need to add 'X' to result FASTA manually
        """.format(sys.argv[0])
if len(sys.argv) != 4: sys.exit(msg)

def FASTA_Gen(pdb_file, pdb_id, output_file):

  m =PDBParser(PERMISSIVE=1).get_structure(pdb_id, pdb_file)

  w = open(output_file, 'w')
  for chain in m.get_chains():
    residues = list(chain.get_residues())
    chain_id = chain.get_id()

    s = ''
    peptides = PPBuilder().build_peptides(chain, aa_only=False)
    for peptide in peptides:
      s = s + peptide.get_sequence() + '\n\n'

    w.write('>{0}_{1}|{2}-{3}\n'.format(
                  pdb_id, chain_id, 
                  residues[0].get_id()[1], residues[-1].get_id()[1] ) )
    w.write(str(s))

  w.close()

##############################################
if __name__ == "__main__":
  FASTA_Gen(sys.argv[1], sys.argv[2], sys.argv[3])
