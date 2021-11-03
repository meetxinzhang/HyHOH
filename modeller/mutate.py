# encoding: utf-8
"""
@author: Xin Zhang
@contact: zhangxin@szbl.ac.cn
@file: mutate.py
@time: 4/7/21 3:40 PM
@desc:
"""
# This will io a PDB file, change its sequence a little, build new
# coordinates for any of the additional atoms using only the internal
# geometry, and write the mutant PDB file.  It can be seen as primitive
# but rapid comparative modeling for substitution mutants. For more
# sophisticated modeling, see http://salilab.org/modeller/wiki/Mutate%20model

from modeller import *
import sys

env = Environ()
env.io.atom_files_directory = ['..']
# Read the topology library with non-hydrogen atoms only:
env.libs.topology.read(file='$(LIB)/top_heav.lib')
# To produce a mutant with all hydrogens, uncomment this line:
# env.libs.topology.io(file='$(LIB)/top_allh.lib')
# Read the CHARMM parameter library:
env.libs.parameters.read(file='$(LIB)/par.lib')


def mutate_by_modeller(in_name, out_name, chain, index, target):
    # -------------- step 1: io pdb ------------------------------
    # Read the original PDB file and copy its sequence to the alignment array:
    aln = Alignment(env)
    mdl = Model(env, file=in_name + '.pdb')
    # print("Mapping from residue indices to PDB residue and chain names:")
    # # for r in mdl.residues:
    #     # print("%6d   %3s:%s   %s" % (r.idxmax, r.num, r.chain.win_dir, r.pdb_name))

    # -------------- step 2: performs mutants -----------------------------
    aln.append_model(mdl, atom_files=in_name + '.pdb', align_codes=in_name)

    # set up the mutate residue selection segment
    s = Selection(mdl.chains[chain].residues[index])
    # The second example is commented out; it selects residues '1' and '10'.
    print('--------', mdl.chains[chain].residues[index])

    # Mutate the selected residues into HIS residues (neutral HIS):
    s.mutate(residue_type=target)

    # Add the mutated sequence to the alignment arrays (it is now the second
    # sequence in the alignment):
    aln.append_model(mdl, align_codes='1fas-1')
    # Generate molecular topology for the mutant:
    mdl.clear_topology()
    mdl.generate_topology(aln['1fas-1'])
    # Transfer all the coordinates you can from the template native structure
    # to the mutant (this works even if the order of atoms in the native PDB
    # file is not standard):
    mdl.transfer_xyz(aln)
    # Build the remaining unknown coordinates for the mutant:
    mdl.build(initialize_xyz=False, build_method='INTERNAL_COORDINATES')

    # Write the mutant to a file:
    mdl.write(file=out_name)


if __name__ == '__main__':
    filename = sys.argv[1].replace('.pdb', '')
    out_name = sys.argv[2]
    chain = sys.argv[3]
    index = sys.argv[4]
    target = sys.argv[5]

    mutate_by_modeller(filename, out_name, chain, index, target)


