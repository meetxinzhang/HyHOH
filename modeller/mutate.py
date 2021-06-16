# encoding: utf-8
"""
@author: Xin Zhang
@contact: zhangxin@szbl.ac.cn
@file: mutate.py
@time: 4/7/21 3:40 PM
@desc:
"""
# This will read a PDB file, change its sequence a little, build new
# coordinates for any of the additional atoms using only the internal
# geometry, and write the mutant PDB file.  It can be seen as primitive
# but rapid comparative modeling for substitution mutants. For more
# sophisticated modeling, see http://salilab.org/modeller/wiki/Mutate%20model

from modeller import *

# -------------- step 1: read pdb ------------------------------
env = environ()
env.io.atom_files_directory = ['../outputs']
outputs_dir = '../outputs/'

in_name = '7c8d_fixed_L452R'
out_name = in_name + '_' + 'E484Q'
index = 748
target = 'GLN'

# Read the topology library with non-hydrogen atoms only:
env.libs.topology.read(file='$(LIB)/top_heav.lib')
# To produce a mutant with all hydrogens, uncomment this line:
# env.libs.topology.read(file='$(LIB)/top_allh.lib')
# Read the CHARMM parameter library:
env.libs.parameters.read(file='$(LIB)/par.lib')

# Read the original PDB file and copy its sequence to the alignment array:
aln = alignment(env)
mdl = model(env, file=in_name+'.pdb')
print("Mapping from residue indices to PDB residue and chain names:")
for r in mdl.residues:
    print("%6d   %3s:%s   %s" % (r.index, r.num, r.chain.name, r.pdb_name))

# -------------- step 2: performs mutants -----------------------------
aln.append_model(mdl, atom_files=in_name+'.pdb', align_codes=in_name)

# Select the residues to be mutated: in this case all ASP residues:
# sel = selection(mdl).only_residue_types('ASP')

# The second example is commented out; it selects residues '1' and '10'.
print('--------', mdl.residues[index])
sel = selection(mdl.residues[index])
# Mutate the selected residues into HIS residues (neutral HIS):
sel.mutate(residue_type=target)

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
mdl.write(file=outputs_dir+out_name+'.pdb')
