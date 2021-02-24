from modeller import *

env = environ()
env.io.atom_files_directory = ['./out', '../atom_files']

aln = alignment(env)
mdl = model(env, file='6ws6_origin', model_segment=('FIRST:A','LAST:F'))
aln.append_model(mdl, align_codes='6ws6_origin', atom_files='6ws6_origin.pdb')
aln.append(file='6ws6.ali', align_codes='6ws6')
aln.align2d()
aln.write(file='alignment_6ws6.ali', alignment_format='PIR')
aln.write(file='alignment_6ws6.pap', alignment_format='PAP')
