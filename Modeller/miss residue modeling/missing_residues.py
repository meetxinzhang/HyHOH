from modeller import *
from modeller.automodel import *    # Load the automodel class

log.verbose()
env = environ()

# directories for input atom files
env.io.atom_files_directory = ['.', './']

class MyModel(automodel):
    def select_atoms(self):
        return selection(self.residue_range('1', '5'),
                         self.residue_range('86', '100'),
                         self.residue_range('291', '291'),
                         self.residue_range('401', '404'),
                         self.residue_range('440', '453'),
                         self.residue_range('565', '582'))

a = MyModel(env, alnfile = '6ws6-fix.ali',
            knowns = '6ws6', sequence = 'rcsb_pdb_6WS6',
            assess_methods=(assess.DOPE, assess.normalized_dope, assess.GA341))
a.starting_model = 1
a.ending_model = 2

a.make()