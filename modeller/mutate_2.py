# encoding: utf-8
"""
@author: Xin Zhang
@contact: zhangxin@szbl.ac.cn
@file: mutate_2.py
@time: 4/7/21 5:05 PM
@desc:
"""
import sys
import os
from modeller import *
from modeller.optimizers import molecular_dynamics, conjugate_gradients
from modeller.automodel import autosched


def optimize(atmsel, sched):
    for step in sched:
        step.optimize(atmsel, max_iterations=200, min_atom_shift=0.001)
        refine(atmsel)
        cg = conjugate_gradients()
        cg.optimize(atmsel, max_iterations=200, min_atom_shift=0.001)


def refine(atmsel):
    md = molecular_dynamics(cap_atom_shift=0.39, md_time_step=4.0, md_return='FINAL')
    init_vel = True
    for (its, equil, temps) in (
    (200, 20, (150.0, 250.0, 400.0, 700.0, 1000.0)), (200, 600, (1000.0, 800.0, 600.0, 500.0, 400.0, 300.0))):
        for temp in temps:
            md.optimize(atmsel, init_velocities=init_vel, temperature=temp, max_iterations=its, equilibrate=equil)
            init_vel = False


def make_restraints(mdl1, aln):
    rsr = mdl1.restraints
    rsr.clear()
    s = selection(mdl1)
    for typ in ('stereo', 'phi-psi_binormal'):
        rsr.make(s, restraint_type=typ, aln=aln, spline_on_site=True)
    for typ in ('omega', 'chi1', 'chi2', 'chi3', 'chi4'):
        rsr.make(s, restraint_type=typ + '_dihedral', spline_range=4.0, spline_dx=0.3, spline_min_points=5, aln=aln,
                 spline_on_site=True)


modelname, respos1, restyp1, respos2, restyp2, chain, pdb_filename, out_name = sys.argv[1:]  # change

log.verbose()
env = environ(rand_seed=-49837)
env.io.hetatm = True
env.edat.dynamic_sphere = False
env.edat.dynamic_lennard = True
env.edat.contact_shell = 4.0
env.edat.update_dynamic = 0.39
env.libs.topology.read(file='$(LIB)/top_heav.lib')
env.libs.parameters.read(file='$(LIB)/par.lib')

mdl1 = model(env, file=pdb_filename)
ali = alignment(env)
ali.append_model(mdl1, atom_files=pdb_filename, align_codes=modelname)
# ResiduePositionMutationSelection
s = selection(mdl1.chains[chain].residues[respos1], mdl1.chains[chain].residues[respos2])  # change
s.mutate(residue_type=restyp1)  # change
s.mutate(residue_type=restyp2)  # change
ali.append_model(mdl1, align_codes=modelname)
mdl1.clear_topology()
mdl1.generate_topology(ali[-1])
mdl1.transfer_xyz(ali)
mdl1.build(initialize_xyz=False, build_method='INTERNAL_COORDINATES')
mdl2 = model(env, file=pdb_filename)
mdl1.res_num_from(mdl2, ali)
# WriteAndReadMutation
mdl1.write(file=modelname + restyp1 + respos1 + restyp2 + respos2 + chain + '.tmp')  # change
mdl1.read(file=modelname + restyp1 + respos1 + restyp2 + respos2 + chain + '.tmp')  # change
make_restraints(mdl1, ali)
mdl1.env.edat.nonbonded_sel_atoms = 1
sched = autosched.loop.make_for_model(mdl1)
# MutationOptimization
s = selection(mdl1.protein_atoms['CA:' + respos1 + ':' + chain].select_sphere(5),
              mdl1.protein_atoms['CA:' + respos2 + ':' + chain].select_sphere(5))  # change
mdl1.restraints.unpick_all()
mdl1.restraints.pick(s)
s.energy()
s.randomize_xyz(deviation=4.0)
mdl1.env.edat.nonbonded_sel_atoms = 2
optimize(s, sched)
mdl1.env.edat.nonbonded_sel_atoms = 1
optimize(s, sched)
s.energy()
mdl1.write(file=out_name)
os.remove(modelname + restyp1 + respos1 + restyp2 + respos2 + chain + '.tmp')
