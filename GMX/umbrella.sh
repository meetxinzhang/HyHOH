#!/bin/bash

# gmx make_ndx -f npt.gro
# chose atoms

########################### pull #############################################
gmx grompp -f mdp/pull.mdp -c npt.gro -p topol.top -n index.ndx -t npt.cpt -o pull.tpr -r npt.gro
# or > -r npt.gro
gmx mdrun -s pull.tpr

mkdir pull_out
echo -e 0 \n | gmx trjconv -s pull.tpr -f traj_comp.xtc -o pull_out/conf.gro -sep

perl distances.pl

# ######################### for ############################################
# confs=(0 123 165 183 203 226 242 261 281 298 312 324 339 351 358 375 387 402 411 422 435 445 455 466 476 494 507 521 535 545 561 576 589 612 627 641 664 674 681 692 707 720 730 744 756 769 791)
# confs=(0 165 203 242 281 312 339 358 387 411 435 455 476 507 535 561 589 627 664 681 707 730 756 791)
# confs=(123 183 226 261 298 324 351 375 402 422 445 466 494 521 545 576 612 641 674 692 720 744 769)

# mkdir um_md_out
# mkdir analysis
# mkdir um_npt_out

# cat tpr-files.dat
# cat pullf-files.dat

# for j in ${confs[@]}
# do
#     echo "um_md_out/umbrella$j.tpr" >> tpr-files.dat
#     echo "analysis/pull-umbrella$j.xvg" >> pullf-files.dat
# done

# for i in ${confs[@]}
# do
#     echo "--------------- step: $i -----------------------------"
#     gmx grompp -f mdp/npt_umbrella.mdp -c pull_out/conf$i.gro -p topol.top -n index.ndx -o um_npt_out/npt$i.tpr -r pull_out/conf$i.gro
#     gmx mdrun -deffnm um_npt_out/npt$i

#     gmx grompp -f mdp/md_umbrella.mdp -c um_npt_out/npt$i.gro -p topol.top -n index.ndx -o um_md_out/umbrella$i.tpr -r um_npt_out/npt$i.gro
#     gmx mdrun -deffnm um_md_out/umbrella$i -pf analysis/pull-umbrella$i.xvg -px analysis/pullx-umbrella$i.xvg
# done

# gmx wham -it tpr-files.dat -if pullf-files.dat -o -hist -unit kCal
