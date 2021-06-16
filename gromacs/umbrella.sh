#!/bin/bash

cd ACE2_pull

# gmx make_ndx -f npt.gro
# chose atoms

########################### pull #############################################
gmx grompp -f ../mdp/pull.mdp -c md_0_1.gro -p topol.top -n index.ndx -t md_0_1.cpt -o pull.tpr -r md_0_1.gro
gmx mdrun -s pull.tpr

mkdir pull_out
echo -e 0 \n | gmx trjconv -s pull.tpr -f traj_comp.xtc -o pull_out/conf.gro -sep

cat summary_distances.dat
cd pull_out

count=$(ls | wc -l)
let "count=count-1"
for gro in {0..500}
do
    gmx distance -s ../pull.tpr -f conf$gro.gro -n ../index.ndx -oav dist$gro.xvg -select 'com of group "receptor" plus com of group "ligand"' >/dev/null 2>&1
    
    cat dist$gro.xvg | while read line
    do  
        if [[ ${line:0:1} != "#" ]]&&[[ ${line:0:1} != "@" ]];then
            array=($line)
            echo "$gro   ${array[1]} " >> ../summary_distances.dat
        fi
    done
    rm dist$gro.xvg  # delete 
    echo "processed conf$gro.gro"
done


# perl ../distances.pl

# ######################### for ############################################
# confs=(0 22 41 54 74 84 88 101 105 111 118 122 130 138 148 155 164 170 179 184 187 192 196 202 207 212 215 220 226 229 236 242 245 250 257 262 266 273 275)
# # confs=(0 165 203 242 281 312 339 358 387 411 435 455 476 507 535 561 589 627 664 681 707 730 756 791)

# # mkdir um_em_out
# mkdir um_npt_out
# mkdir um_md_out
# mkdir analysis

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

#     # gmx grompp -f ../mdp/em_umbrella.mdp -c pull_out/conf$i.gro -p topol.top -n index.ndx -o um_em_out/em$i.tpr -r pull_out/conf$i.gro
#     # gmx mdrun -v -deffnm um_em_out/em$i

#     gmx grompp -f ../mdp/npt_umbrella.mdp -c pull_out/conf$i.gro -p topol.top -n index.ndx -o um_npt_out/npt$i.tpr -r pull_out/conf$i.gro
#     gmx mdrun -deffnm um_npt_out/npt$i

#     gmx grompp -f ../mdp/prod_umbrella.mdp -c um_npt_out/npt$i.gro -p topol.top -n index.ndx -o um_md_out/umbrella$i.tpr -r um_npt_out/npt$i.gro
#     gmx mdrun -deffnm um_md_out/umbrella$i -pf analysis/pull-umbrella$i.xvg -px analysis/pullx-umbrella$i.xvg
# done

# gmx wham -it tpr-files.dat -if pullf-files.dat -o -hist -unit kCal
