# Automao -e 13 \n | e some time-consuming operations
# 2020/12/17
# usage:
# >chmod a+x script.sh
# >./script.sh

################### Pre Process ##################
# cd ACE2
# # grep -v HOH 7jw0.BL00010001.pdb > clean.pdb 

# # water: spce tip3p, -ignh 
# echo -e 2 \n | gmx pdb2gmx -f renum.pdb -o processed.gro -water tip3p
# # force field, count of atom*2 due to H

# # gmx editconf -f processed.gro -o newbox.gro -center 11 10 35 -box 20 20 45
# # gmx editconf -f processed.gro -o newbox.gro -center 14 12.5 15 -box 24 24 45
# gmx editconf -f processed.gro -o newbox.gro -c -d 1.5 -bt cubic

# gmx solvate -cp newbox.gro -cs spc216.gro -o solv.gro -p topol.top

# gmx grompp -f ../../mdp/ions.mdp -c solv.gro -p topol.top -o ions.tpr -maxwarn 1

# echo -e 13 \n | gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral -conc 0.15
# # 13 for SOL

##################  EM  #######################
# gmx make_ndx -f solv_ions.gro
# chose atoms
cd interaction/CNGA1_PCG

gmx grompp -f ../../mdp/em.mdp -c solv_ions.gro -p topol.top -o em.tpr 

gmx mdrun -v -deffnm em

gmx grompp -f ../../mdp/em_2.mdp -c em.gro -p topol.top -o em_2.tpr 

gmx mdrun -v -deffnm em_2

# ################## Balance ####################
gmx grompp -f ../../mdp/nvt.mdp -c em_2.gro -r em_2.gro -p topol.top -o nvt.tpr 

gmx mdrun -deffnm nvt

gmx grompp -f ../../mdp/npt.mdp -c nvt.gro -r nvt.gro -p topol.top -o npt.tpr

gmx mdrun -deffnm npt

################## MD ########################
gmx grompp -f ../../mdp/prod.mdp -c npt.gro -t npt.cpt -p topol.top -o md_0_1.tpr

gmx mdrun -deffnm md_0_1

##################### umbrella ######################
# ./umbrella.sh

############### analysis  #####################
# mkdir md_out
# gmx energy -f nvt.edr -o md_out/temperature.xvg
# # 16 0

# gmx energy -f npt.edr -o md_out/pressure.xvg
# # 18 0

# gmx energy -f npt.edr -o md_out/density.xvg 
# # 24 0

# gmx trjconv -s md_0_1.tpr -f md_0_1.xtc -o md_0_1_noPBC.xtc -pbc mol -ur compact
# # 0 for system

# gmx rms -s md_0_1.tpr -f md_0_1_noPBC.xtc -o md_out/rmsd.xvg -tu ns
# # 4 for backbone

# gmx gyrate -s md_0_1.tpr -f md_0_1_noPBC.xtc -o md_out/gyrate.xvg
# # 1 for protein

# gmx_MMPBSA -O -i ../../mmpbsa.in -cs md_0_1.tpr -ci index.ndx -cg 18 19 -ct md_0_1.xtc

####################################
# sudo shutdown -h +30
