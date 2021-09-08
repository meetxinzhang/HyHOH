# Automao -e 13 \n | e some time-consuming operations
# 2020/12/17
# usage:
# >chmod a+x script.sh
# >./script.sh

gpu_id=0
mdp_dir=/media/xin/WinData/ACS/github/BioUtil/gmx

################### Pre Process ##################
# cd interaction/ding/7KFY

# grep -v HOH 7jw0.BL00010001.pdb > clean.pdb 

# water: spce tip3p. select force field 2. count of atom*2 due to H
# echo -e 2 \n | gmx pdb2gmx -f renum.pdb -o processed.gro -water tip3p

"""
# make index 
> echo -e chain H L & 1 \n  gmx make_ndx -f renum.pdb -o index.ndx

# generate weak constrain which used in constraining MD.
1. divide the processed.gro into chain.gro by chains manually
2. run:
echo -e 4 \n | gmx genrestr -fc 100 100 100 -f renum_X.pdb -o posre100_chain_X

# 1. Copy posre_Protein_chain_X.itp and rename > posre100_chain_X.
# 2. Modify all force values from 1000 to 100 carefully, pay attention to 1000th atom if you use [change all] by text editor like VScode.

3. Add following code at last line of each topol_Protein_chain_X.itp respectively  
; Include Prodution Position restraint file 
#ifdef POSRES100 
#include "posre100_chain_X.itp" 
#endif 

4. Modify following para at first line of interesting prod_posres.mdp  
define = -DPOSRES100  
"""

# gmx editconf -f processed.gro -o newbox.gro -center 14 12.5 15 -box 24 24 45
gmx editconf -f processed.gro -o newbox.gro -c -d 1.5 -bt cubic


#################### EM 0 ######################
gmx grompp -f $mdp_dir/mdp/em_0.mdp -c newbox.gro -p topol.top -o em_0.tpr -maxwarn 1

gmx mdrun -v -deffnm em_0


# ################## solvate ####################
# gmx solvate -cp newbox.gro -cs spc216.gro -o solv.gro -p topol.top
gmx solvate -cp em_0.gro -cs spc216.gro -o solv.gro -p topol.top

gmx grompp -f $mdp_dir/mdp/ions.mdp -c solv.gro -p topol.top -o ions.tpr -maxwarn 1

echo -e 13 \n | gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral -conc 0.15
# # 13 for SOL


#################  EM 1-3 ####################
gmx grompp -f $mdp_dir/mdp/em_1.mdp -c solv_ions.gro -p topol.top -o em_1.tpr -r solv_ions.gro

gmx mdrun -v -deffnm em_1

gmx grompp -f $mdp_dir/mdp/em_2.mdp -c em_1.gro -p topol.top -o em_2.tpr 

gmx mdrun -v -deffnm em_2

gmx grompp -f $mdp_dir/mdp/em_3.mdp -c em_2.gro -p topol.top -o em_3.tpr 

gmx mdrun -v -deffnm em_3


################# Balance ####################
gmx grompp -f $mdp_dir/mdp/nvt.mdp -c em_3.gro -r em_3.gro -p topol.top -o nvt.tpr

gmx mdrun -deffnm nvt -update gpu -gpu_id $gpu_id

gmx grompp -f $mdp_dir/mdp/npt.mdp -c nvt.gro -r nvt.gro -p topol.top -o npt.tpr
gmx mdrun -deffnm npt -update gpu -gpu_id $gpu_id


################ MD with posres ###############

# gmx grompp -f /media/xin/WinData/ACS/github/BioUtil/gmx/mdp/prod_posres.mdp -c npt.gro -t npt.cpt -p topol.top -o md_1.tpr -r npt.gro

# gmx mdrun -deffnm md_1 -gpu_id $gpu_id -pme gpu -update gpu -bonded gpu -npme 1 -pin on -ntmpi 3 -ntomp 6

################## MD ########################
gmx grompp -f $mdp_dir/mdp/prod.mdp -c npt.gro -t npt.cpt -p topol.top -o md_0.tpr

gmx mdrun -deffnm md_0 -gpu_id $gpu_id -pme gpu -update gpu -bonded gpu #-pin on -ntmpi 2 -ntomp 6 -npme 1

################ extend MD #################
# gmx convert-tpr -s md_0.tpr -extend 13000 -o md_0.tpr

# gmx mdrun -v -deffnm md_0 -cpi md_0.cpt -pin on -ntmpi 3 -ntomp 12 -gpu_id 012 -pme gpu -npme 1 -update gpu -bonded gpu

##################### umbrella ######################
# ./umbrella.sh

# sudo shutdown -h +30
