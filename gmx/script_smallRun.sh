# 2022/4/13
# usage:
# >chmod a+x script.sh
# >./script.sh path gpu_id 
# for example:
# ./script_smallRun.sh  /home/wurp/workspace/antibody/variant/Omicron/Beta-55   1


#环境参数配置
mdp_dir="/media/xin/WinData/ACS/github/BioUtil/gmx/mdp"
gpu_id=$2  #从输入的第二个参数获取需要运行的 gpu_id


#从第一个参数获取抗体所在目录
thispath=$1
cd $thispath
mkdir -p MD
cd MD


# 循环跑5次
for i in $(seq 1 1)
do 
mkdir -p $i
cd $i

echo -e 2 \n | gmx pdb2gmx -f $thispath/renum.pdb -o processed.gro -water tip3p

gmx editconf -f processed.gro -o newbox.gro -c -d 1.5 -bt cubic


#################### EM 0 ######################
gmx grompp -f $mdp_dir/em_0.mdp -c newbox.gro -p topol.top -o em_0.tpr -maxwarn 1

gmx mdrun -v -deffnm em_0 

################## solvate ####################
gmx solvate -cp em_0.gro -cs spc216.gro -o solv.gro -p topol.top

gmx grompp -f $mdp_dir/ions.mdp -c solv.gro -p topol.top -o ions.tpr -maxwarn 1

echo -e 13 \n | gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral -conc 0.15
# 13 for SOL


#################  EM 1-3 ####################
gmx grompp -f $mdp_dir/em_1.mdp -c solv_ions.gro -p topol.top -o em_1.tpr -r solv_ions.gro

gmx mdrun -v -deffnm em_1 

gmx grompp -f $mdp_dir/em_2.mdp -c em_1.gro -p topol.top -o em_2.tpr 

gmx mdrun -v -deffnm em_2 

gmx grompp -f $mdp_dir/em_3.mdp -c em_2.gro -p topol.top -o em_3.tpr 

gmx mdrun -v -deffnm em_3 


################# Balance ####################
gmx grompp -f $mdp_dir/nvt.mdp -c em_3.gro -r em_3.gro -p topol.top -o nvt.tpr 

gmx mdrun -deffnm nvt -pin on -ntmpi 1 -ntomp 6 -gpu_id $gpu_id -pme gpu -update gpu -bonded gpu

gmx grompp -f $mdp_dir/npt.mdp -c nvt.gro -r nvt.gro -p topol.top -o npt.tpr

gmx mdrun -deffnm npt  -pin on -ntmpi 1 -ntomp 6 -gpu_id $gpu_id -pme gpu -update gpu -bonded gpu


################## MD 不带位置限制 ########################

gmx grompp -f $mdp_dir/prod.mdp -c npt.gro -t npt.cpt -p topol.top -o md_0.tpr

gmx mdrun -deffnm md_0 -gpu_id $gpu_id -pme gpu -update gpu -bonded gpu -pin on -npme 1 -ntmpi 3 -ntomp 12


cd .. # 返回上级目录
done