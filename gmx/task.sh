#!/bin/bash
#PBS -N xins md task
#PBS -q gpu-1-2
#PBS -l nodes=1:ppn=32
#PBS -W x=GRES:gpu@4
#PBS -l walltime=96:00:00
#PBS -V
#PBS -S /bin/bash


#######enviroment
module load intel_mpi/2018_u1
module load gcc/7.5.0
module load gromacs/2019.3-intel_2018_u1-gcc_7.5.0-cuda10.0
module load fftw/3.3.4-intel_2018_u1
# source /lustre/software/gmx/gromacs_2019.3-intel_2018_u1-gcc_7.5.0-cuda10/bin/GMXRC
source /lustre/home/xzhang/env/gmx2020/bin/GMXRC

########Get computing resource variables
echo "process will start at : "
date
echo "++++++++++++++++++++++++++++++++++++++++"
cd $PBS_O_WORKDIR
NP=`cat $PBS_NODEFILE | wc -l`
NN=`cat $PBS_NODEFILE | sort | uniq | tee /tmp/nodes.$$ | wc -l`
cat $PBS_NODEFILE > /tmp/nodefile.$$


#########Command
##echo "++++++++++++++++++++++++++++++++++++++++"
#cd /lustre/home/xzhang/workplace/goto_nssc
#gmx grompp -f mdp/pull.mdp -c npt.gro -p topol.top -n index.ndx -t npt.cpt -o pull.tpr -r npt.gro
#gmx mdrun -s pull.tpr

#mkdir output
#echo -e 0 \n | gmx trjconv -s pull.tpr -f traj_comp.xtc -o output/conf.gro -sep

# perl distances.pl
./script_copy.sh
# ./umbrella.sh

########clean
echo "process end at : "
date
rm -f /tmp/nodefile.$$
rm -f /tmp/nodes.$$
module unload gromacs_2019.3-intel_2018_u1-gcc_7.5.0-cuda10 




