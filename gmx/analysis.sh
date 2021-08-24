
mkdir analysis
cd analysis

############################ nm  
# gmx energy -f nvt.edr -o md_out/temperature.xvg
# # 16 0

# gmx energy -f npt.edr -o md_out/pressure.xvg
# # 18 0

# gmx energy -f npt.edr -o md_out/density.xvg 
# # # 24 0

# gmx gyrate -s md_0.tpr -f md_0_noPBC.xtc -o md_out/gyrate.xvg
# # 1 for protein

# ################# trajectoty
echo -e 0 \n | gmx trjconv -s ../md_0.tpr -f ../md_0.xtc -o md_0_whole.xtc -pbc whole #-ur compact
# 0 for system
echo -e 0 \n | gmx trjconv -s ../md_0.tpr -f md_0_whole.xtc -o md_0_nojump.xtc -pbc nojump
# 0 for system
echo -e 1 \n 0 \n | gmx trjconv -s ../md_0.tpr -f md_0_nojump.xtc -o md_0_mol.xtc -pbc mol -center
# 1 for protein and 0 out for system
echo -e 1 \n 0 \n | gmx trjconv -s ../md_0.tpr -f md_0_mol.xtc -o md_0_fit.xtc -fit rot+trans
# 1 for protein and 0 out for system


############################## RMSD
echo -e 4 \n 4 \n | gmx rms -s ../md_0.tpr -f md_0_fit.xtc -o rmsd.xvg #-tu ns
# 4 for backbone
# xmgrace -nxy rmsd.xvg


############################### mostFre

output=$(python /pygmx/most_freq.py rmsd.xvg)
boundaries=($output)
rd_min=${boundaries[0]}
rd_max=${boundaries[1]}

gmx trjconv -f md_0_fit.xtc -o mostFre.xtc -drop rmsd.xvg -dropunder $rd_min -dropover $rd_max

./gmx_mmpbsa_normal_skip.sh -f ../analysis/mostFre.xtc -s ../md_0.tpr -n ../index.ndx -com Protein -pro receptor -lig ligand -cou dh -ts ie -b 1000 -e 10000 -i 200


# R
# library(bio3d)

# dcdfile <- "md_0_nojump.dcd"
# pdbfile <- "npt.pdb"

# dcd <- read.dcd(dcdfile)
# pdb <- read.pdb(pdbfile)

# # ca.inds <- atom.select(pdb, elety="CA")
# ca.inds <- atom.select(pdb, "backbone")

# xyz <- fit.xyz(fixed=pdb$xyz, mobile=xdcd, fixed.inds=ca.inds$xyz, mobile.inds=ca.inds$xyz)

# dim(xyz) == dim(dcd)
# > [1]  TRUE TRUE 

# # RMSD 
# rd <- rmsd(xyz[1,ca.inds$xyz], xyz[,ca.inds$xyz])

# plot(rd, typ="l", ylab="RMSD (A)", xlab="Frame No.")
# points(lowess(rd), typ="l", col="red", lty=2, lwd=2)

# hist(rd, breaks=40, freq=TRUE, main="RMSD Histogram", xlab="RMSD")
# lines(density(rd), col="gray", lwd=3)

################# RMSD refers to average.pdb 
# gmx rmsf -f md_0_nojump.xtc -s ../npt.gro -o rmsf-perdue.xvg -ox average.pdb  -b 1000 -res
# # Selects backbone for root mean qsuare calculatuion.

# gmx rms -f md_0_nojump.xtc -s average.pdb -o rmsd-all-atom-vs-avg.xvg
# # Selects backbone 

# xmgrace -nxy rmsd-all-atom-vs-avg.xvg

# python ../../../../xvgReorder.py rmsd-all-atom-vs-avg.xvg

# xmgrace -nxy rmsd-all-atom-vs-avg-reorder.xvg

# gmx trjconv -f md_0_nojump.xtc -o cloestAve.xtc -drop rmsd-all-atom-vs-avg.xvg -dropunder 2.80 -dropover 2.82

################ MMPBSA
# $code_dir/gmx_mmpbsa_normal_skip_TdS.sh -f ../analysis/md_0_fit.xtc -s ../md_1.tpr -n ../index.ndx -com Protein -pro receptor -lig ligand -cou dh -ts ie -b 1000 -e 2000 -i 100 -dir ./
# $code_dir/gmx_mmpbsa_normal_skip_TdS.sh -f ../analysis/md_0_fit.xtc -s ../md_0.tpr -n ../index.ndx -com Protein -pro receptor -lig ligand -cou dh -ts ie -b 1000 -e 10000 -i 1000 -dir ./
# # mostFre
# ./../../../../gmx_mmpbsa_ed.bsh -f ../analysis/mostFre.xtc -s ../md_0.tpr -n ../index.ndx -com Protein -pro receptor -lig ligand -cou dh -ts ie -b 1000 -e 10000 -i 200
# # cloestAve
# ./../../../../gmx_mmpbsa_ed.bsh -f ../analysis/cloestAve.xtc -s ../md_0.tpr -n ../index.ndx -com Protein -pro receptor -lig ligand -cou dh -ts ie -b 1000 -e 10000 -i 50

# python ../../plotmmpb.py 


# ################ schlitter
# gmx covar -s ../md_0.tpr -f ../analysis/md_0_nojump.xtc -o eigenvalues.xvg -v eigenvectors.trr
# # backbone

# gmx anaeig -s ../md_0.tpr -f ../analysis/md_0_nojump.xtc -v eigenvectors.trr -eig eigenvalues.xvg -entropy yes


# ################# hyHOH
# python /media/xin/WinData/ACS/github/BioUtil/gromacs/pygmx/short_time_hoh.py /media/xin/Raid0/ACS/gmx/interaction/ding/7KFY/md_0.tpr /media/xin/Raid0/ACS/gmx/interaction/ding/7KFY/analysis/mostFre.xtc 196 632 1 195