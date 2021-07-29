
mkdir analysis
cd analysis

################ nm  
gmx energy -f nvt.edr -o md_out/temperature.xvg
# 16 0

gmx energy -f npt.edr -o md_out/pressure.xvg
# 18 0

gmx energy -f npt.edr -o md_out/density.xvg 
# # 24 0

gmx gyrate -s md_0.tpr -f md_0_noPBC.xtc -o md_out/gyrate.xvg
# 1 for protein

################# trajectoty
gmx trjconv -s ../npt.gro -f ../md_0.xtc -o md_0_noPBC.xtc -pbc mol #-ur compact
# 0 for system


################# RMSD
gmx rms -s ../npt.gro -f md_0_noPBC.xtc -o rmsd.xvg #-tu ns
# 4 for backbone

xmgrace -nxy rmsd.xvg


################# RMSD refers to average.pdb 

gmx rmsf -f md_0_noPBC.xtc -s ../npt.gro -o rmsf-perdue.xvg -ox average.pdb  -b 1000 -res
# Selects backbone for root mean qsuare calculatuion.

gmx rms -f md_0_noPBC.xtc -s average.pdb -o rmsd-all-atom-vs-avg.xvg
# Selects backbone 

xmgrace -nxy rmsd-all-atom-vs-avg.xvg

python ../../../../xvgReorder.py rmsd-all-atom-vs-avg.xvg

xmgrace -nxy rmsd-all-atom-vs-avg-reorder.xvg

gmx trjconv -f md_0_noPBC.xtc -o cloestAve.xtc -drop rmsd-all-atom-vs-avg.xvg -dropunder 2.80 -dropover 2.82


################# mostFre by Bio3D
R

library(bio3d)

dcdfile <- "md_0_noPBC.dcd"
pdbfile <- "npt.pdb"

dcd <- io.dcd(dcdfile)
pdb <- io.pdb(pdbfile)

# ca.inds <- atom.select(pdb, elety="CA")
ca.inds <- atom.select(pdb, "backbone")

xyz <- fit.xyz(fixed=pdb$xyz, mobile=dcd, fixed.inds=ca.inds$xyz, mobile.inds=ca.inds$xyz)

dim(xyz) == dim(dcd)
> [1]  TRUE TRUE 

# RMSD 
rd <- rmsd(xyz[1,ca.inds$xyz], xyz[,ca.inds$xyz])

plot(rd, typ="l", ylab="RMSD (A)", xlab="Frame No.")
points(lowess(rd), typ="l", col="red", lty=2, lwd=2)

hist(rd, breaks=40, freq=TRUE, main="RMSD Histogram", xlab="RMSD")
lines(density(rd), col="gray", lwd=3)

gmx trjconv -f md_0_noPBC.xtc -o mostFre.xtc -drop rmsd.xvg -dropunder 0.38 -dropover 0.40


################ MMPBSA
./../../../../gmx_mmpbsa_ed.bsh -f ../md_1.xtc -s ../md_1.tpr -n ../index.ndx -com Protein -pro receptor -lig ligand -cou dh -ts ie -b 1000 -e 2000 -i 100
./../../../../gmx_mmpbsa_ed.bsh -f ../md_0.xtc -s ../md_0.tpr -n ../index.ndx -com Protein -pro receptor -lig ligand -cou dh -ts ie -b 1000 -e 10000 -i 500
# mostFre
./../../../../gmx_mmpbsa_ed.bsh -f ../analysis/mostFre.xtc -s ../md_0.tpr -n ../index.ndx -com Protein -pro receptor -lig ligand -cou dh -ts ie -b 1000 -e 10000 -i 200
# cloestAve
./../../../../gmx_mmpbsa_ed.bsh -f ../analysis/cloestAve.xtc -s ../md_0.tpr -n ../index.ndx -com Protein -pro receptor -lig ligand -cou dh -ts ie -b 1000 -e 10000 -i 50

python ../../plotmmpb.py 


################ schlitter
gmx covar -s ../md_0.tpr -f ../analysis/md_0_noPBC.xtc -o eigenvalues.xvg -v eigenvectors.trr
# backbone

gmx anaeig -s ../md_0.tpr -f ../analysis/md_0_noPBC.xtc -v eigenvectors.trr -eig eigenvalues.xvg -entropy yes