**Adding explicit interfacial water molecules improves the calculation of binding free energy for the spike RBD – antibody system**

Xin Zhang, Ruiping Wu, Qinglian Liu, Lei Zhou*

Institute of Molecular Physiology, Shenzhen Bay Laboratory, Shenzhen, Guangdong Province, China

*Correspondence to: Lei Zhou，email: leizhou@szbl.ac.cn; ORCID ID: 0000-0002-3838-1836



# Abstract

​	Antibodies that recognize the spike protein of SARS-COV-2, especially the neutralizing antibodies, carry great hope in the treatment and final elimination of COVID-19. Driven by a synchronized global effort, thousands of antibodies against the spike protein have been identified during the past two years. Among them, the corresponding structural information for hundreds of these antibodies became available and thus provided an excellent opportunity for the theoretical investigation to catch up. We chose MM/PBSA, an efficient approach for estimating the binding free energy between the antibodies and the RBD domain of SARS-COV-2 spike proteins. We discovered that adding explicitly-treated water molecules to the interface between RBD and antibody effectively improves the results. These interfacial water molecules, together with surface and internal water molecules, behave drastically from bulk water and exert peculiar impacts on protein dynamics and energy. Our results demonstrate the importance of these structural water molecules and provide an effective treatment for more accurate binding free energy calculations. 



**Keywords**



# Introduction

​	The coronavirus disease 2019 (COVID-19) is caused by the SARS-CoV-2 virus and has become a pandemic for more than two years. The SARS-CoV-2 virus belongs to the family of positive-sense single-stranded RNA viruses, which includes seasonal coronavirus, SARS-CoV, and Middle East Respiratory Syndrome (MERS)-CoV. SARS-CoV-2 attacks cells in the respiratory system, mainly through a specific interaction between the spike glycoprotein on the viral surface and the corresponding receptor on cell surface, the angiotensin converting enzyme 2 (ACE2). Therefore, much of the ongoing global research effort on SARS-CoV-2 has been devoted to the understanding of spike in membrane fusion and viral entry and the development of strategies to block the interaction between spike with ACE2.

​	The spike glycoprotein of SARS-CoV-2 is composed of three presumably identical protomers, with each containing 1273 amino acids. The full-length S protein includes a bulky extracellular ectodomain, a transmembrane domain (a.a. 1213-1237) and a short cytoplasmic tail. The mature form of S protein is extensively glycosylated, with 22 potential N-glycosites and numerous O-glycosites decorating the ectodomain. A distinctive polybasic protease (furin) cleavage site separates the S protein into S1 (a.a. 1-685 or 667???) and S2 (a.a. 686-1273). In addition, within the S2 domain, the S2’ site can be cleaved by TMPRSS2, a membrane-anchoring serine protease. These two post-translational proteolytic processes play important roles in the life cycle of SARS-CoV-2. In berief, the S1 contributes to the initial contact with ACE2 on host cell surface and S2 is responsible for the following membrane fusion and viral entry.

​	Within the S1 domain, two structural motifs can be identified: the more conserved N-terminal domain (NTD, a.a. 13-305) and the receptor-binding domain (RBD, a.a. 319-541). RBD is in direct contact with ACE2 and therefore is the determining factor during the initial viral recognition of host cells. Driven by a concerted effort by the research community, mechanistic understanding of RBD and the interaction with ACE2 has been advancing rapidly. RBD adopts two dramatic conformations, down (closed) and up (open) state. In the open state, RBD swings upward and through a patch composed of ~25 amino acids interacts with the receptor.

​	Antibodies that target the spike protein, especially RBD, can potentially block the interaction between spike and ACE2 and thereby provide an effective treatment for Covid-19. Neutralizing antibodies (nAbs) against SARS-Cov-2 have been identified from diverse sources including convalescent patients recovering from infections of SARS-Cov-2 during the past two years, saved blood cells collected from individuals (S309, more???) infected with SARS-Cov in 2003, naïve human antibody libraries, and genetically humanized mice (VelocImmune) infected with SARS-CoV-2. Hundreds of effective nAbs have been discovered, with many (???) going through clinical trial and several received FDA approval.  

​	The binding affinity between the antibody and the RBD domain is a critical measure of antibody quality. Experimental approaches, including surface plasmon resonance (SPR), biolayer interferometry (BLI), and isothermal titration calorimetry (ITC), have been used to interrogate protein-protein interaction and obtain information on the binding affinity and knetics1,2. On the other hand, computational estimation of the binding free energy has been a major goal of theoretical approaches. For computational methods of binding free energy, a tradeoff exists between the complexity of computation and the accuracy of the result. Methods based on molecular dynamic simulation produce more accurate results than the docking methods and can be separated into pathway-dependent methods, which are rigorous and computational expensive, such as free energy perturbation (FEP)3 and thermodynamic integration (TI)4, and pathway-independent methods, such as linear interaction energy (LIE)5 and molecular mechanics with Poisson-Boltzmann and surface area solvation (MM/PBSA).

​	 is relatively efficient, but not as precise. In contrast, molecular mechanics/Poisson−Boltzmann surface area (MM/PBSA)6 and molecular mechanics/generalized Born surface area (MM/GBSA)7 methods take into account efficiency and accuracy, so it is being used more and more now. But there are still some problems with it, especially when it comes to highly charged protein system. Yan-jing Sheng et al. developed a new method based on MM/PBSA by adding an exponential damping factor to the Coulombic interaction energy, which is called screening MM/PBSA8. Compared to standard MM/PBSA, screening MM/PBSA has largely improved the Pearson correlation coefficient. Furthermore, Interfacial water has long been thought to play an important role in complex binding. The binding of receptor and ligand is the result of multiple interactions at the binding interface. Energetically, these interactions include hydrogen bonding, hydrophobicity, van der Waals interactions and ion pairing. DG et al. find that the interfacial water contributes around 25% of the total calculated binding strength, suggesting water plays a significantly important role in the stabilization of the complex9. Structurally, water has an important complementarity for the binding of receptors and ligands. With the development of X-ray, cryo-EM, high-resolution 3D structures of the complexes reveal that water is a constituent of protein receptors and ligands, which improves complementarity by filling the gaps between amino acids10. Ahmad et al. reveal that water forms an adhesive hydrogen bond network between the interfaces, which stabilizes the complex formation, by molecular dynamics simulations11. There are also some works that consider water when calculating binding free energy. Irene Maffucci et al. improve the correlation by adding explicit ligand hydration shells in four different systems12. Sergio Wong et al. found that when calculating the binding free energy in MMPBSA, the mutant and wild-type could only be distinguished when the explicit water molecule was included13.

​	Herein we selected more than 20 newly reported RBD-antibody structures from PDB database, all of which are accompanied by experimental measurement of binding affinity. We applied the MM/PBSA method to derive the theoretical binding free energy and aligned the results with experimental measurements. We discovered that including explicitly treated interfacial water between RBD and antibody in MM/PBSA calculation dramatically improved the accuracy of the theoretical results. Application of a modified MM/PBSA method, which adopted an exponential damping factor in the Coulomb electrostatic energy function, further improved the results. With these improvements in methodology, we extended our prediction of binding energy to the RBD of recently emerged SARS-COV-2 variants, including Alpha, Beta, Gamma, Delta, and Omicron variants, and further validated our approach. Finally, from the perspective of energy, structure, and dynamics, we investigated the role of interfacial water molecules in the binding between antibodies and the RBD of SAR-COV-2.





# Materials and Methods

## **1.** Starting structures of RBD-antibod complexes.

​	25 (updated?) RBD-antibody structures were downloaded from the Protein Data Bank (PDB), including 5 complex structures containing RBD from SARS-COV-2 variants. Detailed information for these structures as well as the corresponding experimental binding affinity is shown in Table 1.

**Table 1. PDB ID of the Data Set**

| Gourp no. | system                          | PDB ID                                                       |
| --------- | ------------------------------- | ------------------------------------------------------------ |
| 1         | SARS-COV-2 RBD-binding antibody | 7KFY7KFY, 7KFX7KFY, 7KFV7KFY, 7KFW7KFY, 7JVA18, 7KGK19, 7C8D20, 6YZ521, 6ZBP22, 7B2723, 7BWJ24, 7E2325, 7JMO26, 7K8M27, 6W4128, 6YM029, 6ZER30, 7DEO31, 7DPM32, 7MZF33, 7MZG34, 7MZJ35 |
| 2         | SARS-COV-2 variant              | Omicron: 7QNW36;  Gamma: 7NXB37;  Alpha: 7NEH38              |
| 3         | SARS-COV-2 NTD-binding antibody | 7LY3, 7M8J, 7L2C                                             |



## 2. MD Simulations

​	All the MD simulations were carried out by using the GROMACS 2021.3 package40,41 with Amber ff14sb force field42. Detailed settings of the MD simulation are provided in the supplementary file. Briefly, missing residues and atoms in the downloaded structures were fixed by BioExcel (package?), followed by energy minimization (EM; steepest descent) steps. Water molecules embedded in the original experimental structure were preserved in the simulation. The EM step minimization is deemed to be converged when the maximum force is smaller than 900 kJ/mol/nm. TIP3P water was used in the system to fill the simulation box. The minimum distance between the protein molecule and the box was set as 15Å. NaCl at the effective concentration of 150 mM was added to the simulation box. Extra Na+ or Cl- ions were added to neutralize the simulation system. Then three EM steps were carried out to optimize the system. Heavy atoms within protein were position-restrained, with the force constant set as 1000 kJ/mol/nm2 and the maximum force set as 500 kJ/mol/nm in the 1st EM. Position restraints were turned off in the next two EM steps. The maximum force was reduced to 300 and then 200 kJ/mol/nm in the second and third EM steps.

​	After EM steps, two position-restraint MD steps under the condition of NVT and NPT were used to equilibrate the system.  simulation, the temperature was gradually increased from 0 to 298 by simulated annealing, and NPT continued while both the protein was constrained with a spring constant of 1000kj/mol/nm^2. Finally, to sample more adequate conformational variations, 10ns relaxed NPT simulations were performed. 



## 3.Identification of interfacial water molecules between RBD and antibody
​	Interfacial water molecules have distinct structural and dynamic properties from those of bulk water and behave more like the water molecules in the hydration layer on the protein surface. Two metrics were used in the identification of surface water:

1. RMSF. The RMSF of interfacial waters should be comparable to that of protein atoms within a short period of simulation time. The exchange of water molecules at the same position should also be considered.

2. The shortest distance to protein all atoms. The longest length of H-bond is ~3.5 angstrom39. The length of H-O bond inside water molecular is 0.96 angstrom [xxx]. So, interfacial water molecules should be within 3.5 angstrom of the protein atoms. Furthermore, if the shortest distance from a water molecule to the RBD and the shortest distance to the antibody are both less than 3.5, then it is on the binding surface. When considering the difference between these two distances, the water molecule can Recognized as part of RBD or antibody, in subsequent mmpbsa calculations.

	The sliding window algorithm is introduced to use the above two indicators at the same time. In detail, first, a window with a length of 50ps on the time axis is added to the simulation trajectory. For the trajectory segment covered by the window, water molecules are screened according to the RMSF less than 0.3 and the sum of the shortest distance to the RBD and antibody is less than 7 (3.5*2=7). Then, the window moves backward by 50ps, and the above screening process is repeated on the new trajectory segment. Finally, a interfacial water molecule index table that changes dynamically over time is obtained, which will guide the grouping of ligands and receptors for all interfacial water molecules in mmpbsa calculations.
	Since the hydrated water molecule of binding interface maintains a stable RMSF close to the protein atom before escaping, our method of calculating its RMSF and distance in the short term is efficient. In theory, the shorter the length of the sliding window, the more reliable the screening of interfacial water. The visually friendly process is shown as Fig.1.

	![](C:\Users\wurui\Pictures\a论文\screening process.jpg)
	
	In the MMPBSA calculation, each sampling frame can be located by a window, corresponding to a grouping index table. We report a python script to search hydration-water and then perform MMPBSA calculation for automated mass production. It has user-friendly interface and needs only topology and trajectory file from Gromacs40,41 and some few parameters, which can be fined in supplementary.



## 4. Screening MMPBSA
After the MD simulation and selection of interfacial water molecules, the standard procedure of MM/PBSA was applied to extract the binding free energy between RBD and antibody (ref?).



# Result and Discussion
## 1. Compare with non-waters
​	We started from 21 structures of the WT RBD in complex with different antibodies in our binding free energy calculation. We evaluated the impact on the MM/PBSA results by including interfacial water and compared it to the standard procedure of MM/PBSA (without any interfacial water). The selection of interfacial water molecules and the calculation of MM/PBSA were performed every 200ps and thus a total of 50 frames were used in the final analysis. The cross-plot of experiment binding affinity and computational free energy was fitted by Polynomial Regression algorithm [xxx]. Indeed, including interfacial water resulted in a significant increase in the Person’s constant, from 0.654 (standard) to 0.753 (with interfacial water) (Fig. 2). Therefore, for the RBD-Ab system, including explicitly treated interfacial water molecules produced better results compared to the standard MM/PBSA. 


## 2.Contributions to binding free energy by interfacial water 
​	To find out what role interfacial water plays in the binding of RBD to antibodies, we analyzed the interfacial water both energetically and structurally. According the MM energy of waters, Fig. 3(a) shows the change in number of interfacial waters in our calculation of RBD of SARS-CoV-2 Spike glycoprotein in complex with EY6A Fab(PDB ID: 6ZER6ZER). The affinity of 6ZER is ~2nM, and its experimental binding free energy is -12.018 kcal/mol. The result with waters is -16.929 kcal/mol while for non-waters is -9.253 kcal/mol, shown as Fig. 3(c-d). After 4ns, the amount of water decreases dramatically to 0 (Fig. 3a, This phenomenon rarely occurs in our full-data calculations), indicating that the decreasing of interfacial waters meets our search criteria. Those frames without interfacial water will not be calculated in our method. As shown in Fig. 3(d), when excluding all interfacial waters, the binding free energy is increasing and far away to the experimental value.



## 3.H-bond network
