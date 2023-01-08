## Including explicitly treated interfacial water molecules improved the free energy calculation for the binding of SARS-CoV-2 spike RBD and antibodies.



### Highlights：

1) Including explicitly treated interfacial water molecules improved the free energy calculation for the binding of SARS-CoV-2 spike RBD and antibodies.

2) Developed a screening method for interfacial water based on water molecular dynamics features and applied it to MM/PBSA.

3) Structural and energetic analyses confirmed the critical roles played by interfacial water in protein-protein binding.

4) Compromised interaction with interfacial water molecules partly explains the immune escape of the Omicron variant.

5) A standard python API to execute the proposed method with friendly user  interface and documents.


### Typical usage:

`
python /local_path/HyHOH/gmx/pygmx/main.py -tpr ../md_0.tpr -xtc ../md_0.xtc -ri 194 622 -li 1 193 -fm normal -rm hyhoh -t 5000 10000 20
`

-ri The residues index of receptor

-li The residues index of ligand

-fm Sampling method for MD frames, can be 'normal', 'most'(our method). 

-rm Methods to run MM/PBSA, can be 'normal', 'hyhoh'(Our method), 'dsthoh'.

-t Four numeric with spaces separate to specify the sampling interval. [begin, end, interval, frames_per_ps] in ps. This option can be ignored if a frames indexing list was assigned in source code by manual'


