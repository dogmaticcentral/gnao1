#! /usr/bin/python3
story = '''
[* 02] (We see the GPCR in the clump representation in the mebrane with G-tetramer docked.)
Upon activation GNAO1
[* 04] exchanges  ADP ->  ATP
[* 06] dissociates from GPCR
[* 08] dissociates from Gbg
[* 10] undergoes conformational change - the helix domain repositions to the place previously taken by Gbg

Otherwise the helix domain does not have any fine tuned function, 
[* 12] as witnessed by the fact that it is as variable as the ATP-binding domain (ABD), 
[* 14] and yet has no known disease mutations

In the ABD, the side 
[* 16] opposite to the helix domain and behind the catalytic pocket
[* 18] is the business end of the molecule, judged by tha available crystallograpic work
[* 20] in particular, the salt bridge E246-R209 seems to be directly at the interface with  adenylyl cyclase
(though some cautin: this is a model)
[* 22] RGS, the enhancer that speeds up the ATP->ADP catalysis by orders of magnitude also binds here.

It is at this location that we see the majority of 
[* 24] mutations related to MD only phenotype.
[* 26] deeper within this domain of GNAO1 molecule, mostly surrounding the catalytic pocket we encounter
the variants related to  E+MD phenotype
[* 28] somewhat mysteriolusly, away from the catalytic site, the binding interface with effectors, the helix domain
or Gbg in the tetramer conformation, we see the variants causing E -only phenotype
[* 30] the largest cluster of these residues are hydrophobic residues, conserved at least in their hydrophic nature
across the close GNAO1 paralogues in the human genome
[* 32] they are likely to be important for the proper folding of the protein. (see short traj from simulation ?)
It is tempting to hypothesize that such mutations would result in 
the complete misfolding and the loss of the expressed protein.  However no  frameshift or early stop codon alleles
have been reported in the GNAO1-related disorders, and such mutations would also result in the loss of the expressed protein.
4 such cases, though, have been reported in gnomAD. Thus we expect that the complete loss of one allele is an event that
the system can tolerate by regulating the dosage from the othere allele.

[* 34] For some reason the therapy used in these cases has been seldomly reported, except for  one isolated case where
a series of conventional therapties produced no effect (pubmed 29961512)

[* 36] The E+MD cases have one common characteristic that when tractable they respond to therapies targeted
at enhancing the weak baseline signal, such as increasing dopamine availability and enhancing the downstream GABA-A 
receptors.  This can be rationalized by considering the protein location of these varaints which may slow down the rate
of the ATP->ADP catalysis, without completely abolishing the process.

[* 38] Even the  "empty pocket" mutants ... (?)

[* 40] MD casuiin mutations appear to be particularly intractable, given that no  therapy repored in GNAO1 cases so far compensates for
the non-existent interaction with the AC (which should be inhibitory). 

[* 42] While some borad correlation between the variant location on the protein an the reposnse to therapy  proves that
the therapy chioice can be rationalized, the wide diversity of effects that vaariants can induce in GNAO1 point to the
need of undestanding the particular variant in each patinet individually and tayloring the therapy accoringly.)

'''

print(story)
