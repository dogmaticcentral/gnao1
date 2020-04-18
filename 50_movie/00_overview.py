#! /usr/bin/python3
story = '''
[* 02] (We see the GPCR in the clump representation in the mebrane with G-tetramer docked.)
Upon activation GNAO1
[* 04] exchanges  GDP ->  GTP
[* 06] G-hetertrimer dissociates from GPCR, and (then) dissociates into Gbg and Galpha
[* 08] G-alpha undergoes conformational change 
[* 10] and docks to ADCY

The reverse process starts with
[* 12] RGS docking
[* 13] ATP breaks down into ADP and phosphate group
[* 14] the complex at ADCY dissociates, and Gnao re-associates with Gbg
[* 15] when the  GPCR is ready, the G-trimer will dock, open, and the cycle is ready to repeat itself

Switching to GNAO-centric view
[* 16] GNAO1 has two large  lobes or  domains [get transparent while morphing to uncover cartoon structure]
[* 17] the NBD has the following functions: Gbg binding, GPCR binding,  nucleotide binding/ catalytic site,
ADCY binding, RGS binding
[* 18] The helix domain does not seem to have any known interactions


[* 19] taking a more detatiled look into the structure of the proteins

The fucntion of the helix domain is less fien tuned
[* 20] as witnessed by the fact that it is as variable as the GTP-binding domain (GBD), 
[* 22] and yet has no known disease mutations - while they are all in  the GDB


In GDB
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
of the GTP->GDP catalysis, without completely abolishing the process.

[* 38] Even the  "empty pocket" mutants ... (?)

[* 40] MD casuiin mutations appear to be particularly intractable, given that no  therapy repored in GNAO1 cases so far compensates for
the non-existent interaction with the AC (which should be inhibitory). 

[* 42] While some borad correlation between the variant location on the protein an the reposnse to therapy  proves that
the therapy chioice can be rationalized, the wide diversity of effects that vaariants can induce in GNAO1 point to the
need of undestanding the particular variant in each patinet individually and tayloring the therapy accoringly.)

'''

print(story)
