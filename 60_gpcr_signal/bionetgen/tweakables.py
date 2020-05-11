
model_template = '''
begin model

begin compartments
c0	3	1
end compartments

begin parameters
end parameters

begin molecule types
{molecule_types}
end molecule types

begin seed species
{species}
end seed species

begin observables
{observables}
end observables

begin functions
end functions

begin reaction rules
{reaction_rules}
end reaction rules


end model
'''

class Reaction:
	def __init__(self,  partner, subtype, nucleotide,  reaction_from, reaction_to, kf, kr):
		self.subtype = subtype # this can be wt or mutant - we might appropriate mutant to play s-type Galpha
		self.nucleotide = nucleotide # GnP means: any of GDP, GTP; it can also be "none" for empty pocket
		self.partner = partner
		self.reaction_from = reaction_from
		self.reaction_to = reaction_to
		self.direction = "<->"
		self.k_forward = kf
		self.k_reverse = kr

	def name(self):
		return "G_a_{}_{}_{}".format(self.subtype, self.nucleotide, self.partner)

	def prettyprint(self):
		rfrom = self.reaction_from.format(subtype=self.subtype)
		rto =  self.reaction_to.format(subtype=self.subtype)
		if self.direction=="->":
			# f"{self.name()}: {rfrom} {self.direction} {rto} {self.k_forward}\n"
			return "{}: {} {} {} {}".format(self.name(), rfrom, self.direction, rto, self.k_forward)
		elif self.direction=="<-":
			return "{}: {} {} {} {}".format(self.name(), rfrom, self.direction, rto, self.k_reverse)
		else:
			return "{}: {} {} {} {},{}".format(self.name(), rfrom, self.direction, rto, self.k_forward, self.k_reverse)

'''
1   1_Ga_catalysis:	\@c0:Galpha(GPCR,GnP~GTP,p_site,mut~wt) <-> \@c0:Galpha(GPCR,GnP~GDP,p_site,mut~wt)		0.07, 0.001
2   b1_G_trimer_formation:	\@c0:Galpha(GPCR,GnP~GDP,p_site,mut~wt) + \@c0:Gbg(p_site) -> \@c0:Galpha(GPCR,GnP~GDP,p_site!1,mut~wt).Gbg(p_site!1)		6.0
3   d1_Gtrimer_to_GPCR_active:	\@c0:agonist(p_site!1).GPCR(Galpha,agonist!1) + \@c0:Galpha(GPCR,GnP~GDP,p_site!1,mut~wt).Gbg(p_site!1) <-> \@c0:agonist(p_site!2).GPCR(Galpha!3,agonist!2).Galpha(GPCR!3,GnP~GDP,p_site!1,mut~wt).Gbg(p_site!1)		10.0, 0.1
4   e1_Gtrimer_to_GPCR_free:	\@c0:GPCR(Galpha,agonist) + \@c0:Galpha(GPCR,GnP~GDP,p_site!1,mut~wt).Gbg(p_site!1) <-> \@c0:GPCR(Galpha!1,agonist).Galpha(GPCR!1,GnP~GDP,p_site!2,mut~wt).Gbg(p_site!2)		0.3, 0.1
5   g1_GPCR_as_GEF:	\@c0:GPCR(Galpha!1,agonist!+).Galpha(GPCR!1,GnP~GDP,p_site!2,mut~wt).Gbg(p_site!2) -> \@c0:GPCR(Galpha,agonist!+) + \@c0:Galpha(GPCR,GnP~GTP,p_site,mut~wt) + \@c0:Gbg(p_site)		2.0
6   i1_RGS_to_Galpha_T:	\@c0:RGS(Galpha) + \@c0:Galpha(GPCR,GnP~GTP,p_site,mut~wt) <-> \@c0:RGS(Galpha!1).Galpha(GPCR,GnP~GTP,p_site!1,mut~wt)		2.0, 0.2
7   j1_RGS_as_GAP:	\@c0:RGS(Galpha!1).Galpha(GPCR,GnP~GTP,p_site!1,mut~wt) -> \@c0:RGS(Galpha!1).Galpha(GPCR,GnP~GDP,p_site!1,mut~wt)		30.0
8   k1_RGS_to_Galpha_D:	\@c0:RGS(Galpha!1).Galpha(GPCR,GnP~GDP,p_site!1,mut~wt) <-> \@c0:RGS(Galpha) + \@c0:Galpha(GPCR,GnP~GDP,p_site,mut~wt)		100.0, 0.1
9   l1_G_alpha_T_to_effector:	\@c0:Galpha(GPCR,GnP~GTP,p_site,mut~wt) + \@c0:Ga_effector(Galpha) <-> \@c0:Galpha(GPCR,GnP~GTP,p_site!1,mut~wt).Ga_effector(Galpha!1)		4.0, 0.1
f_agonist_to_GPCR_w_Gtrimer:	\@c0:GPCR(Galpha!+,agonist) + \@c0:agonist(p_site) <-> \@c0:GPCR(Galpha!+,agonist!1).agonist(p_site!1)		1.0, 0.062

'''
def set_default_galpha_reaction_rules(subtype):

	reactions = []

	# 1 GTP to GDP catalysis without RGS
	r = Reaction("GTP", subtype, "GTP",
	             "c0:Galpha(GPCR,GnP~GTP,p_site,mut~{subtype})",
	             "c0:Galpha(GPCR,GnP~GDP,p_site,mut~{subtype})",
	             0.07, 0.001)
	reactions.append(r)

	# 2 G-trimer formation, with GDP bound to Ga
	r = Reaction("Gbg",  subtype, "GDP",
	             "c0:Galpha(GPCR,GnP~GDP,p_site,mut~{subtype}) + c0:Gbg(p_site)",
	             "c0:Galpha(GPCR,GnP~GDP,p_site!1,mut~{subtype}).Gbg(p_site!1)",
	             6.0, 0.0)
	r.direction = "->"
	reactions.append(r)

	# 3 G-trimer binding to activated GPCR
	# An incomplete bond may also be specified using “!+”,
	# where the  wild  card  “+”  indicates  that  the  identity  of  the  binding partner
	# is irrelevant  for  purposes  of  matching.
	r = Reaction("GPCR_activated",  subtype, "GDP",
	             "c0:GPCR(Galpha,agonist!+) + c0:Galpha(GPCR,GnP~GDP,p_site!1,mut~{subtype}).Gbg(p_site!1)",
	             "c0:GPCR(Galpha!3,agonist!+).Galpha(GPCR!3,GnP~GDP,p_site!1,mut~{subtype}).Gbg(p_site!1)",
	             10.0, 0.1)
	reactions.append(r)

	# 4 G-trimer binding to  GPCR without agonist
	r = Reaction("GPCR_free", subtype,  "GDP",
	             "c0:GPCR(Galpha,agonist) + c0:Galpha(GPCR,GnP~GDP,p_site!1,mut~{subtype}).Gbg(p_site!1)",
	             "c0:GPCR(Galpha!1,agonist).Galpha(GPCR!1,GnP~GDP,p_site!2,mut~{subtype}).Gbg(p_site!2)",
	             0.3, 0.1)
	reactions.append(r)

	# 5  exchange GDP -> GTP when bound to GPCR
	r = Reaction("GPCR_as_GEF",  subtype, "GDP",
	             "c0:GPCR(Galpha!1,agonist!+).Galpha(GPCR!1,GnP~GDP,p_site!2,mut~{subtype}).Gbg(p_site!2)",
	             "c0:GPCR(Galpha,agonist!+) + c0:Galpha(GPCR,GnP~GTP,p_site,mut~{subtype}) + c0:Gbg(p_site)",
	             2.0, 0.0)
	r.direction = "->"
	reactions.append(r)

	# 6  Galpha with GTP binding to RGS
	r = Reaction("RGS",  subtype, "GTP",
	             "c0:RGS(Galpha) + c0:Galpha(GPCR,GnP~GTP,p_site,mut~{subtype})",
	             "c0:RGS(Galpha!1).Galpha(GPCR,GnP~GTP,p_site!1,mut~{subtype})",
	             2.0, 0.2)

	reactions.append(r)

	# 7  GTP to GDP catalyzed with the help of RGS
	r = Reaction("RGS_as_GAP",  subtype, "GTP",
	             "c0:RGS(Galpha!1).Galpha(GPCR,GnP~GTP,p_site!1,mut~{subtype})",
	             "c0:RGS(Galpha!1).Galpha(GPCR,GnP~GDP,p_site!1,mut~{subtype})",
	             30.0, 0)
	r.direction = "->"
	reactions.append(r)

	# 8 releasing of GDP from Galpha bount to RGS
	r = Reaction("GDP",  subtype, "GDP",
	             "c0:RGS(Galpha!1).Galpha(GPCR,GnP~GDP,p_site!1,mut~{subtype})",
	             "c0:RGS(Galpha) + c0:Galpha(GPCR,GnP~GDP,p_site,mut~{subtype})",
	             100.0, 0.1)
	reactions.append(r)

	# 9 Galpha binding to its effector (presumably adenylate cyclase)
	r = Reaction("effector",  subtype, "GTP",
	             "c0:Galpha(GPCR,GnP~GTP,p_site,mut~{subtype}) + c0:Ga_effector(Galpha)",
	             "c0:Galpha(GPCR,GnP~GTP,p_site!1,mut~{subtype}).Ga_effector(Galpha!1)",
	             4.0, 0.1)
	reactions.append(r)

	return reactions


###################
def empty_pocket_reaction_rules():
	reactions = []

	# empty pocket Ga binding to free GPCR
	# r = Reaction("GPCR_free",  subtype, "none",
	#              "c0:GPCR(Galpha,agonist) + c0:Galpha(GPCR,GnP~none,p_site,mut~{subtype})",
	#              "c0:GPCR(Galpha!1,agonist).Galpha(GPCR!1,GnP~none,p_site,mut~{subtype})",
	#              10.0, 0.1)
	# reactions.append(r)
	#
	# # 3 G-trimer binding to activated GPCR
	# r = Reaction("GPCR_activated",  subtype, "none",
	#              "c0:GPCR(Galpha,agonist!+) + c0:Galpha(GPCR,GnP~none,p_site,mut~{subtype})",
	#              "c0:GPCR(Galpha!1,agonist!+).Galpha(GPCR!1,GnP~none,p_site,mut~{subtype})",
	#              10.0, 0.1)
	# reactions.append(r)

	# language reference: https://www.csb.pitt.edu/Faculty/Faeder/Publications/Reprints/Faeder_2009.pdf
	#  A wild card, “?”, may be used to indicate that a match may occur regardless
	#  of  whether  a  bond  is  present  or  absent.
	#  For example, the two patterns EGFR(Y1068~P) and EGFR(Y1068~P!?) are not equivalent:
	#  The  EGFR(Y1068~P)  pattern selects only EGFR molecules in which the Y1068 component is
	#  phosphorylated  and  unbound, whereas  the EGFR(Y1068~P!?) pattern selects  all  EGFR  molecules
	#  in  which  the  Y1068  component is  phosphorylated, irrespective  of the existence of the bond

	# empty pocket Ga binding to  GPCR irrespective of the agonist
	# 3 G-trimer binding to activated GPCR - note the absence of Gbg
	r = Reaction("GPCR_activated",  "mutant", "none",
	             "c0:GPCR(Galpha,agonist!?) + c0:Galpha(GPCR,GnP~none,p_site,mut~{subtype})",
	             "c0:GPCR(Galpha!1,agonist!?).Galpha(GPCR!1,GnP~none,p_site,mut~{subtype})",
	             10.0, 0.1)
	reactions.append(r)
	# empty pocket does not bind Gbg (N Yu and M simon, 1998)
	# basically, id toes nothing except eliminating itself from the Galpha pool


	return reactions


###################
def set_tweaked_reaction_rules(subtype="wt", tweaks = None):

	reactions = set_default_galpha_reaction_rules(subtype)
	if not tweaks: return reactions

	for partner, [kf, kr] in tweaks.items():
		found = False
		for r in reactions:
			if r.partner != partner: continue
			r.k_forward = kf
			r.k_reverse = kr
			found = True
			break
		if not found:
			print(partner, "not found")
			exit()

	return reactions


###################
def reaction_rules_string(reactions):
	return "\n".join([r.prettyprint() for r in reactions])+"\n"




def galpha_empty_species():
	spec  = "12 @c0:Galpha(GPCR,GnP~none,p_site,mut~mutant) 25.0\n"
	return spec

###################
def galpha_s_species(double=False):
	if double: # when we have Gnao Wt and mutant
		spec  = "13 @c0:Galpha(GPCR,GnP~GDP,p_site,mut~s) 40.0\n"
		spec += "14 @c0:Galpha(GPCR,GnP~GTP,p_site,mut~s) 10.0\n"
	else:
		spec  = "13 @c0:Galpha(GPCR,GnP~GDP,p_site,mut~s) 20.0\n"
		spec += "14 @c0:Galpha(GPCR,GnP~GTP,p_site,mut~s) 5.0\n"

	return spec

