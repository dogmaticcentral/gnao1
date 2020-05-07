
# set the adjustable reaction rates here
def tweakable(tweaked={}):

	tweakable_retstr = "" # return string

	#  GTP->GDP conversion in mutant Galpha wt values are 0.07, 0.001
	if "Ga_catalysis" in tweaked:
		[kf, kr] = tweaked["Ga_catalysis"]
	else:
		[kf, kr] = [0.07, 0.001]
	tweakable_retstr += f"a2_Ga_catalysis:	c0:Galpha(GPCR,GnP~GTP,p_site,mut~mutant) <-> c0:Galpha(GPCR,GnP~GDP,p_site,mut~mutant)	{kf}, {kr} \n"

	# Gbg forming a complex with GDP-bound Galpha wt values are 6.0, 0
	if "G_trimer_formation" in tweaked:
		[kf, kr] = tweaked["G_trimer_formation"]
	else:
		[kf, kr] = [6.0, 0.6]
	tweakable_retstr += "b2_G_trimer_formation:	c0:Galpha(GPCR,GnP~GDP,p_site,mut~mutant) + "
	tweakable_retstr +=	f"c0:Gbg(p_site) -> c0:Galpha(GPCR,GnP~GDP,p_site!1,mut~mutant).Gbg(p_site!1) {kf} \n"

	# Gtrimer binding to activated GPCR wt values 10.0, 0.1
	if "Gtrimer_to_GPCR_active" in tweaked:
		[kf, kr] = tweaked["Gtrimer_to_GPCR_active"]
	else:
		[kf, kr] = [10.0, 0.1]
	tweakable_retstr += "d2_Gtrimer_to_GPCR_active:  c0:agonist(p_site!1).GPCR(Galpha,agonist!1) + c0:Galpha(GPCR,GnP~GDP,p_site!1,mut~mutant).Gbg(p_site!1) <-> "
	tweakable_retstr += f"c0:agonist(p_site!2).GPCR(Galpha!3,agonist!2).Galpha(GPCR!3,GnP~GDP,p_site!1,mut~mutant).Gbg(p_site!1)	{kf}, {kr}\n"

	# Gtrimer binding to  GPCR without agonist wt values 0.3, 0.1\
	if "Gtrimer_to_GPCR_free" in tweaked:
		[kf, kr] = tweaked["Gtrimer_to_GPCR_free"]
	else:
		[kf, kr] = [0.3, 0.1]
	tweakable_retstr += "e2_Gtrimer_to_GPCR_free:	c0:GPCR(Galpha,agonist) + c0:Galpha(GPCR,GnP~GDP,p_site!1,mut~mutant).Gbg(p_site!1) "
	tweakable_retstr += f"<-> c0:GPCR(Galpha!1,agonist).Galpha(GPCR!1,GnP~GDP,p_site!2,mut~mutant).Gbg(p_site!2)   {kf}, {kr}\n"

	# exchange GDP -> GTP in GPCR wt values 2.0. 0.0 
	if "GPCR_as_GEF" in tweaked:
		[kf, kr] = tweaked["GPCR_as_GEF"]
	else:
		[kf, kr] = [2.0, 0.0]
	tweakable_retstr += "g2_GPCR_as_GEF:	c0:GPCR(Galpha!1,agonist!+).Galpha(GPCR!1,GnP~GDP,p_site!2,mut~mutant).Gbg(p_site!2) "
	tweakable_retstr += f"-> c0:GPCR(Galpha,agonist!+) + c0:Galpha(GPCR,GnP~GTP,p_site,mut~mutant) + c0:Gbg(p_site) {kf} \n"

	# Galpha with GTP binding to RGS wt 2.0, 0.2
	if "RGS_to_Galpha_T" in tweaked:
		[kf, kr] = tweaked["RGS_to_Galpha_T" ]
	else:
		[kf, kr] = [2.0, 0.2]
	tweakable_retstr += "i2_RGS_to_Galpha_T:	c0:RGS(Galpha) + c0:Galpha(GPCR,GnP~GTP,p_site,mut~mutant) <-> "
	tweakable_retstr += f"c0:RGS(Galpha!1).Galpha(GPCR,GnP~GTP,p_site!1,mut~mutant)   {kf}, {kr}\n"

	# GTP to GDP catalyzed with the help of RGS wt 30.0, 0.0
	if "RGS_as_GAP" in tweaked:
		[kf, kr] = tweaked["RGS_as_GAP"]
	else:
		[kf, kr] = [30.0, 0.0]
	tweakable_retstr += "j2_RGS_as_GAP:	c0:RGS(Galpha!1).Galpha(GPCR,GnP~GTP,p_site!1,mut~mutant) -> "
	tweakable_retstr += f"c0:RGS(Galpha!1).Galpha(GPCR,GnP~GDP,p_site!1,mut~mutant) {kf}  \n"

	# releasing of GDP bound Galpha from RGS wt 100.0, 0.1
	if "RGS_to_Galpha_D" in tweaked:
		[kf, kr] = ["RGS_to_Galpha_D"]
	else:
		[kf, kr] = [100.0, 0.1]
	tweakable_retstr += "k2_RGS_to_Galpha_D:	c0:RGS(Galpha!1).Galpha(GPCR,GnP~GDP,p_site!1,mut~mutant) <-> "
	tweakable_retstr += f"c0:RGS(Galpha) + c0:Galpha(GPCR,GnP~GDP,p_site,mut~mutant) {kf},{kr}\n"

	# G_alpha binding to its effector (presumably adenylate cyclase)  4.0, 0.1
	if "G_alpha_T_to_effector" in tweaked:
		[kf, kr] = tweaked["G_alpha_T_to_effector"]
	else:
		[kf, kr] = [4.0, 0.1]

	tweakable_retstr += "l2_G_alpha_T_to_effector:	c0:Galpha(GPCR,GnP~GTP,p_site,mut~mutant) + c0:Ga_effector(Galpha) "
	tweakable_retstr += f"<-> c0:Galpha(GPCR,GnP~GTP,p_site!1,mut~mutant).Ga_effector(Galpha!1) {kf}, {kr} \n"

	# G_alpha binding to its effector (presumably adenylate cyclase)  4.0, 0.1
	if "G_alpha_T_to_effector" in tweaked:
		[kf, kr] = tweaked["G_alpha_T_to_effector"]
	else:
		[kf, kr] = [0.0, 0.1]
	tweakable_retstr += "l3_G_alpha_T_to_effector:	c0:Galpha(GPCR,GnP~none,p_site,mut~mutant) + c0:Ga_effector(Galpha) "
	tweakable_retstr += f"<-> c0:Galpha(GPCR,GnP~none,p_site!1,mut~mutant).Ga_effector(Galpha!1)  {kf}, {kr} \n"

	# empty Ga binding to GPCR,  1.0, 0.0
	if "Galpha_to_GPCR_free" in tweaked:
		[kf, kr] = tweaked["Galpha_to_GPCR_free"]
	else:
		[kf, kr] = [2.0, 0.5]
	tweakable_retstr += "e3_Galpha_to_GPCR_free:	c0:GPCR(Galpha,agonist) + c0:Galpha(GPCR,GnP~none,p_site,mut~mutant) <-> "
	tweakable_retstr += f"c0:GPCR(Galpha!1,agonist).Galpha(GPCR!1,GnP~none,p_site,mut~mutant)	  {kf}, {kr}  \n"

	return tweakable_retstr