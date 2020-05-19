
default_molecule_types = '''
Galpha(GPCR,GnP~GTP~GDP~none,p_site,mut~wt~mutant)
GPCR(Galpha,agonist)
Gbg(p_site)
agonist(p_site)
AChE(agonist)
RGS(Galpha)
Ga_effector(Galpha)
Gbg_effector(Gbg)
'''

default_species = '''
1 @c0:Galpha(GPCR,GnP~GDP,p_site,mut~wt) 20.0
2 @c0:Galpha(GPCR,GnP~GTP,p_site,mut~wt) 5.0
3 @c0:GPCR(Galpha,agonist) 5.0
4 @c0:Gbg(p_site) 50.0
5 @c0:agonist(p_site) 0.0
6 @c0:AChE(agonist) 0.0
7 @c0:RGS(Galpha) 30.0
8 @c0:Ga_effector(Galpha) 10.0
9 @c0:Gbg_effector(Gbg) 10.0
10 @c0:Galpha(GPCR,GnP~GDP,p_site,mut~mutant) 20.0
11 @c0:Galpha(GPCR,GnP~GTP,p_site,mut~mutant) 5.0
'''

default_observables = '''
Molecules O0_Galpha_tot @c0:Galpha()
Molecules O0_GPCR_tot @c0:GPCR()
Molecules GPCR_Galpha_GDP_Gbg @c0:GPCR(Galpha!1).Galpha(GPCR!1,GnP~GDP,p_site!2,mut).Gbg(p_site!2)
Molecules O0_Gbg_tot @c0:Gbg()
Molecules Ga_GTP @c0:Galpha(GPCR,GnP~GTP,p_site,mut)
Molecules G_trimer @c0:Galpha(GPCR,GnP~GDP,p_site!1,mut).Gbg(p_site!1)
Molecules Ga_GDP @c0:Galpha(GPCR,GnP~GDP,p_site,mut)
Molecules O0_agonist_tot @c0:agonist()
Molecules O0_AChE_tot @c0:AChE()
Molecules O0_RGS_tot @c0:RGS()
Molecules O0_Ga_effector_tot @c0:Ga_effector()
Molecules O0_Gbg_effector_tot @c0:Gbg_effector()
Molecules Ga_to_effector @c0:Galpha(GPCR,GnP,p_site!1,mut).Ga_effector(Galpha!1)
Molecules Gbg_to_effector @c0:Gbg(p_site!1).Gbg_effector(Gbg!1)
'''

default_reaction_rules = '''
agonist_to_GPCR_free:	    @c0:GPCR(Galpha,agonist) + @c0:agonist(p_site) <-> @c0:GPCR(Galpha,agonist!1).agonist(p_site!1)		1.0, 0.2
agonist_to_GPCR_w_Gtrimer:	@c0:GPCR(Galpha!+,agonist) + @c0:agonist(p_site) <-> @c0:GPCR(Galpha!+,agonist!1).agonist(p_site!1)		1.0, 0.062
AChE_to_agonist:	        @c0:AChE(agonist) + @c0:agonist(p_site) -> @c0:AChE(agonist!1).agonist(p_site!1)		50.0
Gbg_to_effector:	@c0:Gbg(p_site) + @c0:Gbg_effector(Gbg) <-> @c0:Gbg(p_site!1).Gbg_effector(Gbg!1)		4.0, 1.0
'''


equilibration='''
generate_network({max_iter=>10,max_agg=>100,overwrite=>1})
# Equilibration
# Note that the n_steps parameter controls only the reporting interval and not the step size used by the CVODE solver, which uses adaptive time stepping
simulate_ode({t_end=>100000,n_steps=>1000,sparse=>1,steady_state=>1}) 
'''

equilibration_long='''
generate_network({max_iter=>10,max_agg=>100,overwrite=>1})
simulate_ode({t_end=>1000000,n_steps=>10000,sparse=>1,steady_state=>1}) 
'''


agonist_ping = '''
simulate_ode({t_end=>20,n_steps=>20,atol=>1e-8,rtol=>1e-8,sparse=>1})
# now trigger the GPCRs by adding the agonist
setConcentration("@c0:agonist(p_site)", 60.0)
simulate_ode({continue=>1,t_end=>20.1,n_steps=>100,atol=>1e-8,rtol=>1e-8,sparse=>1})
# ... and then, remove the agonist
setConcentration("@c0:AChE(agonist)", 120.0)
simulate_ode({continue=>1,t_end=>220,n_steps=>500,atol=>1e-8,rtol=>1e-8,sparse=>1})
'''
