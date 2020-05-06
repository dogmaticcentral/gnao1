#!/usr/bin/perl

sub set_literals();
sub set_tweakable();

####################################################
my $literal1 = "";
my $literal2 = "";

my $tweakable = "";

set_literals() ;
set_tweakable() ;
# run first with the rootnm base
my $rootnm = "production";
my $svg = 0;

open (OUTF, ">$rootnm.bngl") || die "Cno  $rootnm.bngl: $!\n";
print OUTF $literal1;
print  OUTF $tweakable;
print  OUTF $literal2;
close OUTF;

# I also have alias bing='/home/ivana/third/BioNetGen-2.3.1/BNG2.pl  "$@"'
# in bashrc, but I do not need it here
# run BioNetGen, and print its output
print `/home/ivana/third/bionetgen/BNG2.pl $rootnm.bngl`;

# plot the output
open (OUTF, ">$rootnm.gplt") || die "Cno  $rootnm.gplt: $!\n";
print OUTF "# special character I need: the dot {\267}\n";
print OUTF "set encoding iso_8859_1\n";

if ($svg) {
    print OUTF "# svg\n";
    print OUTF "set terminal svg size 410,250 fname 'Verdana, Helvetica, Arial, sans-serif' fsize '9' rounded dashed \n";
    print OUTF "set output '$rootnm.svg' \n";
    print OUTF "set size ratio 1  \n";
} else {
    print OUTF "# png\n";
    print OUTF "set terminal pngcairo size 615,500 enhanced font 'Verdana,18'\n";
    print OUTF "set output '$rootnm.png' \n";
}

print OUTF "# define axis\n";
print OUTF "# remove border on top and right and set color to gray\n";
print OUTF "set style line 11 lc rgb '#808080' lt 1\n";
print OUTF "set border 3 back ls 11\n";
print OUTF "set tics nomirror\n";
print OUTF "set tics font \", 16\"\n";
print OUTF "# define grid\n";
print OUTF "# set style line 12 lc rgb '#808080' lt 0 lw 1\n";
print OUTF "# set grid back ls 12\n";

print OUTF "# color definitions\n";
print OUTF "set style line 1 lc rgb '#004c91' lt 1 lw 3 # --- blue\n";
print OUTF "set style line 2 lc rgb '#ffc220' lt 1 lw 3 # --- yellow\n";

print OUTF "#set key top right\n";

print OUTF "set xlabel 'seconds'\n";
print OUTF "set ylabel 'molecular population size (%)'\n";

print OUTF "set xrange [0:200]\n";
print OUTF "set yrange [0:100]\n";

print OUTF "labelBG =  'G_{/Symbol b}_{/Symbol g}{\\267}effector'\n";
print OUTF "labelA  =  'G_{/Symbol a}{\\267}effector'\n";


if ($rootnm  ne "base") {

    print OUTF "set style line 3 lc rgb '#004c91'  lt 2 lw 1 dt '.' # --- blue\n";
    print OUTF "set style line 4 lc rgb '#ffc220'  lt 2 lw 1 dt '.' # --- yellow\n";


    print OUTF "plot '$rootnm.gdat' u 1:(\$15/50*100)  t labelBG w lines ls 2, '$rootnm.gdat' u 1:(\$14/50*100)  t labelA w lines ls 1, ";
    print OUTF "'base.gdat' u 1:(\$14/50*100) notitle  w lines ls 3,   'base.gdat' u 1:(\$15/50*100) notitle  w lines ls 4\n";

} else {
  print OUTF "plot '$rootnm.gdat' u 1:(\$15/50*100)  t 'Gbg+effector' w lines ls 2,'$rootnm.gdat' u 1:(\$14/50*100)  t 'Ga+effector' w lines ls 1\n";

}

print `gnuplot $rootnm.gplt`;
$svg ||  print "\nnow run eog $rootnm.png\n";





exit(0);

###################################################
sub set_tweakable() {
	# set the adjustable reaction rates here
	my $kf, $kr;

	#  GTP->GDb converstion in mutant Galpha; wt values are 0.07, 0.001
	$kf = "0.07";
	$kr = "0.001";
	$tweakable .= "a2_Ga_catalysis:	\@c0:Galpha(GPCR,GnP~GTP,p_site,mut~mutant) <-> \@c0:Galpha(GPCR,GnP~GDP,p_site,mut~mutant)	$kf, $kr \n";

	# Gbg forming a complex with GDP-bound Galpha; wt values are 6.0, 0
	$kf = "6.0";
	$kr = "0.0";
	$tweakable .= "b2_G_trimer_formation:	\@c0:Galpha(GPCR,GnP~GDP,p_site,mut~mutant) + \@c0:Gbg(p_site) -> \@c0:Galpha(GPCR,GnP~GDP,p_site!1,mut~mutant).Gbg(p_site!1)	 $kf \n";

	# Gtrimer binding to activated GPCR; wt values 10.0, 0.1
	$kf = "10.0";
	$kr = "0.1";
	$tweakable .= "d2_Gtrimer_to_GPCR_active:	\@c0:agonist(p_site!1).GPCR(Galpha,agonist!1) + \@c0:Galpha(GPCR,GnP~GDP,p_site!1,mut~mutant).Gbg(p_site!1) <-> ";
	$tweakable .= "\@c0:agonist(p_site!2).GPCR(Galpha!3,agonist!2).Galpha(GPCR!3,GnP~GDP,p_site!1,mut~mutant).Gbg(p_site!1)	 $kf, $kr	 \n";

	# Gtrimer binding to  GPCR without agonist; wt values 0.3, 0.1\
	$kf = "0.3";
	$kr = "0.1";
	$tweakable .= "e2_Gtrimer_to_GPCR_free:	\@c0:GPCR(Galpha,agonist) + \@c0:Galpha(GPCR,GnP~GDP,p_site!1,mut~mutant).Gbg(p_site!1) <-> \@c0:GPCR(Galpha!1,agonist).Galpha(GPCR!1,GnP~GDP,p_site!2,mut~mutant).Gbg(p_site!2)    $kf, $kr	\n";

	# exchange GDP -> GTP in GPCR; wt values 2.0. 0.0
	$kf = "2.0";
	$kr = "0.0";
	$tweakable .= "g2_GPCR_as_GEF:	\@c0:GPCR(Galpha!1,agonist!+).Galpha(GPCR!1,GnP~GDP,p_site!2,mut~mutant).Gbg(p_site!2) -> \@c0:GPCR(Galpha,agonist!+) + \@c0:Galpha(GPCR,GnP~GTP,p_site,mut~mutant) + \@c0:Gbg(p_site)	$kf \n";

	# Galpha with GTP binding to RGS; wt 2.0, 0.2
	$kf = "2.0";
	$kr = "0.2";
	$tweakable .= "i2_RGS_to_Galpha_T:	\@c0:RGS(Galpha) + \@c0:Galpha(GPCR,GnP~GTP,p_site,mut~mutant) <-> \@c0:RGS(Galpha!1).Galpha(GPCR,GnP~GTP,p_site!1,mut~mutant)	     $kf, $kr \n";

	# GTP to GDP catalyzed with the help of RGS; wt 30.0, 0.0
	$kf = "30.0";
	$kr = "0.0";
	$tweakable .= "j2_RGS_as_GAP:	\@c0:RGS(Galpha!1).Galpha(GPCR,GnP~GTP,p_site!1,mut~mutant) -> \@c0:RGS(Galpha!1).Galpha(GPCR,GnP~GDP,p_site!1,mut~mutant)	$kf  \n";

	# releasing of GDP bound Galpha from RGS; wt 100.0, 0.1
	$kf = "100.0";
	$kr = "0.1";
	$tweakable .= "k2_RGS_to_Galpha_D:	\@c0:RGS(Galpha!1).Galpha(GPCR,GnP~GDP,p_site!1,mut~mutant) <-> \@c0:RGS(Galpha) + \@c0:Galpha(GPCR,GnP~GDP,p_site,mut~mutant)	$kf, $kr \n";

	# G_alpha binding to its effector (presumably adenylate cyclase);  4.0, 0.1
	$kf = "4.0";
	$kr = "0.1";
	$tweakable .= "l2_G_alpha_T_to_effector:	\@c0:Galpha(GPCR,GnP~GTP,p_site,mut~mutant) + \@c0:Ga_effector(Galpha) <-> \@c0:Galpha(GPCR,GnP~GTP,p_site!1,mut~mutant).Ga_effector(Galpha!1)	$kf, $kr \n";

	# G_alpha binding to its effector (presumably adenylate cyclase);  4.0, 0.1
	$kf = "0.0";
	$kr = "0.1";
	$tweakable .= "l3_G_alpha_T_to_effector:	\@c0:Galpha(GPCR,GnP~none,p_site,mut~mutant) + \@c0:Ga_effector(Galpha) <-> \@c0:Galpha(GPCR,GnP~none,p_site!1,mut~mutant).Ga_effector(Galpha!1)	$kf, $kr \n";

	# empty Ga binding to GPCR,  1.0, 0.0
	$kf = "2.0";
	$kr = "0.5";
	$tweakable .= "e3_Galpha_to_GPCR_free:	\@c0:GPCR(Galpha,agonist) + \@c0:Galpha(GPCR,GnP~none,p_site,mut~mutant) <-> \@c0:GPCR(Galpha!1,agonist).Galpha(GPCR!1,GnP~none,p_site,mut~mutant)	 $kf, $kr"

}


####################################################


sub set_literals() {

	$literal1 = "
begin model

begin compartments
c0	3	1
end compartments

begin parameters
end parameters

begin molecule types
Galpha(GPCR,GnP~GTP~GDP~none,p_site,mut~wt~mutant)
GPCR(Galpha,agonist)
Gbg(p_site)
agonist(p_site)
AChE(agonist)
RGS(Galpha)
Ga_effector(Galpha)
Gbg_effector(Gbg)
end molecule types

begin seed species
1 \@c0:Galpha(GPCR,GnP~GDP,p_site,mut~wt) 20.0
2 \@c0:Galpha(GPCR,GnP~GTP,p_site,mut~wt) 5.0
3 \@c0:GPCR(Galpha,agonist) 50.0
4 \@c0:Gbg(p_site) 50.0
5 \@c0:agonist(p_site) 0.0
6 \@c0:AChE(agonist) 0.0
7 \@c0:RGS(Galpha) 30.0
8 \@c0:Ga_effector(Galpha) 50.0
9 \@c0:Gbg_effector(Gbg) 50.0
10 \@c0:Galpha(GPCR,GnP~GDP,p_site,mut~mutant) 20.0
11 \@c0:Galpha(GPCR,GnP~GTP,p_site,mut~mutant) 5.0
12 \@c0:Galpha(GPCR,GnP~none,p_site,mut~mutant) 0.0
end seed species

begin observables
Molecules O0_Galpha_tot \@c0:Galpha()
Molecules O0_GPCR_tot \@c0:GPCR()
Molecules GPCR_Galpha_GDP_Gbg \@c0:GPCR(Galpha!1).Galpha(GPCR!1,GnP~GDP,p_site!2,mut).Gbg(p_site!2)
Molecules O0_Gbg_tot \@c0:Gbg()
Molecules Ga_GTP \@c0:Galpha(GPCR,GnP~GTP,p_site,mut)
Molecules G_trimer \@c0:Galpha(GPCR,GnP~GDP,p_site!1,mut).Gbg(p_site!1)
Molecules Ga_GDP \@c0:Galpha(GPCR,GnP~GDP,p_site,mut)
Molecules O0_agonist_tot \@c0:agonist()
Molecules O0_AChE_tot \@c0:AChE()
Molecules O0_RGS_tot \@c0:RGS()
Molecules O0_Ga_effector_tot \@c0:Ga_effector()
Molecules O0_Gbg_effector_tot \@c0:Gbg_effector()
Molecules Ga_to_effector \@c0:Galpha(GPCR,GnP,p_site!1,mut).Ga_effector(Galpha!1)
Molecules Gbg_to_effector \@c0:Gbg(p_site!1).Gbg_effector(Gbg!1)
end observables

begin functions
end functions

begin reaction rules
a1_Ga_catalysis:	\@c0:Galpha(GPCR,GnP~GTP,p_site,mut~wt) <-> \@c0:Galpha(GPCR,GnP~GDP,p_site,mut~wt)		0.07, 0.001
e1_Gtrimer_to_GPCR_free:	\@c0:GPCR(Galpha,agonist) + \@c0:Galpha(GPCR,GnP~GDP,p_site!1,mut~wt).Gbg(p_site!1) <-> \@c0:GPCR(Galpha!1,agonist).Galpha(GPCR!1,GnP~GDP,p_site!2,mut~wt).Gbg(p_site!2)		0.3, 0.1
g1_GPCR_as_GEF:	\@c0:GPCR(Galpha!1,agonist!+).Galpha(GPCR!1,GnP~GDP,p_site!2,mut~wt).Gbg(p_site!2) -> \@c0:GPCR(Galpha,agonist!+) + \@c0:Galpha(GPCR,GnP~GTP,p_site,mut~wt) + \@c0:Gbg(p_site)		2.0
b1_G_trimer_formation:	\@c0:Galpha(GPCR,GnP~GDP,p_site,mut~wt) + \@c0:Gbg(p_site) -> \@c0:Galpha(GPCR,GnP~GDP,p_site!1,mut~wt).Gbg(p_site!1)		6.0
c_agonist_to_GPCR_free:	\@c0:GPCR(Galpha,agonist) + \@c0:agonist(p_site) <-> \@c0:GPCR(Galpha,agonist!1).agonist(p_site!1)		1.0, 0.2
h_AChE_to_agonist:	\@c0:AChE(agonist) + \@c0:agonist(p_site) -> \@c0:AChE(agonist!1).agonist(p_site!1)		50.0
d1_Gtrimer_to_GPCR_active:	\@c0:agonist(p_site!1).GPCR(Galpha,agonist!1) + \@c0:Galpha(GPCR,GnP~GDP,p_site!1,mut~wt).Gbg(p_site!1) <-> \@c0:agonist(p_site!2).GPCR(Galpha!3,agonist!2).Galpha(GPCR!3,GnP~GDP,p_site!1,mut~wt).Gbg(p_site!1)		10.0, 0.1
f_agonist_to_GPCR_w_Gtrimer:	\@c0:GPCR(Galpha!+,agonist) + \@c0:agonist(p_site) <-> \@c0:GPCR(Galpha!+,agonist!1).agonist(p_site!1)		1.0, 0.062
i1_RGS_to_Galpha_T:	\@c0:RGS(Galpha) + \@c0:Galpha(GPCR,GnP~GTP,p_site,mut~wt) <-> \@c0:RGS(Galpha!1).Galpha(GPCR,GnP~GTP,p_site!1,mut~wt)		2.0, 0.2
j1_RGS_as_GAP:	\@c0:RGS(Galpha!1).Galpha(GPCR,GnP~GTP,p_site!1,mut~wt) -> \@c0:RGS(Galpha!1).Galpha(GPCR,GnP~GDP,p_site!1,mut~wt)		30.0
k1_RGS_to_Galpha_D:	\@c0:RGS(Galpha!1).Galpha(GPCR,GnP~GDP,p_site!1,mut~wt) <-> \@c0:RGS(Galpha) + \@c0:Galpha(GPCR,GnP~GDP,p_site,mut~wt)		100.0, 0.1
l1_G_alpha_T_to_effector:	\@c0:Galpha(GPCR,GnP~GTP,p_site,mut~wt) + \@c0:Ga_effector(Galpha) <-> \@c0:Galpha(GPCR,GnP~GTP,p_site!1,mut~wt).Ga_effector(Galpha!1)		4.0, 0.1
m_Gbg_to_effector:	\@c0:Gbg(p_site) + \@c0:Gbg_effector(Gbg) <-> \@c0:Gbg(p_site!1).Gbg_effector(Gbg!1)		4.0, 1.0
";

	$literal2 = "
end reaction rules

end model

generate_network({max_iter=>10,max_agg=>100,overwrite=>1})
# Equilibration
# Note that the n_steps parameter controls only the reporting interval and not the step size used by the CVODE solver, which uses adaptive time stepping
simulate_ode({t_end=>100000,n_steps=>100,sparse=>1,steady_state=>1})

# now trigger the GPCRs by adding the agonist
simulate_ode({t_end=>20,n_steps=>20,atol=>1e-8,rtol=>1e-8,sparse=>1})
setConcentration(\"\@c0:agonist(p_site)\", 60.0)
simulate_ode({continue=>1,t_end=>20.1,n_steps=>100,atol=>1e-8,rtol=>1e-8,sparse=>1})
setConcentration(\"\@c0:AChE(agonist)\", 120.0)
simulate_ode({continue=>1,t_end=>220,n_steps=>500,atol=>1e-8,rtol=>1e-8,sparse=>1})
";



}
