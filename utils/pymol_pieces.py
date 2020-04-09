
from math import acos, sqrt, fabs, sin

from pymol import cmd
from math import pi, cos

# from pyquaternion import Quaternion # will not take matrix that is (numerically) non-orthogona;
import numpy as np
import quaternion  # pip3 install numpy-quaternion and numba


def pymol_chdir(outdir):
	cmd.cd(outdir)


def ghost_rep(selection):
	cmd.hide("everything", selection)
	cmd.color("white", selection)
	cmd.set("transparency", 0.5, selection)
	cmd.show("surface", selection)


def interface_outline(reference_selection, interactant, color):
	cmd.hide("everything", interactant)
	name = "{}_if".format(interactant.replace("(","").replace(")","").replace(" ",""))
	# move selection to its own object
	cmd.create(name, "byres {} within 3 of {}".format(interactant, reference_selection))
	cmd.color(color, name)
	cmd.set("transparency", 0.7, name)
	cmd.show("surface", name)


def load_structures(structure_home, structure_filename, structures):


	for structure in structures:
		cmd.load("{}/{}".format(structure_home, structure_filename[structure]), structure)
		cmd.hide("everything", structure)

	if "AC" in structures:
		cmd.remove("AC and chain B")

	if "GPCR" in structures:
		cmd.remove("resi 343-1200 and GPCR")


	return


def residue_color(main_object, res_id, rgb_color_list):
	color_name  = "cnm_{}_{}".format(main_object, res_id)
	cmd.set_color(color_name, rgb_color_list)
	cmd.color(color_name, "{} and resi {}".format(main_object, res_id))


def clump_representation(regions, color, name, transparency=0.0):
	# show regions as clumps or blobs
	cmd.set("surface_quality", 1)
	if transparency>0: cmd.set("transparency", transparency, name)
	cmd.alter("all", "b=50")
	cmd.alter("all", "q=1")
	cmd.set("gaussian_resolution", 7)

	idx = 0
	for region in regions:
		idx += 1
		mapname  = "map_{}_{}".format(name, idx)
		surfname = "surf_{}_{}".format(name, idx)
		cmd.map_new(mapname, "gaussian", 1, region, 2)
		cmd.isosurface(surfname, mapname)
		cmd.color(color, surfname)


def view_string2view(view_string):
	return [float(f) for f in view_string.split(",")]


def view2view_string(view):
	return ",".join([str(f) for f in view])


def view2quat(view):

	array  = np.array(view[:9])
	matrix = np.asmatrix(array.reshape((3,3)))
	return quaternion.from_rotation_matrix(matrix)



def view_interpolate(view_init_str, view_last_str, number_of_frames=5, frameno_offset=0):

	view_init = view_string2view(view_init_str)
	view_last = view_string2view(view_last_str)
	# get quaternions for interpolation
	qstart = view2quat(view_init)
	qend   = view2quat(view_last)

	last_frame = frameno_offset
	for frameno in range(1, number_of_frames+1):
		# (quat start, quat end, time start, time end, time evaluated)
		# of example (qstart, qend, 0, number_of_frames, frameno)
		# (qstart, qend, 0, 1, frameno/number_of_frames) is npt working, even though these
		# numbers ar esupposed to be floats - some numpy shit I would guess
		#  easing: (qstart, qend, 0, number_of_frames, number_of_frames*(1-cos(pi*frameno/number_of_frames))/2
		qcur   = quaternion.slerp(qstart, qend,  0, number_of_frames, number_of_frames*(1-cos(pi*frameno/number_of_frames))/2)
		# some funny results can be obtained with tiny numbers  (2D proteins and such) not sure where that comes from
		intermediate_view = [x if fabs(x)>0.001 else 0 for x in  list(quaternion.as_rotation_matrix(qcur).flatten())]

		for i in range(9, len(view_init)):
			update = view_init[i] + (view_last[i] - view_init[i]) * (frameno / number_of_frames)
			if fabs(update) < 0.001: update = 0
			intermediate_view.append(update)
		cmd.set_view(view2view_string(intermediate_view))
		cmd.png("frm" + str(last_frame).zfill(3), width=1920, height=1080, ray=True)
		last_frame += 1

	return last_frame

