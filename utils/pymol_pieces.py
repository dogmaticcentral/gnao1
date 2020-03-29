from math import acos, sqrt, fabs, sin

from pymol import cmd


# from pyquaternion import Quaternion # will not take matrix that is (numerically) non-orthogona;
import numpy as np
import quaternion  # pip3 install numpy-quaternion and numba


def pymol_chdir(outdir):
	cmd.cd(outdir)


def load_structures(structure_home, structures):

	if "gnao" in structures:
		cmd.load("{}/{}".format(structure_home, "gnao1_after_1azsC.to_3sn6A.pdb"), "gnao")
	if "substrate" in structures:
		cmd.load("{}/{}".format(structure_home, "3c7kA.GDP.to_3sn6A.pdb"), "substrate")
	cmd.hide("all")
	return


def residue_color(main_object, res_id, rgb_color_list):
	color_name  = "cnm_{}_{}".format(main_object, res_id)
	cmd.set_color(color_name, rgb_color_list)
	cmd.color(color_name, "{} and resi {}".format(main_object, res_id))


def clump_representation(regions, color, name):
	# show regions as clumps or blobs
	cmd.set("surface_quality", 1)
	cmd.alter("all", "b=50")
	cmd.alter("all", "q=1")
	cmd.set("gaussian_resolution", 5)

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

	for frameno in range(1, number_of_frames + 1):

		qcur   = quaternion.slerp(qstart, qend, 0, number_of_frames, frameno)
		intermediate_view = list(quaternion.as_rotation_matrix(qcur).flatten())

		for i in range(9, len(view_init)):
			update = view_init[i] + (view_last[i] - view_init[i]) * (frameno / number_of_frames)
			if abs(update) < 0.001: update = 0
			intermediate_view.append(update)

		cmd.set_view(view2view_string(intermediate_view))
		cmd.png("frm" + str(frameno_offset + frameno).zfill(3), width=1920, height=1080, ray=True)
