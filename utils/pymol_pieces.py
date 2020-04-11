
from math import acos, sqrt, fabs, sin

from pymol import cmd
from math import pi, cos

# from pyquaternion import Quaternion # will not take matrix that is (numerically) non-orthogona;
import numpy as np
import quaternion  # pip3 install numpy-quaternion and numba



def style_lipid(lipid_selection_name):
	cmd.show("sticks",lipid_selection_name)
	cmd.color("lightblue", lipid_selection_name)
	cmd.set("stick_transparency", 0.7, lipid_selection_name)


def make_GDP(in_name, new_name):
	cmd.copy(new_name, in_name, zoom=0)
	cmd.remove("{} and not resn GDP".format(new_name))
	return


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


def clump_representation(regions, color, name, transparency=0.0, small_molecule=False):

	# show regions as clumps or blobs
	cmd.set("surface_quality", 1)
	if transparency>0:
		cmd.set("transparency", transparency, name)
	else:
		cmd.set("spec_reflect", 0.0)

	if small_molecule:
		cmd.alter("all", "b=20")
		cmd.alter("all", "q=1")
		cmd.set("gaussian_resolution", 7)
		grid_spacing = 0.2
	else:
		cmd.alter("all", "b=50")
		cmd.alter("all", "q=1")
		cmd.set("gaussian_resolution", 7)
		grid_spacing = 1.5

	idx = 0
	for region in regions:
		idx += 1
		mapname  = "map_{}_{}".format(name, idx)
		surfname = "surf_{}_{}".format(name, idx)
		cmd.map_new(mapname, "gaussian", grid_spacing, region, 2)
		cmd.isosurface(surfname, mapname)
		cmd.color(color, surfname)


def clump_cleanup(regions, name):
	for idx in range(1,len(regions)+1):
		mapname  = "map_{}_{}".format(name, idx)
		surfname = "surf_{}_{}".format(name, idx)
		cmd.delete(surfname)
		cmd.delete(mapname)


def view_string2view(view_string):
	return [float(f) for f in view_string.split(",")]


def view2view_string(view):
	return ",".join([str(f) for f in view])


def view2quat(view):

	array  = np.array(view[:9])
	matrix = np.asmatrix(array.reshape((3,3)))
	return quaternion.from_rotation_matrix(matrix)



def view_interpolate(view_init_str, view_last_str, base_name, number_of_frames=5, frameno_offset=0):

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
		intermediate_view = [x if fabs(x)>0.001 else 0 for x in list(quaternion.as_rotation_matrix(qcur).flatten())]

		for i in range(9, len(view_init)):
			# pymol is sitll on python2; in python2 3/5 = 0; in python3 it is 0.6
			update = view_init[i] + (view_last[i] - view_init[i]) * (float(frameno) / number_of_frames)
			if fabs(update) < 0.001: update = 0
			intermediate_view.append(update)
		cmd.set_view(view2view_string(intermediate_view))
		cmd.png(base_name + str(last_frame).zfill(3), width=1920, height=1080, ray=True)
		last_frame += 1

	return last_frame

'''
   Note that when homogenous is zero, the input matrix is NOT a
    standard homogenous 4x4 transformation matrix.  Instead it is
    something PyMOL-specific which consists of the following:

    1) a 3x3 matrix containing the rotation in the upper-left quadrant

    2) a 1x3 translation to be applied *before* rotation in the bottom row
        (matrix[12],matrix[13],matrix[14]).

    3) a 3x1 translation to be applied *after* rotation in the right-hand
        column (matrix[3],matrix[7],matrix[11])

    In other words, if the matrix is:

    [  m0  m1  m2  m3 \\
       m4  m5  m6  m7 \\
       m8  m9 m10 m11 \\
      m12 m13 m14 m15 ] 

    Atoms will be transformed as follows

    Y = M X

    y0 = m0*(x0+m12) + m1*(x1+m13) +  m2*(x2+m14) + m3 \\
    y1 = m4*(x0+m12) + m5*(x1+m13) +  m6*(x2+m14) + m7 \\
    y2 = m8*(x0+m12) + m9*(x1+m13) + m10*(x2+m14) + m11 

'''
def tfm2quat(tfm):
	array  = np.array(tfm[:3] + tfm[4:7] + tfm[8:11])
	matrix = np.asmatrix(array.reshape((3,3)))
	return quaternion.from_rotation_matrix(matrix)


def intermediate_tfm(source_tfm, target_tfm, number_of_frames, frameno):
	tfm = [0]*16
	qstart = tfm2quat(source_tfm)
	qend = tfm2quat(target_tfm)
	qcur = quaternion.slerp(qstart, qend, 0, number_of_frames, number_of_frames*(1-cos(pi*frameno/number_of_frames))/2)
	intermediate_rot = list(quaternion.as_rotation_matrix(qcur).flatten())
	tfm[0:3]  = intermediate_rot[:3]
	tfm[4:7]  = intermediate_rot[3:6]
	tfm[8:11] = intermediate_rot[6:9]
	for index in [3,7,11]:
		tfm[index] += (target_tfm[index] - source_tfm[index])*float(frameno)/number_of_frames
	return tfm


def object_tfm_interpolate(view, object_properties, base_name, number_of_frames=15, frameno_offset=0):

	# I could not find the way to move the surface, if there's one provided
	# so I'm moving the object and re-creating the mesh

	for objnm in object_properties.keys():
		cmd.hide("everything", objnm)
		clump_cleanup([objnm], objnm)

	last_frame = frameno_offset
	for frameno in range(1, number_of_frames+1):
		tmpnames = []
		for objnm, [source_tfm, target_tfm, reverse, color, small_molecule] in object_properties.items():
			frame = frameno
			if reverse: frame = number_of_frames  - frame
			tfm = intermediate_tfm(source_tfm, target_tfm, number_of_frames, frame)
			tmpnm = objnm + str(frameno)
			cmd.copy(tmpnm, objnm, zoom=0)
			cmd.transform_selection(tmpnm, tfm)
			clump_representation([tmpnm], color, tmpnm, small_molecule=small_molecule)
			tmpnames.append(tmpnm)
		# some of the ops above move camera (why the f do you have to move the camera to create a new rep?)
		cmd.set_view(view)
		cmd.png(base_name + str(last_frame).zfill(3), width=1920, height=1080, ray=True)
		for tmpnm in tmpnames:
			clump_cleanup([tmpnm], tmpnm)
			cmd.remove(tmpnm)
		last_frame += 1



	return last_frame

