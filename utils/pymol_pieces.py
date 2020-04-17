
from math import acos, sqrt, fabs, sin

from pymol import cmd
from math import pi, cos

# from pyquaternion import Quaternion # will not take matrix that is (numerically) non-orthogona;
import numpy as np
import quaternion  # pip3 install numpy-quaternion and numba

#from utils.pymol_constants import * # this is pymol3
from pymol_constants import *

def style_lipid(lipid_selection_name):
	cmd.show("sticks",lipid_selection_name)
	cmd.color( mol_color.get(lipid_selection_name, "white"), lipid_selection_name)
	cmd.set("stick_transparency", 0.7, lipid_selection_name)


def style_substrate(substrate_selection_name, color):
	cmd.show("spheres", substrate_selection_name)
	cmd.color(color, substrate_selection_name)



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
		cmd.load("{}/{}".format(structure_home, structure_filename[structure]), object=structure)
		cmd.hide("everything", structure)

	if "AC" in structures:
		cmd.remove("AC and chain B")

	if "GPCR" in structures:
		cmd.remove("resi 343-1200 and GPCR")

	return


def extract_state_to_object(morph, state, new_object):

	number_of_states = cmd.count_states(morph)
	if number_of_states<state:
		print("in extract_state_to_object() requested state " +
		      "{}  > number of states ({}) available in {}".format(state, number_of_states, morph))
		exit()
	if morph==new_object:
		print("in extract_state_to_object() please picked different " +
		      "name for extraction (currently bouth '{}') ".format(morph))
		exit()

	# after this I will have each state available as its own object, called "{morphname}_000x" (padded to length of 4)
	# x starts with 1 --> not clear - that might depend on the numbering in the input file
	cmd.split_states(morph)
	state_name = morph + "_" + str(state).zfill(4)
	cmd.copy(new_object, state_name)
	for statenum in range(1, number_of_states+1):
		state_name = morph + "_" + str(statenum).zfill(4)
		cmd.delete(state_name)
	cmd.delete(morph)

	return


def residue_color(main_object, res_id, rgb_color_list):
	color_name  = "cnm_{}_{}".format(main_object, res_id)
	cmd.set_color(color_name, rgb_color_list)
	cmd.color(color_name, "{} and resi {}".format(main_object, res_id))


def clump_representation(regions, color, name, transparency=-1.0, small_molecule=False):

	# show regions as clumps or blobs
	cmd.set("surface_quality", 1)
	if transparency<=0:
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
		if transparency>0:
			cmd.set("transparency", transparency, surfname)

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


def intermediate_view(view_init, view_last, qstart, qend, number_of_frames, frameno):
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

	return intermediate_view


def view_interpolate(view_init_str, view_last_str, base_name, number_of_frames=5, frameno_offset=0):

	view_init = view_string2view(view_init_str)
	view_last = view_string2view(view_last_str)
	# get quaternions for interpolation
	qstart = view2quat(view_init)
	qend   = view2quat(view_last)

	last_frame = frameno_offset
	for frameno in range(1, number_of_frames+1):
		intrm_view = intermediate_view(view_init, view_last, qstart, qend, number_of_frames, frameno)
		cmd.set_view(view2view_string(intrm_view))
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


def object_tfm_interpolate(object_properties, number_of_frames, frameno):
	tmpnames = []

	for objnm, properties in object_properties.items():
		[source_tfm, target_tfm, reverse, color, small_molecule] = properties[:5]
		transparency_range = properties[5] if len(properties)>5 else None
		frame = frameno
		if reverse: frame = number_of_frames  - frame
		tfm = intermediate_tfm(source_tfm, target_tfm, number_of_frames, frame)
		tmpnm = objnm + str(frameno)
		cmd.copy(tmpnm, objnm, zoom=0)
		cmd.transform_selection(tmpnm, tfm)
		transparency = -1
		if transparency_range:
			transparency = transparency_range[0] \
			               + float(frameno)/number_of_frames*(transparency_range[1]-transparency_range[0])
		if small_molecule:
			style_substrate(tmpnm, mol_color[objnm])
		else:
			clump_representation([tmpnm], color, tmpnm, small_molecule=small_molecule, transparency=transparency)
		tmpnames.append(tmpnm)
	return tmpnames


def object_cleanup(tmpnames):
	for tmpnm in tmpnames:
		clump_cleanup([tmpnm], tmpnm)
		cmd.remove(tmpnm)



def morph_movie(morph, view, color, base_name, frameno_offset=0, morph_reverse=False):
	number_of_states = cmd.count_states(morph)
	# after this I will have each state available as its own object, called "{morphname}_000x" (padded to length of 4)
	# x starts with 1 --> not clear - that might depend on the numbering in the input file
	cmd.split_states(morph)
	last_frame = frameno_offset
	for statenum in range(0, number_of_states+1):
		if morph_reverse:
			statenum = max(1, cmd.count_states(morph) - statenum)
		else:
			statenum = min(statenum+1, cmd.count_states(morph))
		state_name = morph + "_" + str(statenum).zfill(4)
		clump_representation([state_name], color, state_name),
		cmd.set_view(view)
		cmd.png(base_name + str(last_frame).zfill(3), width=1920, height=1080, ray=True)
		clump_cleanup([state_name], state_name)
		cmd.remove(state_name)
		last_frame += 1

	return last_frame



def scene_interpolate(view_init_str, object_properties, base_name,
                      number_of_frames=15, frameno_offset=0,
                      view_last_str=None, morph_properties=None):

	if morph_properties:
		morph_lengths = set()
		for morph_name in morph_properties.keys():
			morph_lengths.add(cmd.count_states(morph_name))
			cmd.split_states(morph_name)
		if len(morph_lengths)>1:
			print("morphs are all expected to be of the same length; found", morph_lengths)
			exit()
		if len(morph_lengths)==0: # can this happen? morph_properties should be False in thate casem shouldn't it
			morph_properties = None
		else:
			number_of_frames = morph_lengths.pop()

	# I could not find the way to move the surface, if there's one provided
	# so I'm moving the object and re-creating the mesh

	view_init = view_string2view(view_init_str) if view_last_str else None
	view_last = view_string2view(view_last_str) if view_last_str else None
	# get quaternions for interpolation
	qstart = view2quat(view_init) if view_last_str else None
	qend   = view2quat(view_last) if view_last_str else None

	for objnm in object_properties.keys():
		cmd.hide("everything", objnm)
		clump_cleanup([objnm], objnm)

	last_frame = frameno_offset
	for frameno in range(0, number_of_frames+1):
		# object position interpolation
		tmpnames = object_tfm_interpolate(object_properties, number_of_frames, frameno)
		# view interpolation, if called for
		if not view_last_str:
			view = view_init_str
		else:
			view  = view2view_string(intermediate_view(view_init, view_last, qstart, qend, number_of_frames, frameno))
		# morph
		if morph_properties:
			tmp_obj_props = {}
			for morph_name, props in morph_properties.items():
				[morph_color, morph_reverse, tfm_from, tf_to, tfm_reverse] = props
				# change shape
				if morph_reverse:
					stateno = max(1, cmd.count_states(morph_name) - frameno)
				else:
					stateno = min(frameno+1, cmd.count_states(morph_name))
				morph_state_name = morph_name + "_" + str(stateno).zfill(4)
				# [identity_tfm, gbg_tfm, tfm_reverse, color, small_molecule rep]
				tmp_obj_props[morph_state_name] = [tfm_from, tf_to, tfm_reverse, morph_color, False]
			# reposition all morphs and show them as clumps
			tmpnames.extend(object_tfm_interpolate(tmp_obj_props, number_of_frames, frameno))
		# some of the ops above move camera (why the f do you have to move the camera to create a new rep?)
		cmd.set_view(view)
		cmd.png(base_name + str(last_frame).zfill(3), width=1920, height=1080, ray=True)
		object_cleanup(tmpnames)

		last_frame += 1

	return last_frame

