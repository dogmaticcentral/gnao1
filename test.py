#! /usr/bin/python3
# can be run without popping the gui with pymol -cq
# though it becomes ridiculously slow for some reason
# (thus, run from pymol GUI; then )
# note that pymol pieces imports pymol.py
import os
from math import *


from utils.pymol_pieces import *
from utils.pheno_scene_views import *

# from pyquaternion import Quaternion # will not take matrix that is (numerically) non-orthogona;
import numpy as np
import quaternion # pip3 install numpy-quaternion and numba

frames_home = "/home/ivana/projects/gnao1db/movie"



def mat2quat(m):

	tr = m[0] + m[4] + m[8]

	#    m00=0 m01=1 m02=2
	#    m10=3 m11=4 m12=5
	#    m20=6 m21=7 m22=8

	if (tr > 0):
		S = sqrt(tr+1) * 2
		qw = 0.25 * S
		qx = (m[7] - m[5]) / S
		qy = (m[2] - m[6]) / S
		qz = (m[3] - m[1]) / S

	elif ((m[0] > m[4])&(m[0] > m[8])):
		S = sqrt(1 + m[0] - m[4] - m[8]) * 2
		qw = (m[7] - m[5]) / S
		qx = 0.25 * S
		qy = (m[1] + m[3]) / S
		qz = (m[2] + m[6]) / S

	elif (m[4] > m[8]):
		S = sqrt(1 + m[4] - m[0] - m[8]) * 2
		qw = (m[2] - m[6]) / S
		qx = (m[1] + m[3]) / S
		qy = 0.25 * S
		qz = (m[5] + m[7]) / S

	else:
		S = sqrt(1 + m[8] - m[0] - m[4]) * 2
		qw = (m[3] - m[1]) / S
		qx = (m[2] + m[6]) / S
		qy = (m[5] + m[7]) / S
		qz = 0.25 * S

	return [qx,qy,qz,qw]

def quat2mat(Q):

	xx = Q[0]*Q[0]
	xy = Q[0]*Q[1]
	xz = Q[0]*Q[2]
	xw = Q[0]*Q[3]
	yy = Q[1]*Q[1]
	yz = Q[1]*Q[2]
	yw = Q[1]*Q[3]
	zz = Q[2]*Q[2]
	zw = Q[2]*Q[3]

	M = [1.0 - 2*yy - 2*zz,     2*xy - 2*zw,      2*xz + 2*yw,
		 2*xy + 2*zw,       1 - 2*xx - 2*zz,      2*yz - 2*xw,
		 2*xz - 2*yw,           2*yz + 2*xw,  1 - 2*xx - 2*yy]

	return M


#----------------------------------------------------------
# Methods to perform Spherical Linear Interpolation (SLERP)
# The code has been translated to python from sources
# available at http://www.euclideanspace.com/maths/
#          algebra/realNormedAlgebra/quaternions/slerp/
#----------------------------------------------------------
def slerp(qa, qb, t):

	# Calculate angles between quaternions
	cosHalfTheta = qa[0] * qb[0] + qa[1] * qb[1] + qa[2] * qb[2] + qa[3] * qb[3]

	if cosHalfTheta < 0:
		for i in range(4):
			qb[i] = -qb[i]
		cosHalfTheta = -cosHalfTheta

	# if qa=qb or qa=-qb then theta = 0 and we can return qa
	if abs(cosHalfTheta) >= 1.0:
		return qa.copy()

	# Calculate temporary values
	halfTheta    = acos(cosHalfTheta)
	sinHalfTheta = sqrt(1.0 - cosHalfTheta*cosHalfTheta)

	if fabs(sinHalfTheta) < 0.000005: #fabs is floating point absolute
		ratioA = ratioB = 0.5
	else:
		ratioA = sin((1 - t) * halfTheta) / sinHalfTheta
		ratioB = sin(t * halfTheta) / sinHalfTheta

	return [qa[i]*ratioA + qb[i]*ratioB for i in range(4)]


def view_string2view(view_string):
	return [float(f) for f in view_string.split(",")]


def view2view_string(view):
	return ",".join([str(f) for f in view])


def view2quat0(view):
	return mat2quat(view[:9])

def view_interpolate0(view_init_str, view_last_str, number_of_frames=5, frameno_offset=0):

	view_init = view_string2view(view_init_str)
	view_last = view_string2view(view_last_str)
	# get quaternions for interpolation
	qstart = view2quat0(view_init)
	qend   = view2quat0(view_last)

	for frameno in range(1, number_of_frames + 1):
		intermediate_view = []
		qcur   = slerp(qstart, qend, frameno/number_of_frames)
		print(frameno, qcur)
		matcur = quat2mat(qcur)

		for i in range(9):
			update = matcur[i]
			#update = 0 if abs(matcur[i]) < 0.001 else matcur[i]
			intermediate_view.append(update)

		for i in range(9, len(view_init)):
			update = view_init[i] + (view_last[i] - view_init[i]) * (frameno / number_of_frames)
			# if abs(update) < 0.001: update = 0
			intermediate_view.append(update)

		# cmd.set_view(view2view_string(intermediate_view))
		# cmd.png("frm" + str(frameno_offset + frameno).zfill(3), width=1920, height=1080, ray=True)

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

		qcur   = quaternion.slerp(qstart, qend, 0, number_of_frames, frameno/number_of_frames)
		print(frameno, qcur)
		intermediate_view = list(quaternion.as_rotation_matrix(qcur).flatten())

		for i in range(9, len(view_init)):
			update = view_init[i] + (view_last[i] - view_init[i]) * (frameno / number_of_frames)
			if abs(update) < 0.001: update = 0
			intermediate_view.append(update)
		#
		# cmd.set_view(view2view_string(intermediate_view))
		# cmd.png("frm" + str(frameno_offset + frameno).zfill(3), width=1920, height=1080, ray=True)
		print(intermediate_view)


###################################
def main():

	view_interpolate(pheno_view[1], pheno_view[2], number_of_frames=5, frameno_offset=0)

	#cmd.quit()
#########################################
if __name__ == '__main__':
	main()
