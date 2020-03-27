
from pymol import cmd

def load_structures(structure_home, structures):

	if "gnao" in structures:
		cmd.load("{}/{}".format(structure_home,"gnao1_after_1azsC.to_3sn6A.pdb"),"gnao")
	cmd.hide("all")
	return




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
	qm=[]

	#Calculate angles between quaternions
	cosHalfTheta = qa[0] * qb[0] + qa[1] * qb[1] + qa[2] * qb[2] + qa[3] * qb[3]
	#if qa=qb or qa=-qb then theta = 0 and we can return qa

	if (cosHalfTheta < 0):
		for i in range(4):
			qb[i] = -qb[i];
		cosHalfTheta = -cosHalfTheta

	if (abs(cosHalfTheta) >= 1.0):
		for i in range(4):
			qm.insert(i,qa[i])
		return qm

	#Calculate temporary values
	halfTheta = acos(cosHalfTheta)
	sinHalfTheta = sqrt(1.0 - cosHalfTheta*cosHalfTheta)

	if (fabs(sinHalfTheta) < 0.000005):
		#fabs is floating point absolute
		for i in range(4):
			qm.insert(i, qa[i] * 0.5 + qb[i] * 0.5)
		return qm

	ratioA = sin((1 - t) * halfTheta) / sinHalfTheta
	ratioB = sin(t * halfTheta) / sinHalfTheta
	#calculate Quaternion
	for i in range(4):
		qm.insert(i, qa[i] * ratioA + qb[i] * ratioB)
	return qm

# Changes position of the camera
# @frames: number of frames to animate
# @x: view matrix of the starting frame
# @y: view matrix of the target frame
def changeView (frames, x, y):

	setView(x)
	ray()

	# get quaternions for interpolation
	# between starting and target scenes
	qstart = mat2quat(x[:9])
	qend   = mat2quat(y[:9])

	for frame in range(1, frames+1):
		n = []
		qcur   = slerp(qstart, qend, frame/frames)
		matcur = quat2mat(qcur)

		for i in range(9):
			update = 0 if abs(matcur[i]) < 0.001 else matcur[i]
			n.append(update)

		for i in range(9,len(x)):
			update = x[i] + (y[i] - x[i]) * (frame/frames)
			if abs(update) < 0.001: update = 0
			n.append(update)

		setView(n)

		ray()

