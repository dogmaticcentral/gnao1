import os, subprocess

def check_deps(deps):
	for dependency in deps:
		if not os.path.exists(dependency):
			print(dependency, "not found")
			exit()
		if not os.access(dependency, os.X_OK):
			print(dependency, "executable")
			exit()

def cleanup(rootname):
	for fnm in os.listdir():
		if rootname not in fnm: continue
		for ext in [".gplt", ".gplt_out", ".bngl", ".bngl_out", ".gdat", ".cdat", ".net", ".dat"]:
			if fnm[-len(ext):] != ext: continue
			os.remove(fnm)
			break

def run_bngl(bngl, bngl_input):
	bngl_out = bngl_input.replace(".bngl", ".bngl_out")
	#  2>&1 is stderr to stdout
	status = subprocess.call(["bash", "-c", f"{bngl} {bngl_input} > {bngl_out} 2>&1  "])
	if status == 0:
		return "ok"
	return "failure"

def run_gnuplot(gnuplot, gnuplot_input):
	gplt_out = gnuplot_input.replace(".gplt", ".gplt_out")
	subprocess.call(["bash", "-c", f"{gnuplot} {gnuplot_input} > {gplt_out}"])
	return

