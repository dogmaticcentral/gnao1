
def set_gnuplot_outfile(rootname, svg=False):
	retstr = ""
	if svg:
		retstr += "# svg\n"
		retstr += "set terminal svg size 410,250 fname 'Verdana, Helvetica, Arial, sans-serif' fsize '9' rounded dashed \n"
		retstr += f"set output '{rootname}.svg' \n"
		retstr += "set size ratio 1  \n"
	else:
		retstr += "# png\n"
		retstr += "set terminal pngcairo size 615,500 enhanced font 'Verdana,18'\n"
		retstr += f"set output '{rootname}.png' \n"
	return retstr
