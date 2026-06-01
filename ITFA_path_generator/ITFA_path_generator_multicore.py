#!/usr/bin/env python3
import os
import sys
import os.path
import concurrent.futures
import time
from shutil import copyfile

cores=os.cpu_count()
multicores=False
if cores>1:
	multicores=True

def order_search(dir_name):
	j=dir_name.split(os.sep)
	k=j[len(j)-1].split()
	l=k[0].split('_')
	m=k[len(k)-1].split('_')
	if str(l[0]).isdigit():
		return l[0]
	elif str(m[len(m)-1]).isdigit():
		return m[len(m)-1]
	else:
		print("Warning! Files will be copied in a random order.")
		return -1

def file_search(file_end):
	main=os.getcwd()
	path=[]
	temp_path=[]
	for dirpath, dirnames, filenames in os.walk(main):
		for filename in [f for f in filenames if f.endswith(file_end)]:
			temp_path.append([dirpath, filename])
	temp_path.sort(key=lambda x: int(order_search(x[0])))
	for i in temp_path:
		path.append(str(os.path.join(i[0], i[1])))
	temp_path=[]
	return path

def inf_file_gen(filename, path, list_length, min_e, max_e):
	output=open("inf.txt","w")
	print("data input options, data normalization, etc.:", file=output)
	print(path+str(os.sep)+" ! Data files directory", file=output)
	print(str(list_length)+"                                 ! Total number of files n (min 2 max 9999)", file=output)
	print("1                                  ! First file number", file=output)
	print("n                                  ! Prefix ? (if no, next line is ignored)", file=output)
	print("Prefix                             ! Prefix (max 50 char)", file=output)
	print("n                                  ! Suffix ? (if no, next line is ignored)", file=output)
	print("suffix                             ! Suffix (max 50 char)", file=output)
	print(filename+"                                ! Extension (without dot)", file=output)
	print("{:<10} {}".format(str(round(min_e,5)),"                        ! Lower limit of abscissa"),file=output)
	print("{:<10} {}".format(str(round(max_e,5)),"                        ! Upper limit of abscissa"),file=output)
	print("1                                  ! Read in each N-th data point; N = ", file=output)
	print("n                                  ! Add noise to the data ?", file=output)
	print("0.0030                             ! Noise amplituden                                  ", file=output)
	if filename=="chi":
		print("y                                  ! y-weighting (Y=y*x^weighting) ?", file=output)
	else:
		print("n                                  ! y-weighting (Y=y*x^weighting) ?", file=output)
	print("3                                  ! weighting value", file=output)
	print("n                                  ! normalize spectra (y/sqrt(variance(y))) ? ", file=output)
	print("n                                  ! center and normalize spectra (y-mean_y)/data_points/sqrt(variance(y)) ?", file=output)
	print("PRINCIPAL COMPONENT ANALYSIS (PCA) and VARIMAX:", file=output)
	print("{:<10}".format(str(list_length)),"                        ! Number of expected components (PCA)", file=output)
	if filename=="chi":
		print("y                                  ! Fourier transf. of eigenvectors, exp. data and abstr. reproductions ? (PCA)", file=output)
		print("-1                                 ! Graphical y-offset for Fourier transforms (PCA) ", file=output)
		print("-10                                ! Graphical y-offset for abstract reproductions (PCA)", file=output)
	else:
		print("n                                  ! Fourier transf. of eigenvectors, exp. data and abstr. reproductions ? (PCA)", file=output)
		print("-1                                 ! Graphical y-offset for Fourier transforms (PCA) ", file=output)
		print("-5                                 ! Graphical y-offset for abstract reproductions (PCA)", file=output)
	print("1                                  ! Use the abscissa scale from file number n for interpolation on unique x axis (PCA)", file=output)
	print("n                                  !equidistant x-axis? (y/n)", file=output)
	print("ITERATIVE TARGET TEST (ITT):", file=output)
	print("ITT.txt                            ! Concentration test vector filename (ITT)", file=output)
	print("n                                  ! Suppress double concentration maxima ? (ITT)", file=output)
	print("y                                  ! Fourier transform of the component spectra ? (ITT)", file=output)
	print("0.0                                ! Graphical y-offset the component spectra (ITT)", file=output)
	print("n                                  ! adjust sum of relative component concentrations to 1 ? (ITT)", file=output)
	print("n                                  ! de-normalize ? (ITT) (if the spectra were normalized with (y/sqrt(variance(y)), see above)", file=output)
	print("0                                  ! y-weighting of the component spectra (ITT)", file=output)
	output.close()
	return 0

def file_copy(arg_list):
	edge=-1
	filename=arg_list[0]
	source_line=arg_list[1]
	destination_line=arg_list[2]
	file_output=open(destination_line,"w")
	if filename=="nrm":
		temp=0
	with open (source_line,"r") as readl:
		for line in readl:
			if filename=="nrm":
				temp=temp+1
				if temp==3:
					element=line.split()
					e0=element[4].split(',')
					edge=float(e0[0])
			if filename=="chi" and line[0]!="#":
				element=line.split()
				edge=float(element[0])
			if line[0]!="#":
				print(line, end="", file=file_output)
	readl.close()
	file_output.close()
	return edge

def copy_files(filename, path_list):
	pwd=os.getcwd()
	os.mkdir(filename)
	os.chdir(filename)
	output=open("old_"+filename+"_pathes_ref.dat", "w")
	step=0
	global multicores
	edge_min=1000000
	edge_max=0
	print(str(pwd)+str(os.sep), file=output)
	multicore_list=[]
	for files in path_list:
		step=step+1
		source_line=files
		destination_line=str(step)+"."+filename
		multicore_list.append([filename, source_line, destination_line])
		print(source_line+"\t"+"\t"+destination_line, file=output)
	output.close()
	results=[]
	start=0.0
	finish=0.0
	start=time.perf_counter()
	if multicores:
		with concurrent.futures.ProcessPoolExecutor() as executor: 
			results=executor.map(file_copy, multicore_list)
	else:
		for single_core in multicore_list:
			results.append(file_copy(single_core))
	finish=time.perf_counter()
	if multicores:
		print(f'multicore finished in {round(finish-start, 12)} seconds')
	else:
		print(f'single core finished in {round(finish-start, 12)} seconds')
	if filename=="nrm":
		for res in results:
			edge_max=max(edge_max, res)
			edge_min=min(edge_min, res)	
		edge_difference=(edge_max-edge_min)/2
		if(edge_difference<=0.001):
			edge_max=edge_max-edge_difference+0.05
			edge_min=edge_min+edge_difference-0.05
		elif(edge_difference<=0.01):
			edge_max=edge_max-edge_difference+0.1
			edge_min=edge_min+edge_difference-0.1
		elif(edge_difference<=0.1):
			edge_max=edge_max-edge_difference+0.15
			edge_min=edge_min+edge_difference-0.15
		else:
			edge_max=edge_max-edge_difference+0.2
			edge_min=edge_min+edge_difference-0.2
		inf_file_gen(filename, str(os.getcwd()), step, edge_min, edge_max)
	if filename=="chi":
		for res in results:
			edge_min=min(edge_min, res)
		inf_file_gen(filename, str(os.getcwd()), step, 2.0, edge_min)
	os.chdir(pwd)
	return 0

try:
	os.mkdir("XAFS_ITFA")
except OSError:
	print("Error! Please remove/rename old \"XAFS_ITFA\" folder prior to running this script.")
	input()
	sys.exit(-1)
else:
	chi=file_search("chi")
	nrm=file_search("nrm")
	if len(chi)==0:
		print("Error! There is no *.chi files in subfolders.")
		input()
		sys.exit(-2)
	if len(nrm)==0:
		print("Error! There is no *.nrm files in subfolders.")
		input()
		sys.exit(-3)
	if len(chi)!=len(nrm):
		print("Warning! There is not equal amount of *.chi and *.nrm files in subfolders.")
		print("Please press Enter button if you want to proceed. Otherwise, hit another key+Enter.")
		button=input()
		if (button!=""):
			sys.exit(-4)
	os.chdir("XAFS_ITFA")
	copy_files("chi", chi)
	copy_files("nrm", nrm)
	print("Program has finished successfully.")