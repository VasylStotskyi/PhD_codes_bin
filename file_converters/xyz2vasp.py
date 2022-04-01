import math
from os import confstr
def coordinates_collect_and_sort(input_filename):
	step=0
	xyz_coordinates_list=[]
	with open (input_filename) as inp:
		for line in inp:
			element=line.split()
			step+=1
			if step>2 and len(element)==4:
				xyz_coordinates_list.append([element[0], element[1], element[2], element[3]])
	inp.close()
	xyz_coordinates_list.sort()
	for i in xyz_coordinates_list:
		print (i)
	return xyz_coordinates_list

def atom_count(coordinates_list):
	open_loop=True
	atom_type=""
	step=0
	atoms_count_list=[]
	for atom in coordinates_list:
		if open_loop:
			atom_type=atom[0]
			open_loop=False
			step=1
		elif (atom[0]==atom_type):
			step+=1
		else:
			atoms_count_list.append([atom_type, str(step)])
			atom_type=atom[0]
			step=1
	atoms_count_list.append([atom_type, str(step)])
	return atoms_count_list

def atom_count_swap(atom_count_list):
	temp1=[]
	temp2=[]
	for i in atom_count_list:
		temp1.append(i[0])
		temp2.append(i[1])
	atom_count_list=[]
	atom_count_list.append(temp1)
	atom_count_list.append(temp2)
	for i in atom_count_list:
		print(i)
	return atom_count_list

def vasp_coordinates_print_out(input_filename, sorted_coordinates):
	output_filename=input_filename[:-4]+".vasp"
	print (output_filename)
	atoms_number=atom_count_swap(atom_count(coordinates_collect_and_sort(input_filename)))
	with open(output_filename, "w") as output:
		for i in atoms_number:
			print("\t", end="", file=output)
			for j in i:
				print(j, end="\t", file=output)
			print("\n", end="", file=output)
		print("Cartesian", file=output)
		for atom in sorted_coordinates:
			print(" "+atom[1]+"\t"+atom[2]+"\t"+atom[3], file=output)
	output.close()
	return 0
	
def xyz2vasp(input_filename):
	vasp_coordinates_print_out(input_filename, coordinates_collect_and_sort(input_filename))
	return 0

def rotation_matrix(a, b, c, alpha, beta, gamma):
	R11=math.cos(alpha)*math.cos(beta)
	R12=math.cos(alpha)*math.sin(beta)*math.sin(gamma)-math.sin(alpha)*math.cos(gamma)
	R13=math.cos(alpha)*math.sin(beta)*math.cos(gamma)+math.sin(alpha)*math.sin(gamma)
	R21=math.sin(alpha)*math.cos(beta)
	R22=math.sin(alpha)*math.sin(beta)*math.sin(gamma)-math.cos(alpha)*math.cos(gamma)
	R23=math.sin(alpha)*math.sin(beta)*math.cos(gamma)+math.cos(alpha)*math.sin(gamma)
	R31=-math.sin(beta)
	R32=math.cos(beta)*math.sin(gamma)
	R33=math.cos(beta)*math.cos(gamma)
	print(str(R11)+"\t"+str(R12)+"\t"+str(R13))
	print(str(R21)+"\t"+str(R22)+"\t"+str(R23))
	print(str(R31)+"\t"+str(R32)+"\t"+str(R33))


#rotation_matrix(26.45, 27.519, 9.46, 90.46, 98.68,  90.09)

xyz2vasp("Saponite-100_final.xyz")


