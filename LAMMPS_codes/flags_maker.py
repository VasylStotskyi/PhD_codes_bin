import time
import math
#import concurrent.futures
#import multiprocessing as mp
import matplotlib.pyplot as plt
#from memory_profiler import profile

md_water_list=[]



def bond_eval(temp_water_list1, temp_water_list2, x, y, z):   #creates a bonding list from initial coordinates. It's correct as long as water molecules doesn't dissociate or change their position in dump file (sequential number, not coordinates)
    temp_water_list3=[]
    temp_water_list4=[]
    bmax=1.2
    bmin=0.8
    bmax2=bmax**2
    bond_length=0.0
    bx=0.0
    by=0.0
    bz=0.0
    x2=x**2
    y2=y**2
    z2=z**2
    for i in range(0,len(temp_water_list1)):
        temp_water_list4.append(i)
        for j in range(0, len(temp_water_list2)):
            bx=min(abs(temp_water_list1[i][0]-temp_water_list2[j][0]),abs(temp_water_list1[i][0]-temp_water_list2[j][0]+x))
            by=min(abs(temp_water_list1[i][1]-temp_water_list2[j][1]),abs(temp_water_list1[i][1]-temp_water_list2[j][1]+y))
            bz=min(abs(temp_water_list1[i][2]-temp_water_list2[j][2]),abs(temp_water_list1[i][2]-temp_water_list2[j][2]+z))
            bx=bx**2
            by=by**2
            bz=bz**2
            if bx<=bmax2 and by<=bmax2 and bz<=bmax2:
                bond_length=bx+by+bz
                bond_length=math.sqrt(bond_length)
                if bond_length > bmin and bond_length <= bmax:
                    temp_water_list4.append(j)
        if len(temp_water_list4)>3:
            print("Error!")
            print(bond_length)
        temp_water_list3.append(temp_water_list4)
        temp_water_list4=[]
    return temp_water_list3

def read_bonds(filename):
    read_bonds=False
    read_atoms=False
    temp_bond_list=[]
    temp_atoms_list=[]
    with open(filename, "r") as inp:
        for line in inp:
            element=line.split()
            if len(element)==1 and element[0]=="Atoms":
                read_atoms=True
                line=next(inp)
                continue
            elif read_atoms==True:
                if len(element)==0:
                    read_atoms=False
                else:
                    temp_atoms_list.append([float(element[4]), float(element[5]), float(element[6]), int(element[7]), int(element[8]), int(element[9])])
            elif len(element)==1 and element[0]=="Bonds":
                read_bonds=True
                line=next(inp)
                continue
            elif read_bonds==True:
                if len(element)==0:
                    read_bonds=False
                    break
                else:
                    temp_bond_list.append([int(element[2]), int(element[3])])
            else:
                continue
    inp.close()
    for i in temp_bond_list:
        dx=temp_atoms_list[i[0]-1][0]-temp_atoms_list[i[1]-1][0]
        dy=temp_atoms_list[i[0]-1][1]-temp_atoms_list[i[1]-1][1]
        dz=temp_atoms_list[i[0]-1][2]-temp_atoms_list[i[1]-1][2]
        distance=dx**2+dy**2+dz**2
        distance=abs(math.sqrt(distance))
        if distance>3.5:
            temp_atoms_list[i[1]-1][3]=dx
            temp_atoms_list[i[1]-1][4]=dy
            temp_atoms_list[i[1]-1][5]=dz
    for j in temp_atoms_list:
        print(j)
    return 0


bond_filename="output.lammps.dat"
bond_list=[]
bond_list=read_bonds(bond_filename)
"""
start=0.0
finish=0.0
start=time.perf_counter()
output(density_calculation(md_water_list,"z"), "density_list.dat")
finish=time.perf_counter()
print(f'density calculator finished in {round(finish-start, 4)} seconds')
"""