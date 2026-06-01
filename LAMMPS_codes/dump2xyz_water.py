import time
import math
#import concurrent.futures
#import multiprocessing as mp
import matplotlib.pyplot as plt
#from memory_profiler import profile

md_water_list=[]

def coordinates_collect_4(atom_type, input_filename, bond_filename):   #1.5 times slower than coordinates_collect_3, because also calculates ceters of mass of water molecule
    step=0
    counter=-1.0
    temp_md_water_list=[]
    with open (input_filename, "r") as inp:
        step_bool=False
        atom_bool=False
        for line in inp:
            element=line.split()
            if (step_bool==True):    #reads number of atomic coordinates
                step=int(element[0])
                counter=counter+1.0
                line=next(inp)
                element=line.split()
                step_bool=False
            if ((len(element) == 4) and (element[1]=="NUMBER")):   #activates amount of atoms reading
                step_bool=True
                continue
            if counter<1:
                continue
            if (atom_bool==True) and (step>0):   #reads O and H coordinates
                for i in range(1, step):
                    if (element[1]==atom_type):
                        temp_md_water_list.append([str("O"),float(element[2]), float(element[3]), float(element[4])])
                    if (element[1]=="2"):
                        temp_md_water_list.append([str("H"),float(element[2]), float(element[3]), float(element[4])])
                    line=next(inp)
                    element=line.split()
                if (element[1]==atom_type):
                    temp_md_water_list.app0end([str("O"), float(element[2]), float(element[3]), float(element[4])])
                if (element[1]=="2"):
                    temp_md_water_list.append([str("H"), float(element[2]), float(element[3]), float(element[4])])
                atom_bool=False
                md_water_list.extend(temp_md_water_list)  #calculates centers of mass of water molecules and saves them in list; potentially can be parallelized
                temp_md_water_list=[]
                continue
            if len(element) == 7 and element[1]=="ATOMS":   #activates atomic coordinates reading
                atom_bool=True
                continue
        output(md_water_list, "water.dat")
    inp.close()
    return md_water_list

def output(output_list, output_filename):   #prints any list with tabs between elements of list
    outputfile=open(output_filename, "w")
    step=0
    for line in output_list:
        length=0
        length=len(line)
        sentence=""
        if step%480==0:
            print(str(480), file=outputfile)        
            print(str(step//480), file=outputfile)
        for i in range(length):
            if i<length-1:
                sentence=sentence+str(line[i])+"\t"
            else:
                sentence=sentence+str(line[i])
        print(str(sentence), file=outputfile)
        step=step+1
    outputfile.close()

start=0.0
finish=0.0
start=time.perf_counter()
coordinates_collect_4("1","dump.lammps.dat", "output.lammps.dat")
finish=time.perf_counter()
print(f'coordinates_collect_4 finished in {round(finish-start, 4)} seconds')

