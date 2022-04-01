import time
import math
#import concurrent.futures
#import multiprocessing as mp
import matplotlib.pyplot as plt
#from memory_profiler import profile

md_water_list=[]
def coordinates_collect_4(input_filename, alph_list):   #1.5 times slower than coordinates_collect_3, because also calculates ceters of mass of water molecule
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
            if counter<0:
                continue
            if (atom_bool==True) and (step>0):   #reads O and H coordinates
                for i in range(1, step):
                    temp_md_water_list.append([alph_filter(str(element[1]), alph_list),float(element[2]), float(element[3]), float(element[4])])
                    line=next(inp)
                    element=line.split()
                temp_md_water_list.append([alph_filter(str(element[1]), alph_list), float(element[2]), float(element[3]), float(element[4])])
                atom_bool=False
                md_water_list.append(temp_md_water_list)
                temp_md_water_list=[]
                continue
            if len(element) == 7 and element[1]=="ATOMS":   #activates atomic coordinates reading
                atom_bool=True
                continue
    inp.close()
    return md_water_list

def alph_filter(element, alphabeta_list):
    element_name="Unknown"
    for alpha in alphabeta_list:
        if element==alpha[0]:
            element_name=alpha[1]
            break
    return element_name


def alphabet_collect(input_filename):
    alphabet_list=[]
    with open (input_filename, "r") as inp:
        alphabet_bool=False
        for line in inp:
            element=line.split()
            if alphabet_bool==True and len(element)==4:
                alphabet_list.append([str(element[0]),str(element[3])])
            if len(element) == 1 and element[0]=="Masses":
                alphabet_bool=True
                line=next(inp)
            if len(element) == 2 and element[0]=="Pair"and element[1]=="Coeffs":
                break
    return alphabet_list

def alphabet_collect_element(input_filename):
    alphabet_list=[]
    with open (input_filename, "r") as inp:
        alphabet_bool=False
        for line in inp:
            element=line.split()
            if alphabet_bool==True and len(element)==4:
                if str(element[3][0])=="H":
                    alphabet_list.append([str(element[0]),"H"])
                    continue
                if str(element[3][0])=="O":
                    alphabet_list.append([str(element[0]),"O"])
                    continue
                if str(element[3][0])=="S":
                    alphabet_list.append([str(element[0]),"Si"])
                    continue
                if str(element[3][0])=="M":
                    alphabet_list.append([str(element[0]),"Mg"])
                    continue
                if str(element[3][0])=="Lu":
                    alphabet_list.append([str(element[0]),"Lu"])
                    continue
                if str(element[3][0])=="N":
                    alphabet_list.append([str(element[0]),"Ni"])
                    continue
            if len(element) == 1 and element[0]=="Masses":
                alphabet_bool=True
                line=next(inp)
            if len(element) == 2 and element[0]=="Pair"and element[1]=="Coeffs":
                break
    return alphabet_list

def output(output_list, output_filename):   #prints any list with tabs between elements of list
    outputfile=open(output_filename, "w")
    for i in range(len(output_list)):
        print(" "+str(len(output_list[i])), file=outputfile)   
        print(" i = "+str(i*100)+", time = "+str(i*100.0)+" s", file=outputfile)
        for line in output_list[i]:
            length=0
            length=len(line)
            sentence=" "
            for j in range(length):
                if j<length-1:
                    sentence=sentence+str(line[j])+"\t"
                else:
                    sentence=sentence+str(line[j])
            print(str(sentence), file=outputfile)
    outputfile.close()

start=0.0
finish=0.0
start=time.perf_counter()
output(coordinates_collect_4("dump.lammps.dat", alphabet_collect_element("output.lammps.dat")), "LAMMPS_trajectory.xyz")
finish=time.perf_counter()
print(f'coordinates_collect_4 finished in {round(finish-start, 4)} seconds')

