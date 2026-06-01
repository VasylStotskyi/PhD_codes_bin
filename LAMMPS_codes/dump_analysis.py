import time
import math
#import concurrent.futures
#import multiprocessing as mp
import matplotlib.pyplot as plt
#from memory_profiler import profile

md_water_list=[]
md_water_list2=[]

def coordinates_collect_2(atom_type, atom_type2, input_filename): #slightly modified version of coordinates_collect_2, not in use any more. coordinates_collect_2 is significantly faster than coordinates_collect and uses much less RAM.
    step=0
    counter=-1.0
    temp_md_water_list=[]
    temp_md_water_list2=[]
    with open (input_filename) as inp:
        xmin_avg=0.0
        xmax_avg=0.0
        ymin_avg=0.0
        ymax_avg=0.0
        zmin_avg=0.0
        zmax_avg=0.0
        xmin=0.0
        xmax=0.0
        ymin=0.0
        ymax=0.0
        zmin=0.0
        zmax=0.0
        x=0
        y=0
        z=0
        step_bool=False
        atom_bool=False
        cell_iteration=0
        for line in inp:
            element=line.split()
            if (step_bool==True):
                step=int(element[0])
                counter=counter+1.0
                step_bool=False
            if ((len(element) > 1) and (element[1]=="NUMBER")):
                step_bool=True
            if (cell_iteration==3):
                cell_iteration=0
                zmin=float(element[0])
                zmax=float(element[1])
                z=zmax-zmin
                zmin_avg=zmin_avg+zmin
                zmax_avg=zmax_avg+zmax
            if (cell_iteration==2):
                cell_iteration=cell_iteration+1
                ymin=float(element[0])
                ymax=float(element[1])
                y=ymax-ymin
                ymin_avg=ymin_avg+ymin
                ymax_avg=ymax_avg+ymax
            if (cell_iteration==1):
                cell_iteration=cell_iteration+1
                xmin=float(element[0])
                xmax=float(element[1])
                x=xmax-xmin
                xmin_avg=xmin_avg+xmin
                xmax_avg=xmax_avg+xmax
            if len(element) > 1 and element[1]=="BOX" and counter>0.0:
                cell_iteration=cell_iteration+1
            if (atom_bool==True) and (step>0):
                step=step-1
                # if step==1:
                #     print(line)
                    # time.sleep(5)
                if (element[1]==atom_type):
                    temp_md_water_list.append([float(element[2]), float(element[3]), float(element[4])])
                if (element[1]==atom_type):
                    temp_md_water_list2.append([float(element[2]), float(element[3]), float(element[4])])
                if step==0:
                    atom_bool=False
                    # print(line)
                    # time.sleep(5)
                    # temp_md_water_list=normalized_coordinates(xmin, xmax, x, ymin, ymax, y, zmin, zmax, z, temp_md_water_list)
                    # temp_md_water_list2=normalized_coordinates(xmin, xmax, x, ymin, ymax, y, zmin, zmax, z, temp_md_water_list2)
                    # bond_list.extend(bond_eval(temp_md_water_list, temp_md_water_list2, xmax-xmin, ymax-ymin, zmax-zmin))
                    md_water_list.extend(temp_md_water_list)
                    md_water_list2.extend(temp_md_water_list2)
                    temp_md_water_list=[]
                    temp_md_water_list2=[]
                    # time.sleep(10)
            if len(element) > 1 and element[1]=="ATOMS" and counter>0.0:
                atom_bool=True
        # bond_list.insert(0,[counter, xmin_avg/counter, xmax_avg/counter, ymin_avg/counter, ymax_avg/counter, zmin_avg/counter, zmax_avg/counter])
        md_water_list.insert(0,[counter, xmin_avg/counter, xmax_avg/counter, ymin_avg/counter, ymax_avg/counter, zmin_avg/counter, zmax_avg/counter])
    inp.close()
    return md_water_list

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

def center_of_mass_calculate(atoms_list,  bond_list):   #calculate centers of mass of water molecules within the given bonding list
    center_of_mass_mx=0.0
    center_of_mass_my=0.0
    center_of_mass_mz=0.0
    center_of_mass_list=[]
    mass_O=15.9994
    mass_H=1.0008
    for i in bond_list:
        center_of_mass_mx=(atoms_list[i[0]][0]*mass_O+atoms_list[i[1]][0]*mass_H+atoms_list[i[2]][0]*mass_H)/(mass_O+mass_H*2)
        center_of_mass_my=(atoms_list[i[0]][1]*mass_O+atoms_list[i[1]][1]*mass_H+atoms_list[i[2]][1]*mass_H)/(mass_O+mass_H*2)
        center_of_mass_mz=(atoms_list[i[0]][2]*mass_O+atoms_list[i[1]][2]*mass_H+atoms_list[i[2]][2]*mass_H)/(mass_O+mass_H*2)
        center_of_mass_list.append([center_of_mass_mx, center_of_mass_my, center_of_mass_mz])
    return center_of_mass_list

def coordinates_collect_3(atom_type, input_filename):   #2 times faster than coordinates_collect_2 because of optimized reading code and smaller amount of if checks
    step=0
    counter=-1.0
    with open (input_filename, "r") as inp:
        xmin_avg=0.0
        xmax_avg=0.0
        ymin_avg=0.0
        ymax_avg=0.0
        zmin_avg=0.0
        zmax_avg=0.0
        step_bool=False
        atom_bool=False
        cell_iteration=False
        for line in inp:
            element=line.split()
            if (step_bool==True):
                step=int(element[0])
                counter=counter+1.0
                line=next(inp)
                element=line.split()
                step_bool=False
            if ((len(element) == 4) and (element[1]=="NUMBER")):
                step_bool=True
                continue
            if counter<=0:
                continue
            if (cell_iteration==True):
                xmin_avg=xmin_avg+float(element[0])
                xmax_avg=xmax_avg+float(element[1])
                line=next(inp)
                element=line.split()
                ymin_avg=ymin_avg+float(element[0])
                ymax_avg=ymax_avg+float(element[1])
                line=next(inp)
                element=line.split()
                zmin_avg=zmin_avg+float(element[0])
                zmax_avg=zmax_avg+float(element[1])
                line=next(inp)
                element=line.split()
                cell_iteration=False
            if len(element) == 6 and element[1]=="BOX":
                cell_iteration=True
                continue
            if (atom_bool==True) and (step>0):
                for i in range(1, step):
                    if (element[1]==atom_type):
                        md_water_list.append([float(element[2]), float(element[3]), float(element[4])])
                    line=next(inp)
                    element=line.split()
                if (element[1]==atom_type):
                        md_water_list.append([float(element[2]), float(element[3]), float(element[4])])
                atom_bool=False
                continue
            if len(element) == 7 and element[1]=="ATOMS" and counter>0.0:
                atom_bool=True
                continue
        md_water_list.insert(0,[counter, xmin_avg/counter, xmax_avg/counter, ymin_avg/counter, ymax_avg/counter, zmin_avg/counter, zmax_avg/counter])
    inp.close()
    return md_water_list

def coordinates_collect_4(atom_type, atom_type2, input_filename, bond_filename):   #1.5 times slower than coordinates_collect_3, because also calculates ceters of mass of water molecule
    step=0
    counter=-1.0
    temp_md_water_list=[]
    bond_list=[]
    bond_list=read_bonds(bond_filename)
    with open (input_filename, "r") as inp:
        xmin_avg=0.0
        xmax_avg=0.0
        ymin_avg=0.0
        ymax_avg=0.0
        zmin_avg=0.0
        zmax_avg=0.0
        xmin=0.0
        xmax=0.0
        ymin=0.0
        ymax=0.0
        zmin=0.0
        zmax=0.0
        x=0
        y=0
        z=0
        step_bool=False
        atom_bool=False
        cell_iteration=False
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
            if (cell_iteration==True):   #reads cell dimensions
                xmin=float(element[0])
                xmax=float(element[1])
                x=xmax-xmin
                xmin_avg=xmin_avg+xmin
                xmax_avg=xmax_avg+xmax
                line=next(inp)
                element=line.split()
                ymin=float(element[0])
                ymax=float(element[1])
                y=ymax-ymin
                ymin_avg=ymin_avg+ymin
                ymax_avg=ymax_avg+ymax
                line=next(inp)
                element=line.split()
                zmin=float(element[0])
                zmax=float(element[1])
                z=zmax-zmin
                zmin_avg=zmin_avg+zmin
                zmax_avg=zmax_avg+zmax
                line=next(inp)
                element=line.split()
                cell_iteration=False
            if len(element) == 6 and element[1]=="BOX": #activates cell dimensions reading
                cell_iteration=True
                continue
            if (atom_bool==True) and (step>0):   #reads O and H coordinates
                for i in range(1, step):
                    if (element[1]==atom_type) or (element[1]==atom_type2):
                        temp_md_water_list.append([float(element[2]), float(element[3]), float(element[4])])
                    line=next(inp)
                    element=line.split()
                if (element[1]==atom_type) or (element[1]==atom_type2):
                    temp_md_water_list.append([float(element[2]), float(element[3]), float(element[4])])
                atom_bool=False
                md_water_list.extend(center_of_mass_calculate(temp_md_water_list, bond_list))  #calculates centers of mass of water molecules and saves them in list; potentially can be parallelized
                temp_md_water_list=[]
                continue
            if len(element) == 7 and element[1]=="ATOMS":   #activates atomic coordinates reading
                atom_bool=True
                continue
        output(md_water_list, "water.dat")
        md_water_list.insert(0,[counter, xmin_avg/counter, xmax_avg/counter, ymin_avg/counter, ymax_avg/counter, zmin_avg/counter, zmax_avg/counter])   #adds averaged cell dimensions to list on first position. It's required for normalization of atomic coordinates.
    inp.close()
    return md_water_list

#with concurrent.futures.ProcessPoolExecutor() as executor:
#    results=[executor.submit( parallel_collect, atom_type, xmin_avg, xmax_avg, ymin_avg, ymax_avg, zmin_avg, zmax_avg, lines) for lines in range (24)]
#for f in concurrent.futures.as_completed(results):

def coordinates_collect(atom_type, input_filename):   #highly inefficient from the RAM usage point of view
    step=0
    counter=-1.0
    with open (input_filename) as inp:
        lines=inp.readlines()   #remove
        xmin_avg=0.0
        xmax_avg=0.0
        ymin_avg=0.0
        ymax_avg=0.0
        zmin_avg=0.0
        zmax_avg=0.0
        length_lines=len(lines)
        for i in range(0, length_lines):
            line=lines[i]
            element=line.split()
            length=len(element)
            if ((length > 1) and (element[1]=="NUMBER")):
                i=i+1
                line=lines[i]
                element=line.split()
                step=int(element[0])
                counter=counter+1.0
            if length > 1 and element[1]=="BOX" and counter>0.0:
                i=i+1
                line=lines[i]
                element=line.split()
                xmin_avg=xmin_avg+float(element[0])
                xmax_avg=xmax_avg+float(element[1])
                i=i+1
                line=lines[i]
                element=line.split()
                ymin_avg=ymin_avg+float(element[0])
                ymax_avg=ymax_avg+float(element[1])
                i=i+1
                line=lines[i]
                element=line.split()
                zmin_avg=zmin_avg+float(element[0])
                zmax_avg=zmax_avg+float(element[1])
            if length > 1 and element[1]=="ATOMS" and counter>0.0:
                for j in range (step, 0, -1):
                    i=i+1
                    line=lines[i]
                    element=line.split()
                    if element[1]==atom_type:
                        md_water_list.append([float(element[2]), float(element[3]), float(element[4])])
        md_water_list.insert(0,[counter, xmin_avg/counter, xmax_avg/counter, ymin_avg/counter, ymax_avg/counter, zmin_avg/counter, zmax_avg/counter])
    inp.close()
    return md_water_list

def binary_search_closest_index(vmin, vmax, atomic_shift_vector, search_list, vsearch):   #searches the closest right position in list. Uses recursive binary search algorithm. Thererefore, coordinates have to be sorted.
    vmid=0
    vmid=int((vmax+vmin)/2)
    if  vmid==vmin:
        return vmid
    elif (search_list[vmid][atomic_shift_vector]>vsearch):
        vmid=binary_search_closest_index(vmin,vmid, atomic_shift_vector, search_list, vsearch)
    elif (search_list[vmid][atomic_shift_vector]<=vsearch):
        vmid=binary_search_closest_index(vmid,vmax, atomic_shift_vector, search_list, vsearch)
    else:
        print("END")
        return -3
    return vmid

def closest_index(search_list, atomic_shift_vector, vsearch):   #controls binary search. Activates a few fast responces if searched values are smaller or larger than the one existing in the list 
    vmin=0
    vmax=len(search_list)
    if (search_list[0][atomic_shift_vector]>vsearch):
        return 0
    elif(search_list[vmax-1][atomic_shift_vector]<vsearch):
        return vmax
    else:
        result=binary_search_closest_index(vmin, vmax, atomic_shift_vector, search_list, vsearch)
        return result+1
    
def read_temp():   #function for quick read of already sorted and normalized water coordinates for density calculation check
    filename="temp_water.dat"
    temp_water_list=[]
    with open (filename, "r") as readl:
        for line in readl:
            element=line.split()
            temp_water_list.append([float(element[0]),float(element[1]),float(element[2])])
    readl.close()
    return temp_water_list

def read_bonds(filename):
    read_bonds=False
    temp_bond_list=[]
    bond_list=[]
    with open(filename, "r") as inp:
        for line in inp:
            element=line.split()
            if len(element)==1 and element[0]=="Bonds":
                read_bonds=True
                line=next(inp)
                continue
            elif read_bonds==True:
                if len(element)==0:
                    read_bonds=False
                    break
                elif element[1]=="1":
                    temp_bond_list.append([int(element[2]), int(element[3])])
                else:
                    continue
            else:
                continue
    inp.close()
    temp_bond_list.sort()
    O_minimum=1000000
    O_minimum=temp_bond_list[0][0]
    for i in range(1, (len(temp_bond_list)), 2):
        if temp_bond_list[i-1][0]==temp_bond_list[i][0]:
            bond_list.append([temp_bond_list[i][0]-O_minimum, temp_bond_list[i-1][1]-O_minimum, temp_bond_list[i][1]-O_minimum])
            if i+2<(len(temp_bond_list)) and (temp_bond_list[i+2][0]-O_minimum) != (i/2.0+0.5)*3 :
                O_minimum+=(temp_bond_list[i+2][0]-temp_bond_list[i][0]-int((i/2.0-0.5)*3))
        else:
            print("Error!")
    return bond_list

def quick_normalizator(atom, xmin, xmax, x, ymin, ymax, y, zmin, zmax, z):
    if atom[0]<xmin:
        while(atom[0]<xmin):
            atom[0]=atom[0]+x
    if atom[0]>=xmax:
        while(atom[0]>=xmax):
            atom[0]=atom[0]-x
    if atom[1]<ymin:
        while(atom[1]<ymin):
            atom[1]=atom[1]+y
    if atom[1]>=ymax:
        while(atom[1]>=ymax):
            atom[1]=atom[1]-y
    if atom[2]<zmin:
        while(atom[2]<zmin):
            atom[2]=atom[2]+z
    if atom[2]>=zmax:
        while(atom[2]>=zmax):
            atom[2]=atom[2]-z
    return atom

def quick_shifter(atom, xmin, xmax, x, ymin, ymax, y, zmin, zmax, z, shift, nround):
    atom=quick_normalizator(atom, xmin, xmax, x, ymin, ymax, y, zmin, zmax, z)
    if shift==0:
        shift_value_min=xmin
        shift_value_max=xmax
    if shift==1:
        shift_value_min=ymin
        shift_value_max=ymax
    if shift==2:
        shift_value_min=zmin
        shift_value_max=zmax
    if atom[shift]<nround:
        atom[shift]=atom[shift]-shift_value_min
    else:
        atom[shift]=atom[shift]-shift_value_max
    return atom

def normalized_coordinates2(xmin, xmax, x, ymin, ymax, y, zmin, zmax, z, atomic_shift_vector, atomic_coordinates_list):   #normalize and sort atomic coordinates. Shall be parallelized
    for atom in atomic_coordinates_list:
        atom=quick_shifter(atom, xmin, xmax, x, ymin, ymax, y, zmin, zmax, z, atomic_shift_vector, 5)
    atomic_coordinates_list.sort(key= lambda x: x[atomic_shift_vector])   #sorts atomic coordinates
    minimum_value=0
    minimum_value=atomic_coordinates_list[0][atomic_shift_vector]
    for atom in atomic_coordinates_list:    #then it shifts all atomic coordinates to nround value (middle of the cell rounded to 10). New coordinates starts from nround.
        atom[atomic_shift_vector]=atom[atomic_shift_vector]-minimum_value#+nround
    #output(atomic_coordinates_list, "temp_water.dat")   #prints out water coordinates. It's required only when density calculation function is checked/modified
    return atomic_coordinates_list

def normalized_coordinates(xmin, xmax, x, ymin, ymax, y, zmin, zmax, z, atomic_coordinates_list):   #normalize and sort atomic coordinates. Shall be parallelized
    for atom in atomic_coordinates_list:
        if atom[0]<xmin:
            while(atom[0]<xmin):
                atom[0]=atom[0]+x
        if atom[0]>=xmax:
            while(atom[0]>=xmax):
                atom[0]=atom[0]-x
        if atom[1]<ymin:
            while(atom[1]<ymin):
                atom[1]=atom[1]+y
        if atom[1]>=ymax:
            while(atom[1]>=ymax):
                atom[1]=atom[1]-y
        if atom[2]<zmin:
            while(atom[2]<zmin):
                atom[2]=atom[2]+z
        if atom[2]>=zmax:
            while(atom[2]>=zmax):
                atom[2]=atom[2]-z
        if atom[0]<5:
            atom[0]=atom[0]-xmin
        else:
            atom[0]=atom[0]-xmax
    atomic_coordinates_list.sort()   #sorts atomic coordinates
    x_minimum=atomic_coordinates_list[0][0]
    for atom in atomic_coordinates_list:    #then it shifts all atomic coordinates to nround value (middle of the cell rounded to 10). New coordinates starts from nround.
        atom[0]=atom[0]-x_minimum
    #output(atomic_coordinates_list, "temp_water.dat")   #prints out water coordinates. It's required only when density calculation function is checked/modified
    return atomic_coordinates_list

def output(output_list, output_filename):   #prints any list with tabs between elements of list
    outputfile=open(output_filename, "w")
    for line in output_list:
        length=0
        length=len(line)
        sentence=""
        for i in range(length):
            if i<length-1:
                sentence=sentence+str(line[i])+"\t"
            else:
                sentence=sentence+str(line[i])
        print(str(sentence), file=outputfile)
    outputfile.close()

def water_counter(water_list,d_x,x_min,x_max):   #first version of water counter. Extremely slow and inefficient. Has a systematic error: accounts 20.0 (first position as the one, which belongs to next slice).
    count_list=[]
    x_max=round(x_max-x_min+10.0,4)
    x_min=0.0
    count=0.0
    while(x_min < x_max):
        for water in water_list:
            if (water[0]>=x_min) and (water[0]<x_min+d_x):
                count=count+1.0
        count_list.append([x_min, round(float(count),4)])
        x_min=round(x_min+d_x,4)
        count=0.0
    return count_list

def water_counter2(water_list,d_x,x_min,x_max, atomic_shift_vector):   #improved version of water_counter. Speeds up code a lot thanks to recursive binary search
    count_list=[]
    x_max=round(x_max-x_min,4)
    x_min=round(0,4)
    count=0.0
    temp_count1=0
    temp_count2=0
    while(x_min < x_max):
        temp_count2=closest_index(water_list, atomic_shift_vector, x_min)
        count=temp_count2-temp_count1
        count_list.append([x_min, round(float(count),4)])
        count=0.0
        temp_count1=temp_count2
        x_min=round(x_min+d_x,4)
    return count_list

def density_calculation(md_list, shifter="x"):   #function which calculates water density and plot water density profile
    dx=0.2
    y=0.0
    z=0.0
    counter=0.0
    density_list=[]
    for md in md_list:
        length=len(md)
        if length==7:
            counter=md[0]
            xmin=md[1]
            xmax=md[2]
            ymin=md[3]
            ymax=md[4]
            zmin=md[5]
            zmax=md[6]
            x=xmax-xmin
            y=ymax-ymin
            z=zmax-zmin
            del md_list[0]
            break
    if shifter=="x":
        #nround=round(int(x/2)/10)*10
        atomic_shift_vector=0
        dmin=xmin
        dmax=xmax
        volume=dx*y*z
    elif shifter=="y":
        #nround=round(int(y/2)/10)*10
        atomic_shift_vector=1
        dmin=ymin
        dmax=ymax
        volume=dx*x*z
    else:
        #nround=round(int(z/2)/10)*10
        atomic_shift_vector=2
        dmin=zmin
        dmax=zmax
        volume=dx*y*x
    #nround=round(int((xmax-xmin)/2)/10)*10-10
    start_l=0.0
    finish_l=0.0
    start_l=time.perf_counter()
    md_list=normalized_coordinates2(xmin, xmax, x, ymin, ymax, y, zmin, zmax, z, atomic_shift_vector, md_list)
    finish_l=time.perf_counter()
    print(f'normalizator2 finished in {round(finish_l-start_l, 4)} seconds')
    start_l=0.0
    finish_l=0.0
    start_l=time.perf_counter()
    density_list=water_counter2(md_list,dx,dmin,dmax, atomic_shift_vector)
    finish_l=time.perf_counter()
    print(f'counter2 finished in {round(finish_l-start_l, 4)} seconds')
    output(density_list, "temp_counter_2.dat")
    md_list.clear()
    start_l=0.0
    finish_l=0.0
    start_l=time.perf_counter()
    for number in density_list:
        multiplier=(18/(6.02214076*1E23))/(counter*volume*1E-24)
        number[1]=number[1]*multiplier
    finish_l=time.perf_counter()
    print(f'multipl finished in {round(finish_l-start_l, 4)} seconds')
    x_list=[]
    y_list=[]
    step=1
    density=0
    temp_step=0
    temp_density=0
    xv=0
    x_middle=0
    x_middle=(density_list[len(density_list)-1][0]-density_list[0][0])/4
    for number in density_list:
        x_list.append(number[0])
        if (number[1]<0.95):
            if (temp_step*dx)>=7.5:
                step=temp_step
                density=temp_density
            temp_step=0
            temp_density=0
        elif (number[1]>1.05):
            if (temp_step*dx)>=7.5:
                step=temp_step
                density=temp_density
            temp_step=0
            temp_density=0
        else:
            if number[0]<x_middle:
                if number[1]>=0.94 and temp_step==0:
                    xv=number[0]
            temp_step=temp_step+1
            temp_density=temp_density+number[1]
        y_list.append(number[1])
    density=density+temp_density
    step=step+temp_step
    density=round(density/(step),4)
    x_dimens=round((step)*dx,4)
    fig, fx = plt.subplots()
    fx.bar(x_list,y_list, width=dx, align='center', color='black')
    fx.axhline(y=1.0, color='red')
    fx.axvline(x=xv, linewidth=1.0, color='red')
    fx.axvline(x=x_dimens+xv-0.2, linewidth=1.0, color='red')
    #x=round(int((xmax-xmin)/2)/10)*10
    print(density)
    print(x_dimens)
    fx.set_xlabel(fr'x dimension, $\AA$')
    fx.set_ylabel(r'Density, $\rho$')
    fx.set_title(f'Water density distribution profile:\n average density {density} $g/cm^3$, water layer width {x_dimens} $cm$')
    fig.savefig("density_plot.svg", dpi=1200.0)
    return density_list



start=0.0
finish=0.0
start=time.perf_counter()
coordinates_collect_4("1", "3","dump.lammps.dat", "output.lammps.dat")
finish=time.perf_counter()
print(f'coordinates_collect_4 finished in {round(finish-start, 4)} seconds')

# start=0.0
# finish=0.0
# start=time.perf_counter()
# coordinates_collect_3("1","dump.lammps.dat")
# finish=time.perf_counter()
# print(f'coordinates_collect_3 finished in {round(finish-start, 4)} seconds')

# start=0.0
# finish=0.0
# start=time.perf_counter()
# coordinates_collect_2("1","dump.lammps.dat")
# finish=time.perf_counter()
# print(f'coordinates_collect_2 finished in {round(finish-start, 4)} seconds')

start=0.0
finish=0.0
start=time.perf_counter()
output(density_calculation(md_water_list,"x"), "density_list.dat")
finish=time.perf_counter()
print(f'density calculator finished in {round(finish-start, 4)} seconds')
