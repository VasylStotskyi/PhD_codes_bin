#!/usr/bin/env python3
import os
from pickle import FALSE, TRUE
from re import I
import sys
import platform
import os.path

def file_search(file_end):
	main=os.getcwd()
	path=[]
	for dirpath, dirnames, filenames in os.walk(main):
		for filename in [f for f in filenames if f.endswith(file_end)]:
			path.append(str(os.path.join(dirpath, filename)))
	return path

def file_start_search(file_start):
	main=os.getcwd()
	path=[]
	for dirpath, dirnames, filenames in os.walk(main):
		for filename in [f for f in filenames if f.startswith(file_start)]:
			path.append(str(os.path.join(dirpath, filename)))
	return path

def get_energy(filename):
    scf_check_boolean=False
    energy=0.0
    energy_list=[]
    with open (filename,"r") as readl:
        for line in readl:
            elements=line.split()
            if len(elements)==3:
                if elements[0]=="SCF" and elements[1]=="WAVEFUNCTION" and elements[2]=="OPTIMIZATION":
                    scf_check_boolean=TRUE
                    #print(scf_check_boolean)
            if scf_check_boolean:
                if len(elements)==9 and elements[0]=="ENERGY|" and elements[2]=="FORCE_EVAL":
                    energy=float(elements[8])
                    energy_list.append(energy)
                    scf_check_boolean= FALSE

    readl.close()
    return energy_list

def get_mixed_energy(filename):
    energy=0.0
    energy_list=[]
    with open (filename,"r") as readl:
        for line in readl:
            elements=line.split()
            if len(elements)==9 and elements[0]=="ENERGY|" and elements[4]=="MIXED":
                energy=float(elements[8])
                energy_list.append(energy)
    readl.close()
    return energy_list

def gap_calculator(energy1, energy2):
    energy_list=[]
    for i in range(0,len(energy1),1):
        energy=(energy1[i]-energy2[i])
        energy_eV=(energy1[i]-energy2[i])*27.2114
        energy_list.append([i,energy, energy_eV])
        #print(energy)
    return energy_list

def output_gap(energy_l, filename):
    output_file=open(filename, "w")
    sentence=""
    for i in range(0,len(energy_l),1):
        linelen=len(energy_l[i])
        for j in range(0,linelen-1,1):
            if type (energy_l[i][j]) is float:
                sentence+=str(round(energy_l[i][j],8))+"\t"    
            else:
                sentence+=str(energy_l[i][j])+"\t"
        if type (energy_l[i][linelen-1]) is float:
                sentence+=str(round(energy_l[i][linelen-1],8))+"\t"    
        else:
                sentence+=str(energy_l[i][linelen-1])+"\t"
        print(sentence, file=output_file)
        sentence=""
    output_file.close()
    return 0

output_r1_list=file_search("-r-1.out")
output_r2_list=file_search("-r-2.out")
output_r1_list.extend(output_r2_list)
output_r1_list.sort()
step=False
list_length_i=int(len(output_r1_list))
if list_length_i%2 !=0:
    print("Error! Not equal amount of *-r-1.out and *-r-1.out files.")
list_length_i=list_length_i/2
for eachline in output_r1_list:
    #print(eachline)
    if step:
        if temp!=eachline[:-5] :
            print("Error! Files paths doesn't mathch.")
    temp=eachline[:-5]
    step=not step
if len(output_r1_list)==0:
	print("Error! There is no *-r-1.out files in subfolders.")
	input()
	sys.exit(-2)
if len(output_r2_list)==0:
	print("Error! There is no *-r-2.out files in subfolders.")
	input()
	sys.exit(-2)
skipn=200
avg_mix_en_list=[]
m=0
calc_type="geo"
for eachline in output_r1_list:
    #print(eachline)
    if step:
        energy_l1 = get_energy(temp)
        energy_l2 = get_energy(eachline)
        energy_gap_l = gap_calculator(energy_l1, energy_l2)
        path_elem = temp.split(os.sep)
        path_new=""
        for i in range(len(path_elem)-1):
            path_new+=str(os.sep)+str(path_elem[i])
        path_new+=str(os.sep)
        #for i in energy_gap_l:
        #    print(i)
        energy_gap_l.remove(energy_gap_l[0])
        output_gap(energy_gap_l, path_new+"gapGEO.dat")
        if calc_type=="md":
            sum_mix_en=0.0
            for j in range(skipn,len(energy_gap_l)-1,1):
                sum_mix_en+=energy_gap_l[i][(len(energy_gap_l[i])-1)]
            avg_mix_en=sum_mix_en/(len(energy_gap_l)-skipn)
            avg_mix_en_list.append([round((1.0/(list_length_i-1))*m,5),avg_mix_en])
            avg_mix_en=0.0  
        else:
            avg_mix_en_list.append([round((1.0/(list_length_i-1))*m,5),energy_gap_l[len(energy_gap_l)-1][len(energy_gap_l[len(energy_gap_l)-1])-1]])
        #print(path_new)        
        m=m+1
    temp=eachline
    step=not step

def booles_law(avg_mixed_energy_list):
    booles_integral=0.0
    for Nstep in range(0,len(avg_mixed_energy_list),1):
        if Nstep==0 or Nstep==len(avg_mixed_energy_list)-1:
            booles_integral+=avg_mixed_energy_list[Nstep][len(avg_mixed_energy_list[Nstep])-1]*7.0
        elif Nstep%4==0:
            booles_integral+=avg_mixed_energy_list[Nstep][len(avg_mixed_energy_list[Nstep])-1]*14.0
        elif Nstep%2==0:
            booles_integral+=avg_mixed_energy_list[Nstep][len(avg_mixed_energy_list[Nstep])-1]*12.0
        else:
            booles_integral+=avg_mixed_energy_list[Nstep][len(avg_mixed_energy_list[Nstep])-1]*32.0
    booles_integral=booles_integral*2*(1/(len(avg_mixed_energy_list)-1))/45
    return booles_integral
    
def simpsons_law(avg_mixed_energy_list):
    a=avg_mixed_energy_list[0][0]
    b=avg_mixed_energy_list[len(avg_mixed_energy_list)-1][0]
    simpsons_integral=0.0
    for Nstep in range(0,len(avg_mixed_energy_list),1):
        if Nstep==0 or Nstep==len(avg_mixed_energy_list)-1:
            simpsons_integral+=avg_mixed_energy_list[Nstep][len(avg_mixed_energy_list[Nstep])-1]
        elif Nstep%2==0:
            simpsons_integral+=avg_mixed_energy_list[Nstep][len(avg_mixed_energy_list[Nstep])-1]*2.0
        else:
            simpsons_integral+=avg_mixed_energy_list[Nstep][len(avg_mixed_energy_list[Nstep])-1]*4.0
    simpsons_integral=simpsons_integral*((b-a)/(len(avg_mixed_energy_list)-1))/3.0
    return simpsons_integral    

def simpsons_law_38(avg_mixed_energy_list):
    a=avg_mixed_energy_list[0][0]
    b=avg_mixed_energy_list[len(avg_mixed_energy_list)-1][0]
    simpsons_integral=0.0
    for Nstep in range(0,len(avg_mixed_energy_list),1):
        if Nstep==0 or Nstep==len(avg_mixed_energy_list)-1:
            simpsons_integral+=avg_mixed_energy_list[Nstep][len(avg_mixed_energy_list[Nstep])-1]
        elif Nstep%3==0:
            simpsons_integral+=avg_mixed_energy_list[Nstep][len(avg_mixed_energy_list[Nstep])-1]*2.0
        else:
            simpsons_integral+=avg_mixed_energy_list[Nstep][len(avg_mixed_energy_list[Nstep])-1]*3.0
    simpsons_integral=simpsons_integral*3*((b-a)/(len(avg_mixed_energy_list)-1))/8.0
    return simpsons_integral   
"""
output_tmp_list=file_start_search("tmp.")   
output_tmp_list.sort()
for i in range(0,len(output_tmp_list),1):
    energy_mixed=get_mixed_energy(output_tmp_list[i])
    sum_mix_en=0.0
    avg_mix_en=0.0
    for j in range(skipn,len(energy_mixed),1):
        sum_mix_en+=energy_mixed[i]
    avg_mix_en=sum_mix_en/(len(energy_mixed)-skipn)
    avg_mix_en_list.append([(i)/10.0,avg_mix_en])
    avg_mix_en=0.0  
"""
#booles_integral=booles_law(avg_mix_en_list)
booles_integral=0.0
simpsons_integral=0.0
if (len(avg_mix_en_list)-1)%4== 0:
    booles_integral=booles_law(avg_mix_en_list)
elif (len(avg_mix_en_list)-1)%3== 0:
    simpsons_integral=simpsons_law_38(avg_mix_en_list)
elif(len(avg_mix_en_list)-1)%2 == 0:
    simpsons_integral=simpsons_law(avg_mix_en_list)
elif(len(avg_mix_en_list)-1)%2-1 == 0:
    avg_mix_en_list_13=[]
    avg_mix_en_list_38=[]
    avg_mix_en_list_13=avg_mix_en_list.copy()
    for l in range(0,3,1):
        avg_mix_en_list_38.insert(0,avg_mix_en_list_13.pop())
    temp_en_store=avg_mix_en_list_13.pop()
    avg_mix_en_list_38.insert(0,temp_en_store)
    avg_mix_en_list_13.append(temp_en_store)
#    print(avg_mix_en_list_13)
#    print(avg_mix_en_list_38)
    simpsons_integral=simpsons_law(avg_mix_en_list_13)
    simpsons_integral+=simpsons_law_38(avg_mix_en_list_38)
else:
    print("Error! At least 3 points are needed.")
if booles_integral!=0.0:
    avg_mix_en_list.append(["Free energy difference is equal to:", round(booles_integral,8), "eV (Boole's rule)"])
elif simpsons_integral!=0.0:
    avg_mix_en_list.append(["Formation energy is equal to:", round(simpsons_integral,8), "eV (Simpson's rule)"])
else:
    avg_mix_en_list.append(["Free energy difference is equal to:", round(0.0,8), "eV (Check your gaps!)"])
output_gap(avg_mix_en_list, str(os.getcwd())+str(os.sep)+"avg_mixed_energy.dat")
print("Done")