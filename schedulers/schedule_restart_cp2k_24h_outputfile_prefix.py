import os

def file_search(file_end):
	main=os.getcwd()
	path=[]
	temp_path=[]
	for dirpath, dirnames, filenames in os.walk(main):
		for filename in [f for f in filenames if f.startswith(file_end)]:
			temp_path.append([dirpath, filename])
	for i in temp_path:
		path.append(str(os.path.join(i[0], i[1])))
	temp_path=[]
	return path

def read_files_logic(filename):
    result=-2
    line_index=0
    with open(filename, "r") as read_file:
        temp_lines=read_file.readlines()
        lines=temp_lines[-400:]
        temp_lines=[]
        #line_index=int(len(lines)-500)
        result=read_files_worker(line_index, lines)
        print("TEST EXIT RESULT:"+str(result))
        if result==1:
            os.system("sbatch *.job")
    read_file.close()

def read_files_worker(start_index, read_list):
    restart=0
    end_counter=16
    read_times=False
    for i in range(start_index, len(read_list)):
        elements=read_list[i].split()
        if len(elements)==1 and elements[0]=="-"*79:
            end_counter=end_counter-1
        if end_counter==0 and len(elements)==8 and elements[1]=="T" and elements[3]=="M":
            read_times=True
        if end_counter==-1 and read_times==True and len(elements)==7 and elements[0]=="CP2K":
            if float(elements[6])>=86000.0:
                restart=1
                return 1
            else:
                return 0
    return -1

filename=file_search("tmp.")
print(str(filename[-1]))
read_files_logic(str(filename[-1]))