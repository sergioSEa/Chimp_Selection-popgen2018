#Script used for get the Specie-Cluster format from the pop.info file.
#How to use: python change_info_format.py name_of_file_to_transform
#Watch out: Before using the files for the different species should be splitted. Use grep NAME_OF_SPECIE > file_specie.info

import sys
file_o = sys.argv[1]
with open (file_o,"r") as file:
	for line in file:
		line = line.rstrip()
		l = line.split()
		try: 
			output = l[1]+" "+l[1]	
			out = open(file_o+"_good","a+b")
			out.write(output+"\n")
		except:
			pass
