# James Stimson, 17 Dec 2018 to 23 Jan 2019.

import os
import sys

DEP_CUTOFF = 10

REMOVE_HYPER = 1

data_path = "/Users/js6509/indel_prjeb/"
data_path = "D:/James/indel_files/"
vcf_path = "D:/James/indel_calls/"

stats_file_name = "read_stats_" + str(DEP_CUTOFF) + ".txt"
stats_file = open(data_path + stats_file_name, 'w')

indels_loc_file_name = "indel_loc_" + str(DEP_CUTOFF) + ".txt"
indels_loc_file = open(data_path + indels_loc_file_name, 'w')

############### MASKING
hyper_file_name = "nucmer_mask_000962_3.bed"
hyper_file = open(data_path + hyper_file_name, 'r')

pom = hyper_file.readlines()

hyper = {} # simple list of hyper sites CHECK memory inefficient but quick look-up
count = 0
for i in range(len(pom)):
    if i==0: 
        continue
    info = pom[i].split('\t')
    count += 1
    for j in range(int(info[1]), int(info[2].rstrip())+1):
        hyper[j] = 1
        
hyper_file.close()
print ('Read in ' + str(len(hyper)) + ' hyper-variable sites from ' + str(count) + ' lines')
############### END MASKING

# Get the list of things we are working with 
label_file_name = "namelist.txt"
label_file = open(data_path + label_file_name, 'r')
count = 0
err_str = ''
name_list = []

for line in label_file:
    temp = line.split('\t') 
    name_list.append(temp[0].rstrip())
    count += 1
            
label_file.close()  
print ('Read ' + str(count) + ' names')

#sys.exit()


# Tab separated: Column 2 has location, Column 6 has qual, Column 10 is colon separated with 2nd item depth 
class indelcall:
    loc = 0
    depth = 0
    frac = 0
    ref = 'X'
    alt = 'X'
    svtype = 'SHORT'

# Pulls out list of indel calls from vcf file
def listfromvcf(vcf):
    mylist = []
    count = 0
    for line in vcf:
        count += 1
        if line[0] == '#':
            continue
        cols = line.split('\t')
        loc = int(cols[1].rstrip())
        ref = str(cols[3].rstrip())
        alt = str(cols[4].rstrip())

        qual = str(cols[5].rstrip())
        filter_pass = str(cols[6].rstrip())

        if filter_pass != 'PASS':
            continue

        info = (cols[7].rstrip()).split(';')
        
        svtype = 'SHORT'
        infotemp = ((info[0].rstrip()).split('='))[0]
        if infotemp == 'SVTYPE':
            svtype = ((info[0].rstrip()).split('='))[1]
            if svtype == 'DUP':
                continue
            dep = DEP_CUTOFF
            dep_frac = 0
        else:
            dep = int(((info[0].rstrip()).split('='))[1])  
            all_dep = ((info[1].rstrip()).split('='))[1]
            dep_frac = float(dep)/float(all_dep)
            dep_frac = int(dep_frac*100.0)  # Express as percentage
        
        
        if dep >= DEP_CUTOFF:   # Doesnt quite work for my purposes
            item = indelcall()
            item.loc = loc
            item.depth = dep
            item.frac = dep_frac
            item.ref = ref
            item.alt = alt
            item.svtype = svtype
            mylist.append(item)
    
    #print count
    return mylist




# Go through the sample dict and subs back in the ref if read is not deep enough
# Unfortunately this is too simplistic, it creates some false positives
# Before setting back to ref, need some kind of consensus 
# Maybe only set back to ref if no-one calls it 
# Perhaps we can set this up before the following loop
# That is, go through all the variable sites and flag 'keep' if we have at least one deep read
is_deep_enough = {}  # Label -> 1 if depth>cutoff, 0 otherwise
depth_stats = {}
count = 0
class stat:
    HighCount = 0
    LowCount = 0
    HighList = []
    LowList = []
    
    label = ''
 
index_dict = {}

# Label each position with depth flag. Note: this only covers the indexed ones
for i in index_dict:
    is_deep_enough[i] = 0   # Initialise to zero
    depth_stats[i] = stat()


# for each label, create list of calls, add to main list

main_list = []

name_count=0
for label in sorted(name_list):
    # Check each vcf file for deep reads
    file0 = open(vcf_path + 'indels' + label + '.recode.vcf')
    list0 = listfromvcf(file0)  # list of indel calls 
    #break # for TEST
    name_count+=1
    local_list = []
    print(label)
    if name_count >= 500:    # number of results files here
        break
    count = 0 
    for call in list0:
        index_dict[call.loc] = call.loc
        if call.loc not in depth_stats:
            depth_stats[call.loc] = stat()
        
        if call.depth >= DEP_CUTOFF:
            is_deep_enough[call.loc] = 1
            
            depth_stats[call.loc].HighCount += 1
            depth_stats[call.loc].HighList.append(label)

            local_list.append(call.loc)
            
        else:
            depth_stats[call.loc].LowCount += 1
            depth_stats[call.loc].LowList.append(label)
            
        depth_stats[call.loc].label = label
        
    file0.close()
    main_list.append(local_list)

    # output the locations to file
    indels_loc_file.write(label + ',')
    for call in local_list:
        indels_loc_file.write(str(call) + ',')
    indels_loc_file.write('\n')

    
#print (list0[5].loc)
#print (list0[5].svtype)
#print (len(list0))
#print (len(index_dict))



count = 0
for i in is_deep_enough:
    count += is_deep_enough[i]

for i in sorted(index_dict):
    stats_file.write(str(i) + ',' + str(depth_stats[i].HighCount) + ',' + str(depth_stats[i].LowCount) + ',' + str(depth_stats[i].label)+ '\n')

stats_file.close()
indels_loc_file.close()

sys.exit()

main_matrix = []

for local_list in main_list:
    fixed_list = []
    for index in sorted(index_dict):
        if index in local_list:
            fixed_list.append(1)
        else:
            fixed_list.append(0)    
    main_matrix.append(fixed_list)

# now call out pair-wise differences

diff_matrix = [[0 for x in range(len(main_matrix))] for y in range(len(main_matrix))]


for i in range(len(main_matrix)):
    for j in range(len(main_matrix)):
        mysum = 0
        for k in range(len(main_matrix[i])):
            if main_matrix[i][k] != main_matrix[j][k]:
                mysum += 1
        diff_matrix[i][j] = mysum

print(len(diff_matrix))

# output this to file     

matrix_file_name = "indel_matrix_" + str(DEP_CUTOFF) + ".csv"
matrix_file = open(data_path + matrix_file_name, 'w')
matrix_file.write(',')  
for s1 in sorted(name_list): 
    matrix_file.write(s1+',')
matrix_file.write('\n')   

#Then count the differences and output
#order matters
c1 = 0
c2 = 0

for s1 in sorted(name_list):  #order by
    output_line = s1+','
    c1 += 1
    for s2 in sorted(name_list):
        c2 += 1
        diff = 'NA'
        if c2<=len(diff_matrix) and c1<=len(diff_matrix):
            diff = diff_matrix[c1-1][c2-1]
        output_line += str(diff) + ','
    output_line += '\n'
    matrix_file.write(output_line)
    c2 = 0

matrix_file.close()    
    
