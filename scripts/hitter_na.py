#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
if len(sys.argv[1:]) == 0:
    print(
"""
Build an abundance table from a list of clusters and read abundance files.
usage: python3 hitter_na.py hitfile_list.txt clusters.otu project
Miguel Romero github.com/romeromig
"""
)
    sys.exit()

args = sys.argv[1:]

handle = open(args[0],'r')
hitfiles = handle.read().split('\n')[:-1]
handle.close()

#group all mapped CDSs files in a single list
hits = []
for filename in hitfiles:
    handle = open(filename,'r')
    rawfile = handle.read().split('\n')[:-1]
    handle.close()
    for line in rawfile:
        hits.append(line)

# load the cluster tables        
handle = open(args[1],'r')
clusterfileraw = handle.read().split('\n')[:-1]
handle.close()

# load the names of all mapped sequences
mapdseqs = []
for line in hits:
    mapdseqs.append(line.split(' ')[-1])

# recover clusters with mapped sequences
print('Looking for clusters with mapped reads.')

clstrs = []                             
for line in clusterfileraw:             # for each cluster
    linesep = line.split('\t')[1:-1]    # recover the sequence names
    for seq in linesep:                 # for each sequence
        seq = seq[1:]                    
        if seq in mapdseqs:             # if it is mapped
            if line not in clstrs:      
                clstrs.append(line)     # add it to the sequence list
                break                   # go to the next cluster

# each CDS points to its abundance
hit_dict = {}
for line in hits:
    linesep = line.split(' ')[-2:]
    hit_dict[linesep[1]] = int(linesep[0])

# load the sample names
samids = []
for filename in hitfiles:
    samids.append(filename.split('.hits')[0])

print('Building table.')

clstr_dict = {}                                                            # each cluster points to its abundance in each sample
for line in clstrs:                                                        # for each cluster
    linesep = line.split('\t')[:-1]
    clstr_dict['C'+ linesep[0]] = {}                                       # make a dictionary
    for sample in samids:                                                  # for each sample
        clstr_dict['C'+ linesep[0]][sample] = 0                            # begin counting
    for seq in linesep[1:]:                                                # for each sample in the cluster
        seq = seq[1:]                                                      
        samseq = seq.split('_')[0]                                         # get the sample name
        for sample in samids:                                              # then for each sample
            if samseq == sample:                                           # 
                try:                                                       # add its abundance
                    clstr_dict['C'+ linesep[0]][sample] += hit_dict[seq]
                except:                                                    
                    pass

# write the table
outfile = open( args[2] +'_na.tsv','w')
outfile.write('cluster\t'+ '\t'.join(samids) +'\n')

for clstr in clstr_dict.keys():
    columns = []
    columns.append(clstr)
    for sample in samids:
        columns.append(str(clstr_dict[clstr][sample]))
    outfile.write('\t'.join(columns) +'\n')
outfile.close()
print('Done. Final table saved in '+ args[2] +'_na.tsv')
