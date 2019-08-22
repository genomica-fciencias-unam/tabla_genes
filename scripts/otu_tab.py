#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
if len(sys.argv[1:]) == 0:
    print(
"""
Get a tsv table from an otu table
usage: python3 otu_tab.py sample_prefixes.txt clusters.otu
Miguel Romero github.com/romeromig
"""
)
    sys.exit()

args = sys.argv[1:]

# carga nombre de las muestras
handle = open(args[0],'r')
muestras = handle.read().split('\n')[:-1]
handle.close()

# carga la tabla de clusters
handle = open(args[1],'r')
clstrs = handle.read().split('\n')[1:-1]
handle.close()

# es necesario que los nombres de las muestras y secuencias coincidan (usar header.fasta.numbers.2.pl)
print('Generating table')

clstr_dict = {}                                                            # cada cluster apunta a su abundancia
for line in clstrs:                                                        # por cluster
    linesep = line.split('\t')[:-1]
    clstr_dict[linesep[0]] = {}                                       # crea un diccionario
    for sample in muestras:                                                  # por muestra
        clstr_dict[linesep[0]][sample] = 0                            # crea un subdiccionario con abundancia cero
    for seq in linesep[1:]:                                                # por cada secuencia en el cluster
        seq = seq[1:]                                                      # elimina el espacio vacio
        samseq = seq.split('_')[0]                                         # obten el nombre de la muestra
        clstr_dict[linesep[0]][samseq] += 1                                # y suma uno

# escribe la tabla
outfile = open( args[1] +'.tsv','w')
outfile.write('cluster\t'+ '\t'.join(muestras) +'\n')

for clstr in clstr_dict.keys():
    columns = []
    columns.append(clstr)
    for sample in muestras:
        columns.append(str(clstr_dict[clstr][sample]))
    outfile.write('\t'.join(columns) +'\n')
outfile.close()
print('Done, table saved in '+ args[1] +'.tsv')
