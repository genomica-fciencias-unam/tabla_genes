#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
if len(sys.argv[1:]) == 0:
    print(
"""
Genera una tabla de genes a partir de clusters y mapeo de lecturas
usage: python3 hitter_na.py hitfile_list.txt clusters.otu proyecto
Miguel Romero github.com/romeromig
"""
)
    sys.exit()

import os

args = sys.argv[1:]

handle = open(args[0],'r')
hitfiles = handle.read().split('\n')[:-1]
handle.close()

#agrupar todos los transcritos mapeados en un solo archivo
os.system('cat '+ ' '.join(hitfiles) +' > '+ args[2] +'.tmp1')

print('Buscando clusters con secuencias mapeadas (toma timepo)')
#obtener solo los clusters con marcos de lectura mapeados !tarda
os.system("for id in `awk -F "+'" "'+" '{print $2}' "+ args[2] +".tmp1`; do grep -P $id'\\t' "+ args[1] +" >> "+ args[2] +".tmp2; done")
os.system("sort "+ args[2] +".tmp2 | uniq > "+ args[2] +".tmp3")

handle = open(args[2] +'.tmp1','r')
hits = handle.read().split('\n')[:-1]
handle.close()

hit_dict = {}
for line in hits:
    linesep = line.split(' ')[-2:]
    hit_dict[linesep[1]] = int(linesep[0])

handle = open(args[2]+".tmp3",'r')
clstrs = handle.read().split('\n')[:-1]
handle.close()

placeholder = []
for i in range(len(hitfiles)):
    placeholder.append(0)

#a partir de aqui es necesario que los nombres de las muestras y secuencias coincidan
samids = []
for filename in hitfiles:
    samids.append(filename.split('.hits')[0])

print('Creando tabla')

clstr_dict = {}
for line in clstrs:
    linesep = line.split('\t')[:-1]
    clstr_dict['C'+ linesep[0]] = {} #crea un diccionario por cluster
    for sample in samids:
        clstr_dict['C'+ linesep[0]][sample] = 0 #crea un subdiccionario por muestra
    for seq in linesep[1:]: #por cada una de las secuencias
        seq = seq[1:]
        samseq = seq.split('_')[0]
        for sample in samids: #revisa a que muestra pertenece
            if samseq == sample:
                try: #si tiene reads mapeados, suma la abundancia
                    clstr_dict['C'+ linesep[0]][sample] += hit_dict[seq]
                except: #de otra forma, pasa
                    pass

outfile = open( args[2] +'_na.tsv','w')
outfile.write('cluster\t'+ '\t'.join(samids) +'\n')

for clstr in clstr_dict.keys():
    columns = []
    columns.append(clstr)
    for sample in samids:
        columns.append(str(clstr_dict[clstr][sample]))
    outfile.write('\t'.join(columns) +'\n')
outfile.close()
print('Hecho. Tabla final guardada en '+ args[2] +'_na.tsv')
print('Se recomienda borrar los archivos temporales: rm '+ args[2] +'.tmp*')