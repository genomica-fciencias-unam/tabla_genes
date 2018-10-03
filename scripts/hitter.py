#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
if len(sys.argv[1:]) == 0:
    print(
"""
Generar tablas de secuencias anotadas con conteos basados en mapeo de lecturas
usage: python3 hitter.py annotation.txt hitfile.txt faa90.otu muestra
Miguel Romero github.com/romeromig
"""
)
    sys.exit()

import os

infiles = sys.argv[1:]

#buscar solo las secuencias que se mapearon y que se anotaron en m5nr
print('Buscando secuencias mapeadas y anotadas')
os.system("awk '{print $1}' "+ infiles[0] +'> '+ infiles[3] +'.tmp1') #secuencias anotadas
os.system("awk -F ' ' '{print $2}' "+ infiles[1] +'> '+ infiles[3] +'.tmp2' ) #secuencias mapeadas
os.system("cat "+ infiles[3] +".tmp1 "+ infiles[3] +".tmp2 | sort | uniq -c | grep '2 ' | awk -F ' ' '{print $2}' > "+ infiles[3] +".tmp3")

#obtener listas de secuencias en clusters
print('Obteniendo lista de clusters')
os.system("awk '{print $3}' "+ infiles[2] +" | sort | uniq > "+ infiles[3] +".tmp.clstrseeds")
os.system("for id in `cat "+ infiles[3] +".tmp.clstrseeds`; do grep -P $id'\\t' "+ infiles[2] +" >> "+ infiles[3] +".tmp.clstrs; done")

#buscar solo clusters con secuencias consideradas
print('Buscando clusters a conservar')
os.system("for id in `cat "+ infiles[3] +".tmp3`; do grep -P $id'\\t' "+ infiles[3] +".tmp.clstrs >> "+ infiles[3] +".tmp4; done")
os.system("sed -i 's/ //g' "+ infiles[3] +".tmp4") #eliminar espacios para facilitar procesamiento

handle = open(infiles[3] +'.tmp3','r')
seqs = handle.read().split('\n')[:-1]
handle.close()

handle = open(infiles[3] +'.tmp4','r')
clstrs_raw = handle.read().split('\n')[:-1] 
handle.close()

clstrs = [] #lista de listas de secuencias en clusters
for cluster in clstrs_raw:
    clstrs.append(cluster.split('\t')[1:-1])

handle = open(infiles[0],'r')
annot = handle.read().split('\n')[:-1]
handle.close()

annot_dict = {} #cada secuencia apunta a su md5
for line in annot:
    linesep = line.split('\t')
    annot_dict[linesep[0]] = linesep[-1]

handle = open(infiles[1],'r')
hits = handle.read().split('\n')[:-1]
handle.close()

hit_dict = {} #cada secuencia apunta a su abundancia
for line in hits:
    linesep = line.split(' ')[-2:]
    hit_dict[linesep[1]] = int(linesep[0])

print('Sumando conteos por cluster')
#corregir conteos por secuencias presentes en clusters
clstr_dict = {} #cada secuencia considerada apunta a su cluster

for seqid in seqs:
    for cluster in clstrs:
        if seqid in cluster:
            clstr_dict[seqid] = cluster
            break

for seqid in clstr_dict.keys():
    value = 0
    for member in clstr_dict[seqid]:
        try:
            value += hit_dict[member]
        except:
            pass
    hit_dict[seqid] = value #corregir por la suma del cluster

print('Sumando conteos por md5')
md5_dict = {} #sumar la abundancia por cada md5
for seqid in seqs:
    md5_dict[annot_dict[seqid]] = 0

for seqid in seqs:
    md5_dict[annot_dict[seqid]] += hit_dict[seqid]

outfile = open(infiles[3]+'.hout','w') #guardar en un archivo de texto
for md5 in md5_dict.keys():
    outfile.write(md5 +'\t'+ str(md5_dict[md5]) +'\n')
outfile.close()
print('Hecho. Tabla final guardada en '+ infiles[3] +'.hout')
print('Se recomienda borrar los archivos temporales: rm '+ infiles[3] +'.tmp*')