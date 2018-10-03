#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
if len(sys.argv[1:]) == 0:
    print(
"""
Agrupar tablas generadas por hitter.py
usage: python3 hitter_table.py hout_names_list proyecto
Miguel Romero github.com/romeromig
"""
)
    sys.exit()

import os

argvs = sys.argv[1:]

handle = open(argvs[0],'r')
infiles = handle.read().split('\n')[:-1]
handle.close()

os.system('cat '+ ' '.join(infiles) +" | awk '{print $1}' | sort | uniq > "+ argvs[1] +'.tmp1')

placeholder = []
for i in range(len(infiles)):
    placeholder.append('0')

handle = open(argvs[1] +'.tmp1','r')
md5s = handle.read().split('\n')[:-1]
handle.close()

proyect_dict = {} #generar tabla final vacia
for md5 in md5s:
    proyect_dict[md5] = placeholder.copy()

print('Construyendo tabla')
index = 0
for infile in infiles:
    print('Agregando archivo '+ infile)
    handle = open(infile,'r')
    sample = handle.read().split('\n')[:-1]
    handle.close()
    for line in sample:
        linesep = line.split('\t')
        proyect_dict[linesep[0]][index] = linesep[1]
    index += 1

outfile = open(argvs[1]+'.tsv','w') #guardar en un archivo de texto
outfile.write('md5\t'+ '\t'.join(infiles) +'\n')
for md5 in proyect_dict.keys():
    outfile.write(md5 +'\t'+ '\t'.join(proyect_dict[md5]) +'\n')
outfile.close()
print('Hecho. Tabla final guardada en '+ argvs[1] +'.tsv')
print('Se recomienda borrar el archivo temporal: rm '+ argvs[1] +'.tmp1')
