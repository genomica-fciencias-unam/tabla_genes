#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
if len(sys.argv[1:]) == 0:
    print(
"""
Generar tablas de secuencias anotadas con conteos basados en mapeo de lecturas
usage: python3 hitter.py annotation.txt hitfile.txt muestra
Miguel Romero github.com/romeromig
"""
)
    sys.exit()

infiles = sys.argv[1:]

print('Cargando archivos')

# cargar tabla de anotacion
handle = open(infiles[0])
annot = handle.read().split('\n')[:-1]
handle.close()

# cada secuencia punta a su anotacion
annot_dict = {} 
for line in annot:
    linesep = line.split('\t')
    annot_dict[linesep[0]] = linesep[-1]

# cargar tabla de mapeo
handle = open(infiles[1])
hits = handle.read().split('\n')[:-1]
handle.close()

# cada secuencia apunta a su abundancia
hit_dict = {}
for line in hits:
    linesep = line.split(' ')[-2:]
    hit_dict[linesep[1]] = int(linesep[0])

print('Buscando secuencias anotadas y mapeadas')
chidas = []
for seqid in list(annot_dict.keys()):
    if seqid in list(hit_dict.keys()):
        chidas.append(seqid)

print('Creando tabla')

# sumar la abundancia por cada anotacion
annotcount_dict = {} 
for seqid in chidas:
    annotcount_dict[annot_dict[seqid]] = 0

print('Sumando abundancias')

for seqid in chidas:
    annotcount_dict[annot_dict[seqid]] += hit_dict[seqid]

print('Ordenando')

rownames = list(annotcount_dict.keys())
rownames.sort()

print('Escribiendo tabla')

# guardar en un archivo de texto
outfile = open(infiles[2]+'.hout','w') 
for annotid in rownames:
    outfile.write(annotid +'\t'+ str(annotcount_dict[annotid]) +'\n')
outfile.close()

print('Hecho. Tabla final guardada en '+ infiles[2] +'.hout')
