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

args = sys.argv[1:]

handle = open(args[0],'r')
hitfiles = handle.read().split('\n')[:-1]
handle.close()

#agrupar todos los cdss mapeados en una lista
hits = []
for filename in hitfiles:
    handle = open(filename,'r')
    rawfile = handle.read().split('\n')[:-1]
    handle.close()
    for line in rawfile:
        hits.append(line)

# carga la tabla de clusters        
handle = open(args[1],'r')
clusterfileraw = handle.read().split('\n')[:-1]
handle.close()

# carga los nombres de las secuencias mapeadas
mapdseqs = []
for line in hits:
    mapdseqs.append(line.split(' ')[-1])

# obtener solo los clusters con marcos de lectura mapeados !tarda  MODIFICADO
print('Buscando clusters con secuencias mapeadas (toma tiempo)')

clstrs = []                             # lista de clusters mapeados
for line in clusterfileraw:             # por cada cluster
    linesep = line.split('\t')[1:-1]    # recupera los nombres de las secuencias
    for seq in linesep:                 # por cada secuencia
        seq = seq[1:]                   # elimina el espacio vacio 
        if seq in mapdseqs:             # si esta en la lista de cdss mapeados
            if line not in clstrs:      # y si el cluster no esta en la lista de clusters mapeados
                clstrs.append(line)     # agregalo a la lista
                break                   # cambia de cluster

# cada cdss apunta a su abundancia de lecturas mapeadas
hit_dict = {}
for line in hits:
    linesep = line.split(' ')[-2:]
    hit_dict[linesep[1]] = int(linesep[0])

# a partir de aqui es necesario que los nombres de las muestras y secuencias coincidan
# carga los nombres de las muestras i.e. prefijos de los cdss
samids = []
for filename in hitfiles:
    samids.append(filename.split('.hits')[0])

print('Creando tabla')

clstr_dict = {}                                                            # cada cluster apunta a su abundancia por muestra
for line in clstrs:                                                        # por cluster
    linesep = line.split('\t')[:-1]
    clstr_dict['C'+ linesep[0]] = {}                                       # crea un diccionario
    for sample in samids:                                                  # por muestra
        clstr_dict['C'+ linesep[0]][sample] = 0                            # crea un subdiccionario con abundancia cero
    for seq in linesep[1:]:                                                # por cada secuencia en el cluster
        seq = seq[1:]                                                      # elimina el espacio vacio
        samseq = seq.split('_')[0]                                         # obten el nombre de la muestra
        for sample in samids:                                              # y por muestra
            if samseq == sample:                                           # cuando encuentres a cual pertenece
                try:                                                       # si tiene reads mapeados, suma la abundancia
                    clstr_dict['C'+ linesep[0]][sample] += hit_dict[seq]
                except:                                                    # de otra forma, pasa
                    pass

# escribe la tabla
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
