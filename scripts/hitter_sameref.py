#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
if len(sys.argv[1:]) == 0:
    print(
"""
Genera una tabla de conteos a partir de mapeo de lecturas contra una misma referencia
usage: python3 hitter_sameref.py hitfile_list.txt proyecto
Miguel Romero github.com/romeromig
"""
)
    sys.exit()

args = sys.argv[1:]

handle = open(args[0],'r')
hitfiles = handle.read().split('\n')[:-1]
handle.close()

# cargar los archivos de mapeo y crear diccionarios
print('Cargando archivos de mapeo y buscando secuencias únicas.')
hit_dict = {}
uniq_hit = {}
for filename in hitfiles:
    handle = open(filename,'r')
    rawfile = handle.read().split('\n')[:-1]
    handle.close()
    hit_dict[filename] = {}
    for line in rawfile:
        linesep = line.split(' ')[-2:]
        hit_dict[filename][linesep[1]] = int(linesep[0])
        uniq_hit[linesep[1]] = 0

# carga los nombres de las muestras
samids = []
for filename in hitfiles:
    samids.append(filename.split('.hits')[0])


print('Creando tabla.')
clstr_dict = {}                                                    # cada cluster apunta a su abundancia por muestra
for line in uniq_hit.keys():                                              # por secuencia única
    clstr_dict[line] = {}                                          # crea un diccionario
    for sample in samids:                                          # por muestra
        clstr_dict[line][sample] = 0                               # crea un subdiccionario con abundancia cero
        try:                                                       # si tiene reads mapeados, suma la abundancia
            clstr_dict[line][sample] += hit_dict[sample +'.hits'][line]
        except:                                                    # de otra forma, pasa
            pass

# escribe la tabla
print('Escribiendo tabla.')
outfile = open( args[1] +'.tsv','w')
outfile.write('sequence_id\t'+ '\t'.join(samids) +'\n')

for clstr in clstr_dict.keys():
    columns = []
    columns.append(clstr)
    for sample in samids:
        columns.append(str(clstr_dict[clstr][sample]))
    outfile.write('\t'.join(columns) +'\n')

outfile.close()
print('Hecho. Tabla final guardada en '+ args[1] +'.tsv')