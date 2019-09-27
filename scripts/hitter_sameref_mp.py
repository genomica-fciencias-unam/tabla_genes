#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
import itertools
from multiprocessing import Pool, Manager

def greet():
    print(
    """
    Genera una tabla de conteos a partir de mapeo de lecturas contra una misma referencia
    usage: python3 hitter_sameref.py hitfile_list.txt proyecto cores
    Miguel Romero github.com/romeromig
    """
    )
    sys.exit()

def get_hits(filename,hit_dict,uniq_hit):
    sub_dict_hit={}
    sub_dict_uniq={}
    handle = open(filename+'.hits','r')
    rawfile = handle.read().split('\n')[:-1]
    handle.close()
    for line in rawfile:
        linesep = line.split(' ')[-2:]
        sub_dict_hit[(filename,linesep[1])] = int(linesep[0])
        sub_dict_uniq[linesep[1]] = 0
    hit_dict.update(sub_dict_hit)
    uniq_hit.update(sub_dict_uniq)

def func_star_get_hits(a_b):
    """Convert `f([1,2])` to `f(1,2)` call."""
    return get_hits(*a_b)

def build_tab(line,clstr_dict,samids):   
    sub_dict_clstr = {}                                           # por secuencia crea un diccionario
    for sample in samids:                                          # por muestra
        sub_dict_clstr[(line,sample)] = 0                               # crea un subdiccionario con abundancia cero
        try:                                                       # si tiene reads mapeados, suma la abundancia
            sub_dict_clstr[(line,sample)] += hit_dict[(sample,line)]
        except:
            pass
    clstr_dict.update(sub_dict_clstr)

def func_star_build_tab(a_b):
    """Convert `f([1,2])` to `f(1,2)` call."""
    return build_tab(*a_b)

# programa
    
def main():

    if len(sys.argv[1:]) == 0:
        greet()
    
    args = sys.argv[1:]
    
    handle = open(args[0],'r')
    hitfiles = handle.read().split('\n')[:-1]
    handle.close()
    
    # cargar los archivos de mapeo y crear diccionarios
    print('Cargando archivos de mapeo y buscando secuencias únicas.')
    
    # carga los nombres de las muestras
    samids = []
    for filename in hitfiles:
        samids.append(filename.split('.hits')[0])
    
    manager = Manager()
    hit_dict = manager.dict()
    uniq_hit = manager.dict()
    
    pool = Pool(processes=int(args[2]))
    pool.map(func_star_get_hits, zip(samids, itertools.repeat(hit_dict), itertools.repeat(uniq_hit)))
    
    pool.close()
    pool.join()
    pool.terminate()
    
    hit_dict = hit_dict.copy()                                      # necesario para acelerar
    uniq_hit = uniq_hit.copy()
    
    print('Creando tabla')
    
    uniqs = uniq_hit.keys()                                             # lista de las secuencias unicas
    
    tab_dict = manager.dict()                                        # cada cluster apunta a su abundancia por muestra
    
    pool = Pool(processes=int(args[2]))
    pool.map(func_star_build_tab, zip(uniqs, itertools.repeat(tab_dict), itertools.repeat(samids)))
    
    pool.close()
    pool.join()
    pool.terminate()
    
    tab_dict = tab_dict.copy()
    
    # escribe la tabla
    print('Escribiendo tabla.')
    outfile = open( args[1] +'.tsv','w')
    outfile.write('sample_id\t'+ '\t'.join(samids) +'\n')
    
    for line in uniqs:
        columns = []
        columns.append(line)
        for sample in samids:
            columns.append(str(tab_dict[(line,sample)]))
        outfile.write('\t'.join(columns) +'\n')
    
    outfile.close()
    print('Hecho. Tabla final guardada en '+ args[1] +'_na.tsv')

if __name__ == '__main__':
    main()
