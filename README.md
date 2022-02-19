Instrucciones para generar una tabla de genes anotados y no anotados, con conteos basados en mapeo de lecturas.
Los scripts mencionados se encuentran disponibles en este repositorio.

## Obtención de tablas de abundancia de genes

Archivos iniciales:

*acacia.faa* contiene secuencias de aminoácidos renombradas. 
(Los nombres de las secuencias deben tener el nombre de la muestra como prefijo y el número de secuencia como sufijo)

*acacia.fna* contiene secuencias de nucleótidos de los marcos de lectura renombradas.
(Los nombres de las secuencias deben tener el nombre de la muestra como prefijo y el número de secuencia como sufijo y las secuencias deben corresponder con las del archivo de aminoácidos)

*acacia-rn.fasta* contiene secuencias de lecturas de calidad tanto pareadas como sin parear.

### Abundancias
Generar librería de los marcos de lectura para mapear con bowtie2.

`bowtie2-build acacia.fna acacia`

Alinear lecturas a marcos de lectura para asignar abundancias.

`bowtie2 -f -x acacia -U acacia-rn.fasta -S acacia.sam --quiet -p 20 --very-sensitive`

Obtener únicamente alineamientos de calidad y recuperar las secuencias de referencia mapeadas con su frecuencia.

`grep -v '^@' acacia.sam | awk '{if($5 == 3 || $5 == 8 || $5 == 23 || $5 == 24 || $5 == 40 || $5 == 42) print $3}' | sort | uniq -c > acacia.hits`

### Anotar

Partir el archivo de secuencias para paralelizar la anotación

`partefasta 10000 acacia90.faa`

Correr la anotación contra la base de datos m5nr en paralelo por cada archivo de 10,000 secuencias

`diamond blastp -d m5nr -q acacia90.faa.1.fas -f 6 -e 1e-10 -k 10 -p 1 --quiet -o acacia90.faa.1.fas.bout`

Unir todas las salidas del archivo *acacia.faa*

`cat acacia90.faa.*.bout > acacia-m5nr.bout`

Oneliner para ordenar la salida de blast por valor de bitscore, remover duplicados con el mismo valor de bitscore y guardar en el archivo *acacia_best_uniq*

`
cat acacia-m5nr.bout | perl -pe ' $name_col=0; $score_col=11; while(<>) { s/\r?\n//; @F=split /\t/, $_; ($n, $s) = @F[$name_col, $score_col]; if (! exists($max{$n})) { push @names, $n }; if (! exists($max
{$n}) || $s > $max{$n}) { $max{$n} = $s; $best{$n} = () }; if ($s == $max{$n}) { $best{$n} .= "$_\n" }; } for $n (@names) { print $best{$n} } ' >best;  perl -e ' $column=0; $unique=0; while(<>) { s/\r?\n//; @F=s
plit /\t/, $_; if (! ($save{$F[$column]}++)) { print "$_\n"; $unique++ } } ' best > acacia_best_uniq; rm best
`

Simplificar la salida

`awk '{print $1"\t"$3"\t"$11"\t"$12"\t"$2}' acacia_best_uniq > acacia_best.simple.tsv`

### Construir tablas

#### Importante: Para utilizar los scripts que generan las tablas de abundancias (hitter.py, hitter_table.py y hitter_na.py) es imprescindible que los nombres de las secuencias tengan como prefijo el nombre de la muestra (y como sufijo el número de la secuencia) y que los archivos de secuencias mapeadas tengan el nombre de la muestra. Los nombres de las muestras deben ser consistentes  y no deben contener guiones bajos.

Por ejemplo, el archivo *acacia.hits* debe verse así:

```
      1 acacia_10
      2 acacia_100
     15 acacia_10000
     72 acacia_100000
    116 acacia_100001
    121 acacia_100002
     61 acacia_100003
     18 acacia_100004
    147 acacia_100005
      4 acacia_100006 
```


Generar una tabla que contenga cada identificador del m5nr y su abundancia según las lecturas mapeadas a los marcos de lectura. Se usa el script hitter.py y se obtiene el archivo *acacia.hout*

`python3 hitter.py acacia_best.simple.tsv acacia.hits acacia`

Se usó el mismo método para las muestras cca y ccac. Crear una lista con las muestras a juntar en una tabla

`ls *.hout > lista`

Unir las muestras en una sola tabla con el script hitter_table.py se obtiene el archivo *nac.tsv*

`python3 hitter_table.py lista nac`

### Agregar secuencias no anotadas agrupadas por identidad

Obtener lista de secuencias anotadas

`awk '{print $1}' acacia_best.simple.tsv > acacia_anotados.txt`

Obtener lista de todas las secuencias

`grep '>' acacia.faa | sed 's/>//g' > acacia_todos.txt`

Recuperar los nombres de las secuencias no anotadas

`cat acacia_anotados.txt acacia_todos.txt | sort | uniq -c | grep '1 ' | awk '{print $2}' > acacia_na.txt`

Extraer las secuencias no anotadas

`seqtk subseq acacia.faa acacia_na.txt > acacia_na.faa`

Unir las secuencias no anotadas de todas las muestras

`cat *_na.faa > todos_na.faa`

Correr cd-hit para agrupar secuencias por identidad y cobertura

`cd-hit -i todos_na.faa -o todos70 -c 0.70 -n 4 -aL 0.7 -d 0 -M 3000 -T 2 > todos70.cdhit.out`

Convertir la salida a lista de clusters

`perl -pne 's/\t//g;s/^.*,//g;s/\.\.\..*$//g;s/\n/\t/g;s/\>Cluster\ /\n/g;s/\>//g; eof && do{chomp; print "$_ \n"; exit}' todos70.clstr > todos70.otu`

Obtener una lista de los archivos de mapeo de lecturas en marcos de lectura

`ls *.hits > nac_list.txt`

Obtener la tabla de presencia de clusters con la frecuencia convertida en mapeo de lecturas con el script hitter_na.py se obtiene el archivo *nac_na.tsv*

`python3 hitter_na.py nac_list.txt todos70.otu nac`
