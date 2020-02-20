Get a table of annotated and unannotated genes with counts based in read mapping. 
The scripts are included in this repository.

## Get gene abundance tables

Starting files:

*acacia.faa* protein sequences. 
(The names of the sequences must have the sample name as prefix and the sequence number as suffix, separated by an underscore.)

*acacia.fna* CDS sequences.
(The names of the sequences must have the sample name as prefix and the sequence number as suffix, separated by an underscore. The sequence names must correspond with the protein sequences file.)

*acacia-rn.fasta* reads.

### Abundance
Create the CDS sequence database with bowtie2.

`bowtie2-build acacia.fna acacia`

Map the reads to the CDS .

`bowtie2 -f -x acacia -U acacia-rn.fasta -S acacia.sam --quiet -p 20 --very-sensitive`

Recover only unambiguous alignments and get the mapped CDS with their respective number of hits.

`grep -v '^@' acacia.sam | awk '{if($5 == "42") print $3}' | sort | uniq -c > acacia.hits`

### Annotation

Split the protein multi fasta to speed up the annotation process.

`partefasta 10000 acacia90.faa`

Search the protein sequences against the M5NR database.

`diamond blastp -d m5nr -q acacia90.faa.1.fas -f 6 -e 1e-10 -k 10 -p 1 --quiet -o acacia90.faa.1.fas.bout`

Join all search tables in a sigle file *acacia.faa*

`cat acacia90.faa.*.bout > acacia-m5nr.bout`

Oneliner to sort the hits by bitscore and remove bitscore duplicates.

`
cat acacia-m5nr.bout | perl -pe ' $name_col=0; $score_col=11; while(<>) { s/\r?\n//; @F=split /\t/, $_; ($n, $s) = @F[$name_col, $score_col]; if (! exists($max{$n})) { push @names, $n }; if (! exists($max
{$n}) || $s > $max{$n}) { $max{$n} = $s; $best{$n} = () }; if ($s == $max{$n}) { $best{$n} .= "$_\n" }; } for $n (@names) { print $best{$n} } ' >best;  perl -e ' $column=0; $unique=0; while(<>) { s/\r?\n//; @F=s
plit /\t/, $_; if (! ($save{$F[$column]}++)) { print "$_\n"; $unique++ } } ' best > acacia_best_uniq; rm best
`

Get only relevant columns.

`awk '{print $1"\t"$3"\t"$11"\t"$12"\t"$2}' acacia_best_uniq > acacia_best.simple.tsv`

### Build tables

#### Important: In order to use the following scripts (hitter.py, hitter_table.py y hitter_na.py) the sequences must have the sample name as prefix and the mapped sequence abundance must have the sample name. The sample names must be consistent and should not have underscores.

E.g. the mapped sequence file *acacia.hits* should look like this:

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

Get a table containing a M5NR id in each row and its mapped read abundance in the second column.
The resulting file is named *acacia.hout*.

`python3 hitter.py acacia_best.simple.tsv acacia.hits acacia`

Create a list of all table files.

`ls *.hout > lista`

Join all samples in a single table with the hitter_table.py script. The resulting file is named *nac.tsv*.

`python3 hitter_table.py lista nac`

### Add clusters of unannotates sequences

Get a list of annotated sequences.

`awk '{print $1}' acacia_best.simple.tsv > acacia_anotados.txt`

Get a list of all sequences.

`grep '>' acacia.faa | sed 's/>//g' > acacia_todos.txt`

Extract the names of the unannotated sequences.

`cat acacia_anotados.txt acacia_todos.txt | sort | uniq -c | grep '1 ' | awk '{print $2}' > acacia_na.txt`

Recover the unannotated sequences.

`seqtk subseq acacia.faa acacia_na.txt > acacia_na.faa`

Join the multi fasta files of all unannotated sequences from all samples.

`cat *_na.faa > todos_na.faa`

Cluster these sequences using CD-HIT.

`cd-hit -i todos_na.faa -o todos70 -c 0.70 -n 4 -aL 0.7 -d 0 -M 3000 -T 2 > todos70.cdhit.out`

Transform the cluster file to a list.

`perl -pne 's/\t//g;s/^.*,//g;s/\.\.\..*$//g;s/\n/\t/g;s/\>Cluster\ /\n/g;s/\>//g; eof && do{chomp; print "$_ \n"; exit}' todos70.clstr > todos70.otu`

Get a list of the mapped read abundance files.

`ls *.hits > nac_list.txt`

Build a table with a protein cluster in each row and its mapped read abundance along all samples as columns. The resulting fila is named *nac_na.tsv*

`python3 hitter_na.py nac_list.txt todos70.otu nac`

Get the table with both M5NR annotated proteins and unannotated protein clusters and their mapped read abundance in each sample.

`cat nac.tsv nac_na.tsv > nac_complete.tsv`
