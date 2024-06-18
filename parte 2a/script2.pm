#!/usr/bin/perl
use strict;
use warnings;

my $input_file = 'input.fasta';

my $output_file = 'blast2.out';

my $blastp_cmd = "blastp -query $input_file -db ../../../BioInformatica/Blast/ncbi-blast-2.15.0+/data/swissprot -out $output_file";

system($blastp_cmd);

print "BLAST completado. Resultados guardados en $output_file\n";
