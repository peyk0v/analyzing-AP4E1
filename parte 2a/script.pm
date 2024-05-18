#!/usr/bin/perl
use strict;
use warnings;

my $input_file = 'input.fasta';

my $output_file = 'blast.out';

my $blastp_cmd = "blastp -query $input_file -db ../../../Blast/ncbi-blast-2.15.0+/data/swissprot -out $output_file -outfmt 7";

system($blastp_cmd);

print "BLAST completado. Resultados guardados en $output_file\n";
