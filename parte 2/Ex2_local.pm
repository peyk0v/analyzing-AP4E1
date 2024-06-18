#!/usr/bin/perl
use strict;
use warnings;

my ($input_file) = @ARGV; #aminoacid.fasta
die "Uso: $0 <input_file> \n" unless @ARGV == 1;

my $output_file = 'local_blast.out';

my $blastp_cmd = "blastp -query $input_file -db ../../../Blast/ncbi-blast-2.15.0+/data/swissprot -out $output_file";

system($blastp_cmd);

print "BLAST completado. Resultados guardados en $output_file\n";
