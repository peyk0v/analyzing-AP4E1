#!/usr/bin/perl

use strict;
use warnings;

my ($input_file) = @ARGV;
die "Uso: $0 <seq_nucleotidos> \n" unless @ARGV == 1;

my $output_file = "orfs.fasta";
my $output_domains = "dominios.txt";

# Generate ORFs
my $getorf_cmd = "getorf -sequence $input_file -outseq $output_file";
system($getorf_cmd) == 0 or die "Error ejecutando getorf: $!";
print "ORFs obtenidos correctamente en $output_file\n";

# Generate domains
my $patmatmotifs_cmd = "patmatmotifs -sequence $input_file -outfile $output_domains -full";
system($patmatmotifs_cmd) == 0 or die "Error ejecutando patmatmotifs: $!";
print "Dominios funcionales encontrados en $output_domains\n";
