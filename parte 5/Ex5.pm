#!/usr/bin/perl

use strict;
use warnings;

my ($seq_nucleotidos_file, $aminoacid_file) = @ARGV;
die "Uso: $0 <seq_nucleotidos.fasta> <aminoacid.fasta> \n" unless @ARGV == 2;

my $output_file = "orfs.fasta";
my $output_domains = "dominios.txt";

# Generate ORFs from nucleotids file
my $getorf_cmd = "getorf -sequence $seq_nucleotidos_file -outseq $output_file";
system($getorf_cmd) == 0 or die "Error ejecutando getorf: $!";
print "ORFs obtenidos correctamente en $output_file\n";

# Generate domains from aminoacid seq (protein)
my $patmatmotifs_cmd = "patmatmotifs -sequence $aminoacid_file -outfile $output_domains -full";
system($patmatmotifs_cmd) == 0 or die "Error ejecutando patmatmotifs: $!";
print "Dominios funcionales encontrados en $output_domains\n";
