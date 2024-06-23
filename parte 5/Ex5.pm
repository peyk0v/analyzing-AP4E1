#!/usr/bin/perl

use strict;
use warnings;
use Bio::Factory::EMBOSS;
use Bio::SeqIO;

my ($seq_nucleotidos_file, $aminoacid_file) = @ARGV;
die "Uso: $0 <seq_nucleotidos.fasta> <aminoacid.fasta> \n" unless @ARGV == 2;

my $output_file = "orfs.fasta";
my $output_domains = "dominios.txt";

# Generate ORFs from nucleotids file
my $getorf_cmd = "getorf -sequence $seq_nucleotidos_file -outseq $output_file";
system($getorf_cmd) == 0 or die "Error ejecutando getorf: $!";
print "ORFs obtenidos correctamente en $output_file\n";

# Generate domains from aminoacid seq (protein)
my $factory = Bio::Factory::EMBOSS->new();
my $prosextract = $factory->program('prosextract');
$prosextract->run({});

my $patmatmotifs = $factory->program('patmatmotifs');
$patmatmotifs->run({-sequence => $aminoacid_file, -outfile => $output_domains, -full => 'Y'});
print "Dominios funcionales encontrados en $output_domains\n";
