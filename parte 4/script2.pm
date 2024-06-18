#!/usr/bin/perl

use strict;
use warnings;
use Bio::SearchIO;
use Bio::DB::SwissProt;
use Bio::SeqIO;

my ($blast_file, $pattern) = @ARGV;
die "Uso: $0 <blast_output> <pattern>\n" unless @ARGV == 2;

print "Patrón: $pattern\n";
print "Archivo de BLAST: $blast_file\n";

my $searchio = Bio::SearchIO->new(-format => 'blast', -file => $blast_file);

my $sp = Bio::DB::SwissProt->new(); # La razón por la que 'funciona' es porque use SwissProt y no GenBank

my $output_file = 'hits.fasta';
my $seq_out = Bio::SeqIO->new(-file => ">$output_file", -format => 'fasta');

my @hits;

while (my $result = $searchio->next_result) {
    while (my $hit = $result->next_hit) {
        if ($hit->description =~ /$pattern/i) {
            my $accession = $hit->accession;
            print "Encontrado: $accession - ", $hit->description, "\n";
            push @hits, $accession;
        }
    }
}

foreach my $hit (@hits) {
    my $seq;
    eval { $seq = $sp->get_Seq_by_acc($hit) };
    if ($@) {
        warn "Error obteniendo la secuencia para $hit: $@\n";
        next;
    }
    $seq_out->write_seq($seq);
}

print "Se han guardado las secuencias completas en $output_file.\n";
