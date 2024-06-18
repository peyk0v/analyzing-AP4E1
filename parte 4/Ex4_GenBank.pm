#!/usr/bin/perl

use strict;
use warnings;
use Bio::SearchIO;
use Bio::DB::GenBank;

# Leer los argumentos de la línea de comandos
my ($blast_file, $pattern) = @ARGV;
die "Uso: $0 <blast_output> <pattern>\n" unless @ARGV == 2;

print "Patrón: $pattern\n";
print "Archivo de BLAST: $blast_file\n";

open my $blast_fh, '<', $blast_file or die "No se pudo abrir el archivo $blast_file: $!";

my @hits;

my ($identifier, $specie);

while (my $line = <$blast_fh>) {
    chomp $line;

    if ($line =~ /^>(\S+)/) {
        if (defined $identifier && defined $specie) {
            if ($specie eq $pattern) {
                print "Accession: $identifier\n";
                print "Species: $specie\n";
                print "\n";
                push @hits, $identifier;
            }
        }

        $identifier = $1;
        $specie = undef;
    } elsif ($line =~ /\[(.*?)\]/) {
        $specie = $1;
    }
}

if (defined $identifier && defined $specie && $specie eq $pattern) {
    print "Accession: $identifier\n";
    print "Species: $specie\n";
    push @hits, $identifier;
}

close $blast_fh;

# Hasta acá tenemos todos los 'Accessions' que coinciden con el 'bast.out'

my $gb = Bio::DB::GenBank->new;

my $output_file = 'hits.fasta';

my $seqio_out = Bio::SeqIO->new(-file => ">$output_file", -format => 'fasta');

foreach my $hit (@hits) {
    my $seq = $gb->get_Seq_by_acc(@hits); # Estos son los accessions, por ejemplo Q8L7A9. Pero tira error, ya que no está en GenBank
    if ($seq) {
        $seqio_out->write_seq($seq);
    } else {
        warn "Failed to retrieve sequence for accession: $hit\n";
    }
}


print "Se han guardado las secuencias completas en $output_file.\n";
