#!/usr/bin/perl
use strict;
use warnings;
use Bio::SeqIO;

my $genbank_file = './AP4E1.gb';
my $output_file= './output.fasta';

my $seqio_obj = Bio::SeqIO->new(-file => $genbank_file, -format => "genbank");
my $fasta_out = Bio::SeqIO->new(-file => $output_file, -format => 'fasta');

while (my $seq_obj = $seqio_obj->next_seq) {
  # Traducción para marcos de lectura 0, 1 y 2
  for my $frame (0..2) {
    my $pro_seq_obj = $seq_obj->translate(-frame => $frame);
    $fasta_out->write_seq($pro_seq_obj);
  }

  # Traducción para marcos de lectura reverseados 
  my $rev_seq = $seq_obj->revcom;
  for my $frame (0..2) {
    my $pro_seq_obj = $rev_seq->translate(-frame => $frame);
    $fasta_out->write_seq($pro_seq_obj);
  }
}

print "Secuencias de aminoácidos traducidas y escritas en output.fasta\n";
