#!/usr/bin/perl
use strict;
use warnings;
use Bio::SeqIO;

my $genbank_file = 'AP4E1.gb';

my $seqio_obj = Bio::SeqIO->new(-file => $genbank_file, -format => "genbank");

my $fasta_out = Bio::SeqIO->new(-file => '>output.fasta', -format => 'fasta');

while (my $seq_obj = $seqio_obj->next_seq) {
    #para los frames 0, 1 y 2
	for my $frame (0..2) {
    	my $pro_seq_obj = $seq_obj->translate(-frame => $frame);
    	$fasta_out->write_seq($pro_seq_obj);
	}

    #para los frames -1 y -2
    my $rev_seq = $seq_obj->revcom;
    for my $frame (0..2) {
        my $pro_seq_obj = $rev_seq->translate(-frame => $frame);
        $fasta_out->write_seq($pro_seq_obj);
    }
}

print "Secuencias de aminoácidos traducidas\n";
