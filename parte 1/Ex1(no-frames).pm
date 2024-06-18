use strict;
use warnings;
use Bio::SeqIO;

# Nombre del archivo GenBank
my $file = 'AP4E1.gb';

my $protein_file = 'mRNA_protein_encoding.fasta';

my $aminoacid_file = 'protein_aminoacid.fasta';

# Crear un objeto SeqIO para leer el archivo GenBank
my $seqio = Bio::SeqIO->new(-file => $file, -format => 'genbank');

# Iterar sobre las características (features) en busca de CDS
while (my $seq = $seqio->next_seq) {
    for my $feat_object ($seq->get_SeqFeatures) {
        if ($feat_object->primary_tag eq 'CDS') {
            my $cds_seq = $feat_object->spliced_seq->seq; # Obtener la secuencia de nucleótidos del CDS
            my $rna_seq = transcribe($cds_seq); # Transcribir la secuencia de nucleótidos a ARN
            createFile("Encoding Protein of AP4E1", $rna_seq, $protein_file);
            my $protein_seq = translate($rna_seq); # Traducir la secuencia de ARN a aminoácidos
            createFile("Aminoacid Translation of AP4E1's mRNA", $protein_seq, $aminoacid_file);
            print "CDS Translation: $protein_seq\n";
        }
    }
}

sub createFile {
    my ($title, $seq, $file) = @_;
    open(my $fh, '>', $file) or die "Cannot open file $file: $!";
    print $fh ">$title\n$seq\n";
    close $fh or die "Cannot close file $file: $!";
}

# Subrutina para transcribir una secuencia de ADN a ARN
sub transcribe {
    my ($dna_seq) = @_;
    $dna_seq =~ tr/T/U/;
    return $dna_seq;
}

# Subrutina para traducir una secuencia de ARN a aminoácidos
sub translate {
    my ($rna_seq) = @_;

    my %codons_table = (
        'UUU' => 'F', 'UUC' => 'F',    # Fenilalanina (Phe)
        'UUA' => 'L', 'UUG' => 'L',    # Leucina (Leu)
        'CUU' => 'L', 'CUC' => 'L', 'CUA' => 'L', 'CUG' => 'L',    # Leucina (Leu)
        'AUU' => 'I', 'AUC' => 'I', 'AUA' => 'I',    # Isoleucina (Ile)
        'AUG' => 'M',    # Metionina (Met)
        'GUU' => 'V', 'GUC' => 'V', 'GUA' => 'V', 'GUG' => 'V',    # Valina (Val)
        'UCU' => 'S', 'UCC' => 'S', 'UCA' => 'S', 'UCG' => 'S',    # Serina (Ser)
        'CCU' => 'P', 'CCC' => 'P', 'CCA' => 'P', 'CCG' => 'P',    # Prolina (Pro)
        'ACU' => 'T', 'ACC' => 'T', 'ACA' => 'T', 'ACG' => 'T',    # Treonina (Thr)
        'GCU' => 'A', 'GCC' => 'A', 'GCA' => 'A', 'GCG' => 'A',    # Alanina (Ala)
        'UAU' => 'Y', 'UAC' => 'Y',    # Tirosina (Tyr)
        'UAA' => '*', 'UAG' => '*', 'UGA' => '*',    # Codones de parada
        'CAU' => 'H', 'CAC' => 'H',    # Histidina (His)
        'CAA' => 'Q', 'CAG' => 'Q',    # Glutamina (Gln)
        'AAU' => 'N', 'AAC' => 'N',    # Asparagina (Asn)
        'AAA' => 'K', 'AAG' => 'K',    # Lisina (Lys)
        'GAU' => 'D', 'GAC' => 'D',    # Ácido aspártico (Asp)
        'GAA' => 'E', 'GAG' => 'E',    # Ácido glutámico (Glu)
        'UGU' => 'C', 'UGC' => 'C',    # Cisteína (Cys)
        'UGG' => 'W',    # Triptófano (Trp)
        'CGU' => 'R', 'CGC' => 'R', 'CGA' => 'R', 'CGG' => 'R',    # Arginina (Arg)
        'AGU' => 'S', 'AGC' => 'S',    # Serina (Ser)
        'AGA' => 'R', 'AGG' => 'R',    # Arginina (Arg)
        'GGU' => 'G', 'GGC' => 'G', 'GGA' => 'G', 'GGG' => 'G'     # Glicina (Gly)
    );
    my $translation = '';
    for (my $i = 0; $i < length($rna_seq) - 2; $i += 3) {
        my $codon = substr($rna_seq, $i, 3);
        my $aa = exists $codons_table{$codon} ? $codons_table{$codon} : '-';
        $translation .= $aa;
    }
    return $translation;
}
