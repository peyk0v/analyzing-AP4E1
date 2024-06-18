use strict;
use warnings;
use Bio::SeqIO;

# Nombre del archivo GenBank
my $file = 'AP4E1.gb';
my $output_prefix = 'output_';

# Crear un objeto SeqIO para leer el archivo GenBank
my $seqio = Bio::SeqIO->new(-file => $file, -format => 'genbank');

# Iterar sobre las secuencias en el archivo GenBank
while (my $seq = $seqio->next_seq) {
    my $seq_id = $seq->display_id;
    my $cds_count = 0;
    
    for my $feat_object ($seq->get_SeqFeatures) {
        if ($feat_object->primary_tag eq 'CDS') {
            my $cds_seq = $feat_object->spliced_seq->seq; # Obtener la secuencia de nucleótidos del CDS

            # Generar las seis posibles secuencias de aminoácidos
            for my $frame (0..2) {
                my $protein_seq1 = translate_frame($cds_seq, $frame);
                my $protein_seq2 = translate_frame(reverse_complement($cds_seq), $frame);

                my $title1 = ">Aminoacid Translation of ${seq_id} frame ${frame} (direct)";
                my $title2 = ">Aminoacid Translation of ${seq_id} frame ${frame} (reverse)";

                create_fasta_file("${output_prefix}${seq_id}_CDS${cds_count}_frame${frame}_1.fasta", $title1, $protein_seq1);
                create_fasta_file("${output_prefix}${seq_id}_CDS${cds_count}_frame${frame}_2.fasta", $title2, $protein_seq2);
            }
            
            $cds_count++;
        }
    }
}

sub translate_frame {
    my ($seq, $frame) = @_;
    my %codons_table = (
        'TTT' => 'F', 'TTC' => 'F',    # Fenilalanina (Phe)
        'TTA' => 'L', 'TTG' => 'L',    # Leucina (Leu)
        'CTT' => 'L', 'CTC' => 'L', 'CTA' => 'L', 'CTG' => 'L',    # Leucina (Leu)
        'ATT' => 'I', 'ATC' => 'I', 'ATA' => 'I',    # Isoleucina (Ile)
        'ATG' => 'M',    # Metionina (Met)
        'GTT' => 'V', 'GTC' => 'V', 'GTA' => 'V', 'GTG' => 'V',    # Valina (Val)
        'TCT' => 'S', 'TCC' => 'S', 'TCA' => 'S', 'TCG' => 'S',    # Serina (Ser)
        'CCT' => 'P', 'CCC' => 'P', 'CCA' => 'P', 'CCG' => 'P',    # Prolina (Pro)
        'ACT' => 'T', 'ACC' => 'T', 'ACA' => 'T', 'ACG' => 'T',    # Treonina (Thr)
        'GCT' => 'A', 'GCC' => 'A', 'GCA' => 'A', 'GCG' => 'A',    # Alanina (Ala)
        'TAT' => 'Y', 'TAC' => 'Y',    # Tirosina (Tyr)
        'TAA' => '*', 'TAG' => '*', 'TGA' => '*',    # Codones de parada
        'CAT' => 'H', 'CAC' => 'H',    # Histidina (His)
        'CAA' => 'Q', 'CAG' => 'Q',    # Glutamina (Gln)
        'AAT' => 'N', 'AAC' => 'N',    # Asparagina (Asn)
        'AAA' => 'K', 'AAG' => 'K',    # Lisina (Lys)
        'GAT' => 'D', 'GAC' => 'D',    # Ácido aspártico (Asp)
        'GAA' => 'E', 'GAG' => 'E',    # Ácido glutámico (Glu)
        'TGT' => 'C', 'TGC' => 'C',    # Cisteína (Cys)
        'TGG' => 'W',    # Triptófano (Trp)
        'CGT' => 'R', 'CGC' => 'R', 'CGA' => 'R', 'CGG' => 'R',    # Arginina (Arg)
        'AGT' => 'S', 'AGC' => 'S',    # Serina (Ser)
        'AGA' => 'R', 'AGG' => 'R',    # Arginina (Arg)
        'GGT' => 'G', 'GGC' => 'G', 'GGA' => 'G', 'GGG' => 'G'     # Glicina (Gly)
    );

    my $translation = '';
    for (my $i = $frame; $i < length($seq) - 2; $i += 3) {
        my $codon = substr($seq, $i, 3);
        my $aa = exists $codons_table{$codon} ? $codons_table{$codon} : '-';
        $translation .= $aa;
    }
    return $translation;
}

sub reverse_complement {
    my ($dna_seq) = @_;
    $dna_seq = reverse($dna_seq);
    $dna_seq =~ tr/ACGTacgt/TGCAtgca/;
    return $dna_seq;
}

sub create_fasta_file {
    my ($file, $title, $seq) = @_;
    open(my $fh, '>', $file) or die "No se pudo abrir el archivo '$file' $!";
    print $fh "$title\n";
    print $fh "$seq\n";
    close $fh;
    print "Archivo '$file' creado con éxito.\n";
}
