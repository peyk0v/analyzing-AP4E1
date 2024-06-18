#!/usr/bin/perl
use strict;
use warnings;
use Bio::SeqIO;
use Bio::Tools::Run::RemoteBlast;

# Configuración de BLAST remoto
my $prog   = 'blastp';         # Tipo de BLAST que quieres ejecutar
my $db     = 'swissprot';      # Base de datos a usar
my $outfile = 'remote_results.out';  # Nombre del archivo de salida

# Crear objeto de RemoteBlast
my $factory = Bio::Tools::Run::RemoteBlast->new(
    -prog       => $prog,
    -data       => $db,
    -outfile    => $outfile,
);

# Leer la secuencia desde un archivo fasta
my $seqio = Bio::SeqIO->new(-file => 'protein_aminoacid.fasta', -format => 'fasta');
my $seq_obj = $seqio->next_seq;  # Obtener la primera secuencia del archivo

# Ejecutar BLAST para la secuencia obtenida
$factory->submit_blast($seq_obj);
print "Submitted BLAST job\n";

# Esperar hasta que el BLAST se complete
my @rids = $factory->each_rid;
while (my $rid = shift @rids) {
    my $rc = $factory->retrieve_blast($rid);
    if (!ref($rc)) {
        if ($rc < 0) {
            $factory->remove_rid($rid);
        }
        print "Waiting for results... (", $rid, ")\n";
        sleep 5;
        push(@rids, $rid);  # Volver a agregar el RID para seguir esperando
    } else {
        my $result = $rc->next_result;
        my $filename = $result->query_name . "_blast_result.xml";
        $factory->save_output($filename);
        $factory->remove_rid($rid);
        print "BLAST search completed for ", $result->query_name, "\n";
    }
}

print "BLAST completado.\n";
