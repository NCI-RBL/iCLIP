use strict;
use Bio::SeqIO;
my %unique;

#read args from command line
my $file = $ARGV[0];
my $output_file = $ARGV[1];

#file handling
my $seqio = Bio::SeqIO->new(-file => $file, -format => "fasta");
my $outseq = Bio::SeqIO->new(-file => ">$output_file.uniq", -format => "fasta");

while(my $seqs = $seqio->next_seq) {
    my $id = $seqs->display_id;
    my $seq = $seqs->seq;
    unless(exists($unique{$seq})) {
        $outseq->write_seq($seqs);
        $unique{$seq} +=1;
    }
}