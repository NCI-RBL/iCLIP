use strict;
use Bio::SeqIO;
my %unique;
my $file = "gencode.v32.chr_patch_hapl_scaff.annotation.gtf.SplicedTransc.Rad46Num100kMin10Splices.fasta";
my $seqio = Bio::SeqIO->new(-file => $file, -format => "fasta");
my $outseq = Bio::SeqIO->new(-file => ">$file.uniq", -format => "fasta");

while(my $seqs = $seqio->next_seq) {
    my $id = $seqs->display_id;
    my $seq = $seqs->seq;
    unless(exists($unique{$seq})) {
        $outseq->write_seq($seqs);
        $unique{$seq} +=1;
    }
}

use strict;
use Bio::SeqIO;
my %unique;
my $file = "gencode.v32.chr_patch_hapl_scaff.annotation.gtf.SplicedTransc.Rad46Num100kMin10Transcripts.fasta";
my $seqio = Bio::SeqIO->new(-file => $file, -format => "fasta");
my $outseq = Bio::SeqIO->new(-file => ">$file.uniq", -format => "fasta");

while(my $seqs = $seqio->next_seq) {
    my $id = $seqs->display_id;
    my $seq = $seqs->seq;
    unless(exists($unique{$seq})) {
        $outseq->write_seq($seqs);
        $unique{$seq} +=1;
    }
}

