#!/usr/bin/env perl

# ========================================================================================
# --- Edit aTRAM output files and make gene list
# --- Modified version of exon_stitching's editfasta.pl
# ---   based on original script by Julie Allen, et al.
# ========================================================================================

# ----------------------------------------------------------------------------------------
#  Input:
#     assembly files from aTRAM
#     list of reference genes
#  Output:
#     cleaned FASTA files, *.ed.fasta, in same passed $fasta_dir directory
#     FASTA files for all reference genes, *.reference.fasta
# ----------------------------------------------------------------------------------------

use strict;
use warnings;

use File::Glob;
use Parallel::ForkManager;

my $fasta_dir      = shift;  # E.g. ../results/atram/atram_output_$PROTEOME_NAME/combined/
my $gene_list_out  = shift;  # E.g. ../results/atram/atram_output_$PROTEOME_NAME/gene_list.txt
my $proteome_query = shift;  # E.g. ../genomes/schistosoma_mansoni/[...].fa_PART000
chomp $proteome_query;

# Delete empty FASTA files
`find $fasta_dir -name "*.fasta" -empty -type f -delete`;

# --- Create version of each file without weird characters, add unique number to each contig

# Get list of FASTA files
my @fasta_files =  <$fasta_dir/*>;

my $pm = Parallel::ForkManager->new(20);

sub process_fasta {
    my ($file) = @_;

    if ($file =~ /(\S+).fasta/) {

        my $fasta_prefix = $1;
        print STDERR "Processing FASTA file with prefix [$fasta_prefix]\n";

        open IN_FA,  "< $fasta_prefix.fasta";
        open OUT_FA, "> $fasta_prefix.ed.fasta";

        my $num = 0;
        while (<IN_FA>) {
            if (/^>(.*)$/) {
                my $contig = $1;
                $num++;
                $contig =~ s/[-,=@+\[\]:!]//g;
                $contig =~ s/\s+//g;
                print OUT_FA ">$contig\_$num\n";
            } elsif (/^\S+/ && ! /Command|Hostname|completed/) {
            	print OUT_FA;
            }
        }
    }
}

foreach my $file (@fasta_files) {
    my $pid = $pm->start and next;
    process_fasta($file);
    $pm->finish(0, { input => $file });
}
$pm->wait_all_children;


# --- Create gene list and make copy of reference genes' FASTA files

open GENE_LIST,  "> $gene_list_out";
open PROT_QUERY, "< $proteome_query";

while (<PROT_QUERY>) {

    if (/^>(.*)/) {
        my $gene = $1;
        $gene =~ s/\s+/_/g;
        open GENE, ">$gene.reference.fasta";
        print GENE_LIST "$gene.reference.fasta\n";
        print GENE ">$gene\n";
    } else {
    	print GENE;
    }
}

exit;
