#!/usr/bin/env perl
#!/usr/bin/env perl

use strict;
use warnings;

use File::Basename;

# ========================================================================================
# --- Compute summary stats by taxon
# --- (Summary stats by gene output in exon stitching step.)
# --- Modified version of exon_stitching's summary_stats.pl
# ---   based on original script by Julie Allen, et al.
# ========================================================================================

# ----------------------------------------------------------------------------------------
#  Input:
#     Summary stats by gene file (from exon stitching step)
#         [out_dir]/final/summary_stats.per.taxon.csv
#  Output:
#     Summary stats by taxon file
#         [out_dir]/final/summary_stats.per.taxon.csv
# ----------------------------------------------------------------------------------------


my $exon_out_dir = shift;  # E.g. ../exon_stitching_schistosoma
my $exon_out_dir_final = "$exon_out_dir/final";

my @taxarray  = ();
my %taxhash   = ();
my %full      = ();

my %ninefive  = ();
my %ninezero  = ();
my %eightzero = ();
my %sevenzero = ();
my %fivezero  = ();
my %ten       = ();
my %lessten   = ();

my $header = "Taxon,NumGenes,FullExons,95%,90%,80%,70%,50%,10%\n";

open SUMSTATS, "> $exon_out_dir_final/summary_stats.per.taxon.csv";
print SUMSTATS $header;

print STDERR "Summary Stats Per Taxon\n";
print STDERR $header;

my $gene_stats = "$exon_out_dir_final/summary_stats.per.gene.csv";

open GENE, "< $gene_stats";

while (<GENE>) {

    if (/\S+?,(\S+?),(\d+?),(\d+)$/) {

        my $tax          = $1;
        my $querylength  = $2;
        my $targetlength = $3;

        my $proportion = $targetlength / $querylength;

        if (! exists $taxhash{$tax}) {
            $taxhash{$tax}   = 1;
            $full{$tax}      = 0;
            $ninefive{$tax}  = 0;
            $ninezero{$tax}  = 0;
            $eightzero{$tax} = 0;
            $sevenzero{$tax} = 0;
            $fivezero{$tax}  = 0;
            $ten{$tax}       = 0;
            $lessten{$tax}   = 0;
            push @taxarray, $tax;

        } else {
        	$taxhash{$tax}++;
        }

        if ($proportion == 1) {
        	$full{$tax}++;
        } elsif ($proportion >= 0.95 ) {
        	$ninefive{$tax}++;
        } elsif ($proportion >= 0.90 ) {
        	$ninezero{$tax}++;
        } elsif ($proportion >= 0.80 ) {
        	$eightzero{$tax}++;
        } elsif ($proportion >= 0.70 ) {
        	$sevenzero{$tax}++;
        } elsif ($proportion >= 0.50 ) {
        	$fivezero{$tax}++;
        } elsif ($proportion >= 0.10 ) {
        	$ten{$tax}++;
        } elsif ($proportion <  0.10 ) {
        	$lessten{$tax}++;
        }
    }
}

my $tot95 = 0;
my $tot90 = 0;
my $tot80 = 0;
my $tot70 = 0;
my $tot50 = 0;
my $tot10 = 0;
my $rest  = 0;

for my $tax (@taxarray) {

    $tot95 = ( $ninefive{$tax}  + $full{$tax} );
    $tot90 = ( $ninezero{$tax}  + $tot95      );
    $tot80 = ( $eightzero{$tax} + $tot90      );
    $tot70 = ( $sevenzero{$tax} + $tot80      );
    $tot50 = ( $fivezero{$tax}  + $tot70      );
    $tot10 = ( $ten{$tax}       + $tot50      );
    $rest  = ( $lessten{$tax}   + $tot10      );

    my $info = "$tax,$taxhash{$tax},$full{$tax},";
    $info = $info . "$tot95,$tot90,$tot80,$tot70,$tot50,$tot10\n";

    print SUMSTATS $info;
    print STDERR $info;
}

exit;
