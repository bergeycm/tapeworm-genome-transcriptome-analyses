#!/usr/bin/env perl

use strict;
use warnings;

use File::Basename;

# ========================================================================================
# --- Stich together contigs identified by exon_stitching_get_contigs.pl
# --- Modified version of exon_stitching's stitch.contigs.pl
# ---   based on original script by Julie Allen, et al.
# ========================================================================================

# ----------------------------------------------------------------------------------------
#  Input:
#     CSV files with info on which contigs to be stitched together
#         {gene}.overlap.{overlap}.contig_list.csv
#  Output:
#     Files of stitched together contigs from exonerate, separated by NNN
#         [out_dir]/final/[gene].stitched_exons.fasta
#     Summary file
#         [out_dir]/final/summary_stats.per.gene.csv
# ----------------------------------------------------------------------------------------

my $exon_out_dir = shift;  # E.g. ../exon_stitching_schistosoma
my $overlap      = shift;  # E.g. 10

my $verbose = 1;

my $exon_out_dir_final = "$exon_out_dir/final";
mkdir $exon_out_dir_final;

open SUMMARY, "> $exon_out_dir_final/summary_stats.per.gene.csv";
print SUMMARY "Locus,Taxon,Query_Length,Target_Length\n";

# ----------------------------------------------------------------------------------------
# --- Determine the number of loci for each locus
# ----------------------------------------------------------------------------------------

# Get list of input CSV files
my @csv_files =  <$exon_out_dir/*.contig_list.csv>;

foreach (@csv_files) {

    if (/(\S+).overlap.*.contig_list.csv/) {

        my $gene_file_full = $1;
        my $gene = basename($gene_file_full);
        if ($verbose) {
            print STDERR "Processing gene [$gene]\n";
        }

        # --- Print out FASTA file for each gene

        open FH,     "< $gene.overlap.$overlap.contig_list.csv";
        open OUT_FA, "> $exon_out_dir_final/$gene.stitched_exons.fasta";

        if ($verbose) {
            print STDERR "Writing to target ";
            print STDERR "[$exon_out_dir_final/$gene.stitched_exons.fasta]\n";
        }

        while (<FH>) {

            # --- Look for gene name

            if (/(\S+)/ && !/Inputfile/ && !/There\s+are/ && !/Number\s+of\s+Contigs/) {

                my $line = $1;
                print STDERR "\tLine from stats file: $line\n" if $verbose;
                my @gene_info = split (/,/, $line);

                # --- Get number of contigs

                my $infile             = $gene_info[1];
                my $querylength        = $gene_info[3];
                my $numcontigs         = $gene_info[4];
                my $combinedexonlength = $gene_info[6];

                print STDERR "\tNumber of contigs: $numcontigs\n" if $verbose;
                print STDERR "\tInput file: $infile\n" if $verbose;

                print SUMMARY "$gene,$infile,$querylength,$combinedexonlength\n";

                if ($numcontigs == 1) {

                    my $start  = $gene_info[7];
                    my $end    = $gene_info[8];
                    my $contig = $gene_info[9];

                    chomp $contig;

                    print STDERR "\tOne contig found for $contig\n" if $verbose;

                    $contig =~ s/,//g; 
                    $contig =~ s/\|//g; 

                    #$contig = " $contig";

                    my $seqflag = 0;
                    my $sequence = ();

                    open FASTA,  "< $gene.exons.fasta";
                    print OUT_FA ">$infile.$contig.$numcontigs\n";

                    print STDERR "\tWriting to $exon_out_dir_final/$infile.$contig.$numcontigs\n";   

                    while (my $line = <FASTA>) {
                        chomp $line;

                        # --- Remove interleavedness

                        if ($line =~ m/^>/) {

                            $seqflag = 0;

                            if ($line =~ /$infile/) {

                                my $linecontig = $line;
                                $linecontig =~ s/>$infile\,(.*),\d+$/$1/g;
                                $linecontig =~ s/,//g;
                                $linecontig =~ s/\|//g;
                                $linecontig = " $linecontig";

                                if ($linecontig =~ /$contig/) {
                                    $seqflag = 1;
                                }
                            }

                        } elsif ($seqflag == 1 ) {
                            $sequence = $sequence.$line;
                        }
                    }
                    close FASTA;

                    # --- Print out NNN at the beginning of the gene
                    # ---  if contig does not start at beginning

                    if ($start > 0) {
                        for (0..$start) {
                            print OUT_FA "NNN";
                        }
                    }

                    # --- Print gene sequence
                    print OUT_FA "$sequence";

                    # --- Print out NNN at the end of the gene
                    # ---  if contig does not cover the full length

                    if ($end < $querylength) {
                        for ($end .. ($querylength-1)) {
                            print OUT_FA "NNN";
                        }
                    }

                    print OUT_FA "\n";

                # --- Otherwise, if there is more than one contig,
                # --- stitch them together with NNN in between
                } else {

                    print STDERR "$infile has >1 contig, " if $verbose;
                    print STDERR "stitching together $numcontigs contigs\n" if $verbose;

                    my $count = 0;
                    my $contigstart = ($numcontigs * 2) + 7;
                    my $gapstart = 0;

                    # Start at 8 and go to 8 + number of contigs
                    for ($contigstart..($contigstart + $numcontigs - 1)) {

                        my $contignumber = $_;
                        print STDERR "contig number is $contignumber\n" if $verbose;

                        my $start = $gene_info[7 + $count];
                        my $end   = $gene_info[8 + $count];
                        my $contig = $gene_info[$contignumber];
                        #$contig = " $contig";
                        $contig =~ s/\|//g;

                        my $last = 0;
                        my $gapend = 0;
                        my $sequence = '';
                        my $seqflag = '';

                        # --- First contig

                        if ($contignumber == $contigstart) {

                            $gapstart = $end + 1;

                            open FASTA, "< $gene.exons.fasta";

                            while (my $line = <FASTA>) {

                                chomp $line;

                                if ($line =~ m/^>/) {

                                    $seqflag = 0;

                                    if ($line =~ /$infile/) {

                                        my $linecontig = $line;
                                        $linecontig =~ s/>$infile\,(.*),\d+$/$1/g;
                                        $linecontig = " $linecontig";
                                        $linecontig =~ s/\|//g;

                                        if ($linecontig =~ /$contig/) {
                                            $seqflag = 1;
                                            if ($verbose) {
                                                print STDERR "Writing header ";
                                                print STDERR "[$infile.$contig.$numcontigs]\n";
                                            }
                                            print OUT_FA ">$infile.$contig.$numcontigs\n";
                                        }
                                    }
                                } elsif ($seqflag == 1) {
                                    $sequence = $sequence.$line;
                                }
                            }
                            close FASTA;

                            # --- Print out NNN if necessary at the beginning
                            if ($start != 0) {
                                for (0..$start) {
                                    print OUT_FA "NNN";
                                }
                            }
                            print OUT_FA "$sequence";
                            my $nextstart = $end;
                            $gapstart = $end + 1;

                        # --- Last contig

                        } elsif ($contignumber == (($numcontigs - 1) + $contigstart)) {

                            # --- Make sure overlap is accounted for
                            for ($gapstart ..($start-1)) {
                                print OUT_FA "NNN";
                            }

                            open FASTA, "< $gene.exons.fasta";

                            while (my $line = <FASTA>) {

                                chomp $line;

                                if ($line =~ m/^>/) {

                                    $seqflag = 0;

                                    if ($line =~ /$infile/) {

                                        my $linecontig = $line;
                                        $linecontig =~ s/>$infile\,(.*),\d+$/$1/g;

                                        $linecontig = " $linecontig";
                                        $linecontig =~ s/\|//g;

                                        if ($linecontig =~ /$contig/) {
                                            $seqflag = 1;
                                        }
                                    }
                                } elsif ($seqflag == 1) {
                                    $sequence = $sequence.$line;
                                }
                            }
                            print OUT_FA "$sequence";

                            # --- Print out NNN at end if sequence doesn't span full length

                            if ($end < $querylength) {
                                for ($end .. ($querylength - 1)) {
                                    print OUT_FA "NNN";
                                }
                            }
                            print OUT_FA "\n";
                            close FASTA;

                        # --- Middle contigs

                        } elsif ($contignumber > $contigstart &&
                                 $contignumber < ($numcontigs-1 + $contigstart)) {

                            if ($verbose) {
                                print STDERR "Processing middle contig $contignumber\n";
                            }

                            # --- Take into account overlap

                            for ($gapstart ..($start-1)) {
                                print OUT_FA "NNN";
                            }

                            open FASTA, "< $gene.exons.fasta";

                            $gapstart = $end + 1;

                            while (my $line = <FASTA>) {

                                chomp $line;

                                if ($line =~ m/^>/) {

                                    $seqflag = 0;

                                    if ($line =~ /$infile/) {

                                        my $linecontig = $line;
                                        $linecontig =~ s/>$infile\,(.*),\d+$/$1/g;

                                        $linecontig = " $linecontig";
                                        $linecontig =~ s/\|//g;

                                        if ($linecontig =~ /$contig/) {
                                            $seqflag = 1;
                                        }
                                    }
                                } elsif ($seqflag == 1) {
                                    $sequence = $sequence.$line;
                                }
                            }

                            print OUT_FA "$sequence";
                        }

                        $count = $count + 2;
                    }
                }
           }
       }
   }
}

exit
