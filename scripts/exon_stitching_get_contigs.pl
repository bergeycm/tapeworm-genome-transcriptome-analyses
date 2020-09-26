#!/usr/bin/env perl

use strict;
use warnings;

use File::Glob;

# ========================================================================================
# --- Figure out which contigs (grabbed by exonerate and exon_stitching_edit_fastas.pl)
# ---   need to be stitched together
# --- Modified version of exon_stitching's getcontigs.pl
# ---   based on original script by Julie Allen, et al.
# ========================================================================================

# ----------------------------------------------------------------------------------------
#  Input:
#     output files from exonerate and exon_stitching_edit_fastas.pl script
#  Output:
#     CSV files with info on which contigs to be stitched together
#         {gene}.overlap.{overlap}.contig_list.csv
# ----------------------------------------------------------------------------------------

# Notes:
# - Contigs need to be in order they appear in the gene,
#    so sorted by the beginning of the contig.
# - Input will be all files in passed folder ($exon_out_dir)
#    that match *results.sorted.csv
# - For each input file, this creates $input.stats.csv
#    containing information about the contigs from exonerate

my $exon_out_dir = shift;   # E.g. exon_stitching_schistosoma/
my $overlap      = shift;   # E.g. 10

my $verbose = 1;

# Get list of input CSV files
my @csv_files =  <$exon_out_dir/*results.sorted.csv>;

foreach (@csv_files) {

    if (/(\S+).results.sorted.csv/) {

        my $inputfile = $1;
        my $gene = $inputfile;
        $gene =~ s/\S+\///g;
        print STDERR "Processing gene [$gene]\n" if $verbose;

        # --------------------------------------------------------------------------------
        # --- Create output file for contig names and where to stitch
        # --------------------------------------------------------------------------------

        my $statoutput = "$inputfile.overlap.$overlap.contig_list.csv";
        open STATS, "> $statoutput";
        my $date = localtime(time);
        print STATS "Inputfile\t$inputfile.csv   ";
        print STATS "Allowing overlap $overlap\n";
        my @taxarray = ();
        my $counttax = 0;
        my %taxhash = ();

        # --------------------------------------------------------------------------------
        # --- Get list of taxa from file
        # --------------------------------------------------------------------------------

        open FH, "< $inputfile.results.sorted.csv";
        while (<FH>) {
            if (/^\S+?,(\S+?),/) {
                my $tax=$1;
                if (! exists $taxhash{$tax}) {
                    $taxhash{$tax} = 1;
                    push @taxarray, $tax;
                    $counttax++;
                }
            }
        }

        print STATS "There are $counttax Libraries in $inputfile\n";
        print STATS "Gene,Taxon,Number of Contigs,GeneLength,";
        print STATS "Contigs to Keep,Total Overlap,Combined Exon Length,";
        print STATS "Beginning,End,Beginning,End,ContigName\n";

        print STDERR "There are $counttax Libraries in $inputfile\n" if $verbose;
        close FH;

        # --------------------------------------------------------------------------------
        # --- For each taxon in the file loop through and find the best contig(s)
        # --------------------------------------------------------------------------------

        for my $tax (@taxarray) {

            print STATS "$gene,$tax,";

            my @contigarray = ();
            my $numcontigs = 0;
            my @bigarray = ();
            my %contighash = ();

            # --- Count the number of contigs for each taxon,
            # ---  push the line to array of arrays, and contig into an array.

            open FH1, "< $inputfile.results.sorted.csv";
            my $pos = 0;
            while (<FH1>) {
                if (/\S+?,$tax,(\S+?,\S+?,\S+?,\S+?),(.*)$/) {
                    my $line   = $1;
                    my $contig = $2;

                    push @contigarray, $contig;
                    $numcontigs++;
                    my @linearray = split(/,/, $line);
                    for my $each (@linearray) {
                        push @{$bigarray[$pos]}, $each;
                    }
                    $pos++;
                }
            }

            print STDERR "The number of contigs is $numcontigs\n" if $verbose;
            print STDERR "Number of contigs $numcontigs, "        if $verbose;
            print STDERR "length of full gene $bigarray[0][0]\n"  if $verbose;

            close FH1;

            # Print the number of contigs and the length of the full gene
            print STATS "$numcontigs,$bigarray[0][0],";

            # --- Get statistics on each contig and determine if there is a full gene

            my $flag = 0;
            my ($beg, $end, $querylength, $queryoverlap) = (0) x 4;

            for my $cont (0..($numcontigs-1)) {

                print STDERR "contig number $cont\n" if $verbose;

                my $targetlength = $bigarray[$cont][0];
                print STDERR "target exon length $targetlength\n" if $verbose;

                $querylength = $bigarray[$cont][1];
                print STDERR "query exon length $querylength\n" if $verbose;

                # --- If there is a contig with the whole length of the gene save it
                if ($targetlength == $querylength) {
                    print STDERR "Target length == query length, full gene found\n" if $verbose;
                    $contighash{$tax} = $contigarray[$cont];
                    $flag = 1;
                    $beg = $bigarray[$cont][2];
                    print STDERR "\tBeginning location $beg\n" if $verbose;
                    $end = $bigarray[$cont][3];
                    print STDERR "\tEnding location $end\n" if $verbose;
                }
            }

            # --- If we have a full gene, print stats

            if ($flag == 1 ) {

                print STATS "1,$queryoverlap,$querylength,$beg,$end,$contighash{$tax}";

            } else {

                # No contig found that is the whole length

                my ($sum, $numlastcontig,  $numtokeep, $lastcontig, $totoverlap) = (0) x 5;
                my (@begarray,  @endarray, @keepcontigarray, %keepcontighashBeg, %keepcontighashEnd);
                my $total = $numcontigs - 1;

                # --- Process for when number of contigs == 1
                if ($numcontigs == 1) {

                    print STDERR "\tOne contig found, but not long enough to be full gene\n" if $verbose;

                    $numtokeep = 1;

                    for my $pos (0..($total)) {

                        $sum = $bigarray[$pos][1];

                        my $contig = $contigarray[$pos];
                        push @keepcontigarray, $contig;

                        my $beg = $bigarray[$pos][2];
                        push @begarray, $beg;
                        $keepcontighashBeg{$contig} = $beg;

                        my $end = $bigarray[$pos][3];
                        push @endarray, $end;
                        $keepcontighashEnd{$contig} = $end;

                        print STATS "$numcontigs,$queryoverlap,$sum,$beg,$end,$contig";

                        print STDERR "$tax one contig total $sum beginning $beg end $end contig $contig\n" if $verbose;
                    }
                }

                # --- Process for when number of contigs > 1

                my @lengtharray;

                if ($numcontigs > 1) {

                    print STDERR "\tFor $tax, number of contigs is >1. $numcontigs contigs found\n" if $verbose;

                    my @comparecontigarray = ();

                    for $pos (0..($total)) {

                        my $contig = $contigarray[$pos];
                        print STDERR "\t\tContig is $contig\t" if $verbose;

                        my $beg = $bigarray[$pos][2];
                        print "\t\tBeginning is $beg\t" if $verbose;

                        my $end = $bigarray[$pos][3];
                        print "\t\tEnd is $end\n" if $verbose;

                        print STDERR "\t\tLength is $bigarray[$pos][1]\n" if $verbose;

                        push @begarray, $beg;
                        push @endarray, $end;
                        push @comparecontigarray, $contig;
                        push @lengtharray, $bigarray[$pos][1];
                    }

                    my $pos = 0;
                    my $next = 1;

                    for (0..$total) {
                        if ($next <= ($total)) {

                            my $nextbeg = $begarray[$next];
                            print STDERR "beginning is $nextbeg\t" if $verbose;
                            my $endpos = $endarray[$pos];
                            my $nextend = $endarray[$next];
                            print STDERR "end is $nextend\n" if $verbose;

                            if ($nextbeg > ($endpos-$overlap) && ($nextend > $endpos)) {
                                print "Does not overlap.\n" if $verbose;

                                # --- Calculate overlap

                                if ($nextbeg < $endpos) {
                                    my $this = $endpos - $nextbeg;
                                    $totoverlap = $totoverlap + $this;
                                    print STDERR "overlap calculated is $totoverlap\n" if $verbose;
                                }

                                # --- If no overlap, stitch together

                                if ($sum == 0) {
                                    $sum = $lengtharray[$pos]+$lengtharray[$next];
                                } elsif ($sum > 0) {
                                    $sum = $sum +$lengtharray[$next];
                                }

                                print STDERR "sum of query is $sum\n" if $verbose;

                                my $contF = $comparecontigarray[$pos];
                                my $contS = $comparecontigarray[$next];

                                if (! exists $keepcontighashBeg{$contF} ) {
                                    push @keepcontigarray, $contF;
                                    $keepcontighashBeg{$contF} = $begarray[$pos];
                                    $keepcontighashEnd{$contF} = $endarray[$pos];
                                    $numtokeep = $numtokeep + 1;
                                }

                                if (! exists $keepcontighashBeg{$contS} ) {
                                    push @keepcontigarray, $contS;
                                    $keepcontighashBeg{$contS} = $begarray[$next];
                                    $keepcontighashEnd{$contS} = $endarray[$next];
                                    $numtokeep = $numtokeep + 1;
                                }

                                $pos = $next;
                                $next++;

                            # --- Overlaps, so find the longest

                            } elsif ($nextbeg <= ($endpos - $overlap)) {

                                if ($sum == 0) {

                                    my $firstsum = $lengtharray[$pos];
                                    my $secondsum = $lengtharray[$next];

                                    if ($firstsum >= $secondsum) {

                                        $sum = $firstsum;
                                        my $cont = $comparecontigarray[$pos];
                                        push @keepcontigarray, $cont;
                                        $keepcontighashBeg{$cont} = $begarray[$pos];
                                        $keepcontighashEnd{$cont} = $endarray[$pos];
                                        $numtokeep = $numtokeep+1;

                                    } elsif ($firstsum < $secondsum) {

                                        # Keep second contig
                                        $sum = $secondsum;
                                        my $cont = $comparecontigarray[$next];
                                        push @keepcontigarray, $cont;
                                        $keepcontighashBeg{$cont} = $begarray[$next];
                                        $keepcontighashEnd{$cont} = $endarray[$next];
                                        $numtokeep = $numtokeep+1;
                                        $pos = $next;
                                    }
                                    $next++;

                                } elsif ($sum != 0) {

                                    my $secondsum = $lengtharray[$next];

                                    if ($sum >= $secondsum) {

                                        $next++;
                                        for my $contig (@keepcontigarray) {}  

                                    } elsif ($sum < $secondsum) {

                                        $sum = $secondsum;
                                        my $cont = $comparecontigarray[$next];

                                        undef (@keepcontigarray);
                                        undef (%keepcontighashBeg);
                                        undef (%keepcontighashEnd);
                                        $totoverlap = 0;

                                        push @keepcontigarray, $cont;

                                        $keepcontighashEnd{$cont} = $endarray[$next];
                                        $keepcontighashBeg{$cont} = $begarray[$next];
                                        $numtokeep = 1;
                                        $pos = $next;
                                        $next++;
                                    }
                                }
                            }
                        }
                    }

                    $numtokeep = scalar(@keepcontigarray);
                    for (0..$numtokeep-1) {
                        $pos=$_;
                    }

                    # --- Do final print

                    print STATS "$numtokeep,$totoverlap,$sum,";

                    print STDERR "num to keep $numtokeep,total overlap $totoverlap, sum $sum," if $verbose;

                    for my $contig  (@keepcontigarray) {
                        print STATS "$keepcontighashBeg{$contig},$keepcontighashEnd{$contig},";
                        print "$keepcontighashBeg{$contig},$keepcontighashEnd{$contig}," if $verbose;
                    }

                    for my $contig  (@keepcontigarray) {
                        print STATS "$contig,";
                        print STDERR "$contig," if $verbose;
                    }
                }
            }

            print STATS "\n";
            print STDERR "\n" if $verbose;

            undef(@contigarray);
            undef(@bigarray);
        }
    }
}

exit;
