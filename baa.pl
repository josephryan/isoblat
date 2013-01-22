#!/usr/bin/perl

# baa.pl 
#
# BAA stands for use BLAT to ASSESS an ASSEMBLY
# Copyright (C) 2012,2013 Joseph F. Ryan 
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

$|++;
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

our $VERSION = 0.07;
our $AUTHOR  = 'Joseph F. Ryan <jfryan@yahoo.com>';

# run "baa.pl --help" for detailed info on these parameters
our $DEFAULT_MIN_TO_COUNT_AS_COVERAGE = 5;
our $DEFAULT_MAX_GAP_TO_CONSIDER_MISSING = 5; # set to 0 to turn this off

# if the bam is sorted it is unnecessary to read the whole
# blat file into memory and sort
# This is not tested, do not use unless you want to test it :-)
our $BAM_PRESORTED = 0; 

# suffix of file to put rearrangments
# it creates a new file that has your blat file name + $REARRANGEMENTS_SUF
# e.g. transcripts_v_assembly.blat.rearrangements
# set to 0 if you don't want to print the rearrangements
# or better, use --do_not_print_rearrangements flag on command line
our $REARRANGEMENTS_SUF = 'rearrangements';

MAIN: {
    my $opt_version = 0;
    my $opt_help = 0;
    my $max_gap_to_consider_missing = $DEFAULT_MIN_TO_COUNT_AS_COVERAGE;
    my $min_to_count_as_coverage = $DEFAULT_MIN_TO_COUNT_AS_COVERAGE; 
    my $do_not_print_rearrangements = 0;

    my $opt_results = Getopt::Long::GetOptions(  "version" => \$opt_version,
               "max_gap_to_consider_missing" => \$max_gap_to_consider_missing,
               "do_not_print_rearrangements" => \$do_not_print_rearrangements,
                  "min_to_count_as_coverage" => \$min_to_count_as_coverage,
                                      "help" => \$opt_help);

    $opt_version && version();
    pod2usage({-exitval => 0, -verbose => 2}) if $opt_help;
    my $blat   = $ARGV[0] or die usage();
    my $fasta  = $ARGV[1] or die usage();
    $REARRANGEMENTS_SUF = 0 if $do_not_print_rearrangements;
    print "# running version $VERSION of baa.pl\n";
    print "# run with this command: $0 @ARGV\n";
    print "# \$max_gap_to_consider_missing = $max_gap_to_consider_missing\n";
    print "# \$min_to_count_as_coverage    = $min_to_count_as_coverage\n";
    print "# \$BAM_PRESORTED = $BAM_PRESORTED\n" if ($BAM_PRESORTED);
    print "\n";

    my $fasta_count  = get_fasta_count($fasta);
    my $rh_dat = get_data($blat,$max_gap_to_consider_missing,$min_to_count_as_coverage);
    print_stats($rh_dat,$fasta_count,"$blat.$REARRANGEMENTS_SUF");
}


sub get_fasta_count {
    my $file = shift;
    my $count = 0;
    open IN, $file or die "cannot open $file:$!";
    while (my $line = <IN>) {
        $count++ if ($line =~ m/^>/);
    }
    return $count;
}

sub print_stats {
    my $rh_dat = shift;
    my $tot = shift;
    my $r_file = shift;

    my $hits = scalar(keys(%{$rh_dat}));
    my $perc = $hits / $tot;
    print "Percentage of transcripts with a BLAT entry ($hits/$tot): $perc\n";

    my %perc_cov_windows = ();
    my %counts = (); 
    my $tot_mapped = scalar(keys %{$rh_dat});
    my $tot_cov = 0;
    my $tot_len = 0;
    my $tot_ctgs = 0;
    my $potential_rearrangements = 0;
    open OUT, ">$r_file" or die "cannot open >$r_file:$!" if ($REARRANGEMENTS_SUF);
    foreach my $id (keys %{$rh_dat}) {
        my $coverage = $rh_dat->{$id}->{'covered'} / $rh_dat->{$id}->{'len'};
        $perc_cov_windows{'90'}++ if ($coverage >= 0.9);
        $perc_cov_windows{'80'}++ if ($coverage >= 0.8);
        $perc_cov_windows{'70'}++ if ($coverage >= 0.7);
        $perc_cov_windows{'less'}++ if ($coverage < 0.7);
        if ($rh_dat->{$id}->{'potential_rearrangements'}) {
            $potential_rearrangements += $rh_dat->{$id}->{'potential_rearrangements'};
            print OUT "$id\n" if ($REARRANGEMENTS_SUF);
        }
        $tot_cov += $rh_dat->{$id}->{'covered'};
        $tot_len += $rh_dat->{$id}->{'len'};
# Error Check
#print "\$id = $id\n";
#if ($tot_cov > $tot_len) {
#    print "$rh_dat->{$id}->{'covered'}\n";
#    print "$rh_dat->{$id}->{'len'}\n";
#    die "$id: $tot_cov > $tot_len";
#}

        my $num_ctgs = $rh_dat->{$id}->{'coverage_count'};
        $counts{'one_ctg'}++     if ($num_ctgs == 1);
        $counts{'two_ctg'}++     if ($num_ctgs == 2);
        $counts{'three_ctg'}++   if ($num_ctgs == 3);
        $counts{'four_ctg'}++    if ($num_ctgs == 4);
        $counts{'five_ctg'}++    if ($num_ctgs == 5);
        $counts{'p_five_ctg'}++  if (($num_ctgs > 5) && ($num_ctgs <= 10));
        $counts{'p_ten_ctg'}++   if ($num_ctgs > 10);
        $tot_ctgs += $num_ctgs;
    }

    my $tmp_p = $tot_cov / $tot_len;
    print "Total % coverage of all positions ($tot_cov / $tot_len): $tmp_p\n";

    $tmp_p = $tot_ctgs / $tot_mapped;
    print "Average number of contigs per mapped transcript: $tmp_p\n";

    print "Number of potential rearrangements = $potential_rearrangements\n\n";

}

sub get_data {
    my $file = shift;
    my $max_gap_to_consider_missing = shift;
    my $min_to_count_as_coverage = shift;

    my %data = ();
    open IN, $file or die "cannot open $file:$!";
    my $first_line = <IN>;
    unless ($first_line =~ m/^psLayout version 3\s*$/) {
        die "expecting version 3 got: $first_line"
    }
    for (0..3) { <IN> }; # remove first 4 lines
    my @current = ();
    my $curr_que = '';
    if ($BAM_PRESORTED) {
        while (my $line = <IN>) {
            chomp $line;
            my @fields = split/\t/, $line;
            if ($fields[9] ne $curr_que && $data{$fields[9]}) {
                die "blat file is not presorted: $fields[9]\n";
            }
            $curr_que = process_fields(\@fields,\@current,$curr_que,\%data,$max_gap_to_consider_missing,$min_to_count_as_coverage);
        }
    } else {
        my @all_lines = ();
        while (my $line = <IN>) {
            chomp $line;
            my @fields = split/\t/, $line;
            push @all_lines, \@fields;
        }
        foreach my $ra_f (sort {$a->[9] cmp $b->[9]} @all_lines) {
            $curr_que = process_fields($ra_f,\@current,$curr_que,\%data,$max_gap_to_consider_missing,$min_to_count_as_coverage);
        }
    }
    process_query(\@current,\%data,$max_gap_to_consider_missing,$min_to_count_as_coverage);
    return \%data;
}

sub process_fields {
    my $ra_fields = shift;
    my $ra_current = shift;
    my $curr_que = shift;
    my $rh_data = shift;
    my $max_gap_to_consider_missing = shift;
    my $min_to_count_as_coverage = shift;

    if (!$curr_que || $curr_que eq $ra_fields->[9]) {
        push @{$ra_current}, $ra_fields;
        $curr_que = $ra_fields->[9];
    } else {
        $curr_que = $ra_fields->[9];
        process_query($ra_current,$rh_data,$max_gap_to_consider_missing,$min_to_count_as_coverage);
        @{$ra_current} = ();
        push @{$ra_current}, $ra_fields;
    }
    return $curr_que;
}

sub process_query {
    my $ra_curr = shift;
    my $rh_dat  = shift;
    my $max_gap_to_consider_missing = shift;
    my $min_to_count_as_coverage = shift;

    my $que = $ra_curr->[0]->[9];
    $rh_dat->{$que}->{'len'} = $ra_curr->[0]->[10];
    $rh_dat->{$que}->{'potential_rearrangements'} = 0;
    my @covered = ();
    foreach my $ra_f (sort {($b->[12] - $b->[11]) <=> ($a->[12] - $a->[11])} @{$ra_curr}) {
        my $adds_to_coverage = 0;
        my @blk_sizes = split ',', $ra_f->[18]; # blockSizes
        my @q_st    = split ',', $ra_f->[19];   # queryStarts
        
        for (my $i = 0; $i < scalar(@blk_sizes); $i++) {
            for (my $j = $q_st[$i]; $j < ($q_st[$i] + $blk_sizes[$i]); $j++) {
                $adds_to_coverage++ unless ($covered[$j]);
                $covered[$j] = 1;
            }
        }
        next unless ($adds_to_coverage >= $min_to_count_as_coverage);
        $rh_dat->{$que}->{'coverage_count'}++;
        if ($rh_dat->{$que}->{'dup_check'}->{$ra_f->[13]}) {
            $rh_dat->{$que}->{'potential_rearrangements'}++;
        } else {
            push @{$rh_dat->{$que}->{'contigs'}}, $ra_f->[13];
        }
        $rh_dat->{$que}->{'dup_check'}->{$ra_f->[13]}++;
    }        
    fill_in_gaps(\@covered,$max_gap_to_consider_missing) if ($max_gap_to_consider_missing);
    foreach my $cov (@covered) {
        $rh_dat->{$que}->{'covered'}++ if ($cov);
    }
}

sub fill_in_gaps {
    my $ra_cov = shift;
    my $max_gap_to_consider_missing = shift;
    my $in_a_gap = 0;
    for (my $i = 0; $i < @{$ra_cov}; $i++) {
        if ($in_a_gap) {
            if ($ra_cov->[$i]) {
                if ($in_a_gap < $max_gap_to_consider_missing) {
                    for (my $j = 1; $j <= $in_a_gap; $j++) {
                        $ra_cov->[$i -$j] = 1;
                    }
                }
                $in_a_gap = 0;
            } else {
                $in_a_gap++;
            }
        } else {
            $in_a_gap = 1 unless ($ra_cov->[$i]);
        } 
    }
}

sub _set_length {
    my $rh_data = shift;
    my $id = shift;
    my $len = shift;
    if ($rh_data->{$id}) {
        unless ($rh_data->{$id}->{'len'} == $len) {
            die "unexpected transcript length: $len";
        }
    } else {
        $rh_data->{$id}->{'len'} = $len;
    }
}

sub version {
    die "baa.pl $VERSION\n";
}

sub usage {
    die "usage: $0 [--version] [--help] [--max_gap_to_consider_missing=INT] [--min_to_count_as_coverage=INT] [--do_not_print_rearrangements] BLAT_FILE FASTA_QUERY_USED_IN_BLAT\n";
}

__END__

=head1 NAME

B<baa.pl> - use Blat to Assess an Assembly

=head1 AUTHOR

Joseph F. Ryan <josephryan@yahoo.com>

=head1 SYNOPSIS

baa.pl [--version] [--help] [--max_gap_to_consider_missing=INT] [--min_to_count_as_coverage=INT] [--do_not_print_rearrangements] BLAT_FILE FASTA_QUERY_USED_IN_BLAT

=head1 BLAT_FILE

=over 2

This program depends on Blat. An alignment tool written by Jim Kent.
Blat is freely available for academic, nonprofit and personal use.  See:
http://genome.ucsc.edu/FAQ/FAQblat.html#blat3

Run Blat using your assembled transcripts/ests as a query and the assembly you
want to evaluate as the subject.  For example:

    blat assembly_1.fasta ests.fasta ests_v_assembly_1.blat
    blat assembly_2.fasta ests.fasta ests_v_assembly_2.blat

=back

=head1 FASTA_QUERY_USED_IN_BLAT

=over 2

this file is the fasta file used in the Blat search 
(in the 2 Blat examples above, 'ests.fasta' is the FASTA_QUERY_USED_IN_BLAT)

=back

=head1 OPTIONS

=over 2

=item B<--min_to_count_as_coverage>

(default: 5)
in the case where there is a match and less than this many positions
match they will not count towards the "Average number of contigs per
mapped transcript"

=item B<--default_max_gap_to_consider_missing>

(default: 5)
if there is a gap at the beginning, end, or in between pairs
only gaps less than this value will be considered missing
NOTE: this value does not effect Total "% coverage of all positions"

=item B<--do_not_print_rearrangements>
by default a file with the suffix '.rearrangements' is printed; this
file contains a list of transcript identifiers that have strange
alignments and may indicate mis-arrangements in the assembly

=item B<--help>

Print this manual
 
=item B<--version>

Print the version. Overrides all other options.

=back

=head1 DESCRIPTION

This script will parse the output of a Blat (transcriptome v. assembly)
and generate values which can be used to compare assemblies

=head1 POTENTIAL REARRANGEMENTS

Exons of a transcript are expected to align sequentially to a single
chromosomal region. Potential rearrangements are regions of an assembly
where the pattern of alignments of transcripts suggests that the
alignment is not sequential.

The current algorithm for estimating rearrangements is far from perfect
One problem is the sorting in the process_query subroutine Blat
hits are sorted by the coverage of the hit (col 12 - col 13 in blat)
This can cause the algorithm to evaluate the wrong hit and can subsequently
inflate the number of potential rearrangements (plus cause errors
in the estimate of coverage).  However, since the main purpose of the tool
is to compare assemblies, most of these false-positive rearrangements will
show up in both assemblies being compared.

This value is also susseptible to the "maxIntron" setting of Blat.  If 
Blat is run with the default introns longer than 750000 will be flagged
as a potential rearrangement. 

=head1 REARRANGEMENTS FILE

By default the program will print potential rearrangements to a file
(with the suffix '.rearrangements). This can be suppressed with the
--do_not_print_rearrangements flag. It is useful to compare the number
of common reported rearrangements between assemblies.  Also, it can be
useful to look closely at the distinct rearrangements in each assembly
to identify possible misassemblies.

Here is a recipe for counting common rearrangments from 2 '.rearrangements'
files that are the result of 2 baa.pl runs on 2 different assemblies

=over 2

sort ests_v_assembly_1.blat.rearrangements > 1.sorted
sort ests_v_assembly_2.blat.rearrangements > 2.sorted
join 1.sorted 2.sorted | wc -l

=back

=head1 BUGS

Please report them to <josephryan@yahoo.com>

=head1 COPYRIGHT

Copyright (C) 2012 Joseph F. Ryan 

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

=cut
