baa.pl
======

use BLAT to ASSESS an ASSEMBLY - this program parses the output of a BLAT run of transcriptome vs. a genome

The output looks something like:

    Percentage of transcripts with a BLAT entry (15661/15752): 0.994222955815135
    Total % coverage of all positions (10303278 / 10491059): 0.982100853688841
    Number of transcripts mapping to a single contig/scaffold: (14750/15661) 0.941830023625567
    Average number of contigs/scaffolds per mapped transcript: 1.08715918523721
    Number of potential rearrangements = 153

INSTALLATION
------------

To install this script and documentation type the following:

    perl Makefile.PL
    make
    make install

RUN
---

    baa.pl [--version] [--help] [--max_gap_to_consider_missing=INT] [--min_to_count_as_coverage=INT] [--do_not_print_rearrangements] BLAT_FILE FASTA_QUERY_USED_IN_BLAT

for detailed documenation

    perldoc baa.pl

DEPENDENCIES
------------

This module requires Perl

You will also want BLAT:

    http://users.soe.ucsc.edu/~kent/src/

COPYRIGHT AND LICENSE
------------

Copyright (C) 2012,2013 Joseph F. Ryan

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program in the file gpl.txt.  If not, see
http://www.gnu.org/licenses/.

