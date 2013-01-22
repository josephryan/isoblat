baa.pl
======

use BLAT to ASSESS an ASSEMBLY - this program parses the output of a BLAT run of transcriptome vs. a genome

The output looks something like:

    Percentage of transcripts with a BLAT entry (27164/27315): 0.994471901885411
    Total % coverage of all positions (14238616 / 14485870): 0.982931366911342
    Average number of contigs per mapped transcript: 1.07377411279635
    Number of potential rearrangements = 122

INSTALLATION
------------

If your perl is not located in /usr/bin/perl 
you need to change the first line of baa.pl
for example to:

    #!/usr/local/bin/perl

You may need to make baa.pl executable

    chmod 755 baa.pl

Next copy baa.pl to a directory in your path
for example: /usr/local/bin:

    sudo cp baa.pl /usr/local/bin/

RUN
---

    baa.pl [--version] [--help] [--max_gap_to_consider_missing=INT] [--min_to_count_as_coverage=INT] [--do_not_print_rearrangements] BLAT_FILE FASTA_QUERY_USED_IN_BLAT

for detailed documenation

    baa.pl --help

DEPENDENCIES
------------

This module requires Perl

COPYRIGHT AND LICENCE
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

