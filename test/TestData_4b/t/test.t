#!/usr/bin/env perl
use v5.10;
use strict;
use warnings;
use autodie;

use Test::More;
use File::Slurp qw(slurp);

for my $file ( qw(
                     gene_de.txt
                     de_gene_data.txt
                     de_transcript_data.txt
                     transcript_de.txt
                     transcripts.gtf
                     Chr19.gtf
                   )
)
{
    
    my $result   = slurp $file;
    my $expected = slurp "t/expected/$file";
    
    is $result, $expected, "$file is as expected";
}

done_testing();