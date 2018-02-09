#!/usr/bin/env perl
use v5.10;
use strict;
use warnings;
use autodie;

use Test::More;

for my $file ( qw(
                     de_data.txt
                     gene_de.txt
                     transcript_de.txt
                     transcripts.gtf
                   )
)
{
    
    my $result   = `diff $file t/expected_hisat2-2.0.4/$file`; 
    chomp $result;
    my $expected = '';
    
    is $result, $expected, "$file is as expected";
}

done_testing();
