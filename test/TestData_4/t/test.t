#!/usr/bin/env perl
use v5.10;
use strict;
use warnings;
use autodie;

use Test::More;
use File::Slurp qw(slurp);

my $vers  = shift // '2.0.4';

my $test_log = 'test.log';

warn "Running RNA-Seq-Toolkit tests (may run for a long time)\n"; #Ends in newline so warning omits line number

system("module load HISAT2-$vers; t/versionless_setup_and_test.sh |& tee $test_log"); 

for my $file ( qw(
                     de_data.txt
                     gene_de.txt
                     transcript_de.txt
                     transcripts.gtf
                   )
)
{
    
    my $result   = slurp $file;
    my $expected = slurp "t/expected_hisat2-$vers/$file";
    
    is $result, $expected, "$file is as expected";
}

done_testing();
