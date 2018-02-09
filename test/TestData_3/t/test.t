#!/usr/bin/env perl
use v5.10;
use strict;
use warnings;
use autodie;

use Test::More;
use File::Slurp qw(slurp);

my $vers  = shift // '2.0.4';
my $DEBUG = shift // 1;

my $test_log = 'test.log';

warn "Running RNA-Seq-Toolkit tests (this may run for ten minutes or so)\n"; #Ends in newline so warning omits line number

if ($DEBUG) {
    system("module load HISAT2-$vers; t/versionless_setup_and_test.sh |& tee $test_log"); 
}
else {
    system("t/versionless_setup_and_test.sh &> /dev/null"); 
}

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

if ($DEBUG) {
    say "To reset the test files, please run 'reset_test'";
}
else {
    system "./reset_test";
}

done_testing();
