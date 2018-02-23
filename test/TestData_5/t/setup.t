#!/usr/bin/env perl
use v5.10;
use strict;
use warnings;
use autodie;

use Test::More;
use File::Slurp qw(write_file);

use Data::Show;

my @direction = qw(R1 R2);
my @sample_names = ( 'C1' .. 'C5', 'D6' .. 'D9', 'D10', 'E11' .. 'E15');
my @sample_names_with_repeated_controls =
    ( 'C1' .. 'C5', 'D6' .. 'D9', 'D10', 'C1' .. 'C5', 'E11' .. 'E15');

system('./reset_test');

{ # CREATE TEMP FILES FOR TESTING
    my $fastq_dir = 'fastq.dir';
    mkdir $fastq_dir;
    
    for my $sample_name (@sample_names) {
        state $sample_id = 1;
        for my $direction (@direction) {
            system("touch $fastq_dir/${sample_name}_${direction}_001.fastq.fz");
        }
        $sample_id++;
    }
    
    my $reference_dir = 'reference.dir';
    
    mkdir $reference_dir;
    
    my $species = 'Genus_species';
    
    system("touch $reference_dir/$species.fa"); 
    system("touch $reference_dir/$species.gtf"); 
    
    my $JSON_text = <<"END";
{
    "SAMPLE_DIR" : "$fastq_dir",
    "SAMPLES": {
        "CONTROL": [ "C1", "C2", "C3", "C4", "C5" ],
        "EXPERIMENTS": {    "drop": [ "D6", "D7", "D8", "D9", "D10" ],
                         "running": [ "E11", "E12", "E13", "E14", "E15" ]
        }
    },
    "REFERENCE" : {
        "GTF": "$reference_dir/$species.gtf",
        "FA":  "$reference_dir/$species..fa"
    },
    "FILTER": {
        "FA": "reference/filter.fa"
    }
}
END
    
    write_file('config.json', $JSON_text);
}


my @result = `ls [XY]/*`;

my @result = sort @result;

is(@result, expected(), 'Reassembled files in the correct directory structure');

# system("rm -rf fastq.dir");
# system("rm -rf reference.dir");
# system("rm config.json");

sub expected {

    my @expected;
    my @generic_dirs = ('s1' .. 's9', 's10');
    my @links = ('set1.fq', 'set2.fq');
    
    my $index_sample_name=0;
    
    for my $exp ('X', 'Y') {
        push @expected, "$exp/transcripts.gtf";
        #TODO: Index files
        for my $generic_dir (@generic_dirs) { 
           for my $link (@links) {
                push @expected, "$exp/$generic_dir/$link";
           }
           my $sample_name
                = $sample_names_with_repeated_controls[$index_sample_name];
           for my $direction (@direction) {
               push @expected, "$exp/$generic_dir/${sample_name}_${direction}_001.fastq";
           }
           $index_sample_name++;
        }
    }

    @expected = sort @expected;
    return @expected;
}
