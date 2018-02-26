#!/usr/bin/env perl
use v5.10;
use strict;
use warnings;
use autodie;

use Test::More;
use File::Slurp qw(write_file);

use Data::Show;

my $DEBUG = shift // 0;

my @direction = qw(R1 R2);
my @sample_names = ( 'C1' .. 'C5', 'D6' .. 'D9', 'D10', 'E11' .. 'E15');
my @sample_names_with_repeated_controls =
    ( 'C1' .. 'C5', 'D6' .. 'D9', 'D10', 'C1' .. 'C5', 'E11' .. 'E15');

system('./reset_test') unless $DEBUG;

my $fastq_dir = 'fastq.dir';
mkdir $fastq_dir;

for my $sample_name (@sample_names) {
    state $sample_id = 1;
    for my $direction (@direction) {
        system("touch $fastq_dir/${sample_name}_${direction}_001.fastq");
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
    "FASTQ_DIR" : "$fastq_dir",
    "SAMPLES": {
        "CONTROL": [ "C1", "C2", "C3", "C4", "C5" ],
        "EXPERIMENTS": { "X": [  "D6",  "D7",  "D8",  "D9", "D10" ],
                         "Y": [ "E11", "E12", "E13", "E14", "E15" ]
        }
    },
    "REFERENCE" : {
        "GTF": "$reference_dir/$species.gtf",
        "FA":  "$reference_dir/$species.fa"
    },
    "FILTER": {
        "FA": "reference/filter.fa"
    }
}
END
    
write_file('config.json', $JSON_text);

# Make index directory (for now) 
system('mkdir index');

# let file system get caught up
sleep 2;

# WARNING: Don't do this without an index directory (for now)
system('./setup_exp_dirs');

my @result_files = `ls experiment_[XY]/*/* experiment_[XY]/*.gtf`;
my @result_dirs  = `ls -d experiment_[XY]/*index*`;

my @result = grep {$_} sort (@result_files, @result_dirs);
chomp @result;

show @result if $DEBUG;

is_deeply(\@result, expected(), 'Reassembled files in the correct directory structure');

system('./reset_test') unless $DEBUG;

done_testing;

sub expected {

    my @expected;
    my @generic_dirs = qw( s_1 s_2 s_3 s_4 s_5 s_6 s_7 s_8 s_9 s_10);
    my @links = ('set1.fq', 'set2.fq');
    
    my $index_sample_name=0;
    
    for my $exp ('X', 'Y') {

        my $exp_dir = "experiment_$exp";

        push @expected, "$exp_dir/transcripts.gtf";
        push @expected, "$exp_dir/index";
        push @expected, "$exp_dir/hisat_index";

        #TODO: Index files
        for my $generic_dir (@generic_dirs) { 
           for my $link (@links) {
                push @expected, "$exp_dir/$generic_dir/$link";
           }
           my $sample_name
                = $sample_names_with_repeated_controls[$index_sample_name];
           for my $direction (@direction) {
               push @expected, "$exp_dir/$generic_dir/${sample_name}_${direction}_001.fastq";
           }
           $index_sample_name++;
        }
    }

    @expected = sort @expected;
    show @expected if $DEBUG;
    return [@expected];
}
