#!/usr/bin/env perl 
#===============================================================================
#
#         FILE:  ballgown_setup.pl
#
#        USAGE:  ./ballgown_setup.pl  
#
#  DESCRIPTION:  Script to setup directory for running ballgown.Rscript
#
#      OPTIONS:  ---
# REQUIREMENTS:  ---
#         BUGS:  ---
#        NOTES:  ---
#       AUTHOR:  Dr. Scott A. Givan (sag), givans@missouri.edu
#      COMPANY:  University of Missouri, USA
#      VERSION:  1.0
#      CREATED:  10/21/16 12:06:51
#     REVISION:  ---
#===============================================================================

use 5.010;      # Require at least Perl version 5.10
use autodie;
use Getopt::Long; # use GetOptions function to for CL args
use warnings;
use strict;

my ($debug,$verbose,$help,$controldir,$expdir);

my $result = GetOptions(
    "cont:s" =>  \$controldir,
    "exp:s"     =>  \$expdir,
    "debug"     =>  \$debug,
    "verbose"   =>  \$verbose,
    "help"      =>  \$help,
);

if ($help) {
    help();
    exit(0);
}

sub help {

    say <<HELP;


HELP

}

open(PHENO,">","pheno.txt");
say PHENO "ids,replicate,sample";

my $cnt = 1;
my $repcnt = 1;
for my $dir (`ls -d $controldir`) {
    chomp($dir);
    say "'" . $dir . "'";
    my $label = "sample";
    if ($cnt < 10) {
        $label .= "0$cnt";
    } else {
        $label .= $cnt;
    }
    symlink($dir,$label);
    say PHENO "$label,$repcnt,control";
    ++$cnt;
    ++$repcnt;
}

#say "identifying experimental sample directories '$controldir'";

$repcnt = 1;
for my $dir (`ls -d $expdir`) {
    chomp($dir);
    say "'" . $dir . "'";
    my $label = "sample";
    if ($cnt < 10) {
        $label .= "0$cnt";
    } else {
        $label .= $cnt;
    }
    symlink($dir,$label);
    say PHENO "$label,$repcnt,experimental";
    ++$cnt;
    ++$repcnt;
}

close(PHENO);

