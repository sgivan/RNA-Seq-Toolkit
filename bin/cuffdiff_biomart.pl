#!/bin/env perl

# copyright Scott Givan, The University of Missouri, July 6, 2012
#
#    This file is part of the RNA-seq Toolkit, or RST.
#
#    RST is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    RST is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with RST.  If not, see <http://www.gnu.org/licenses/>.
#
use 5.010;      # Require at least Perl version 5.8
use strict;     # Must declare all variables before using them
use warnings;   # Emit helpful warnings
use autodie;
use Getopt::Long; # use GetOptions function to for CL args
use lib '/ircf/ircfapps/share/biomart-perl/lib';# <-- should point to directory of BioMart perl modules
use BioMart::Initializer;
use BioMart::Query;
use BioMart::QueryRunner;

my ($debug,$verbose,$help);
my ($db,$cuffdiff_file,$annots,$fast,$outfile,$tempfile);

my $result = GetOptions(
    "db=s"          =>  \$db,
    "cuff=s"        =>  \$cuffdiff_file,
    "annots=i"      =>  \$annots,
    "fast"          =>  \$fast,
    "outfile=s"     =>  \$outfile,
    "tempfile=s"    =>  \$tempfile,
    "debug"         =>  \$debug,
    "help"          =>  \$help,
    "verbose"       =>  \$verbose,
);

if ($help) {
    help();
    exit(0);
}

sub help {

say <<HELP;
Script to retrieve annotation from BioMart

Option      Description
--db        database to retrieve annotation; ie, phytozome
--cuff      cuffdiff output file [default = gene_exp.diff]
--annots    number of annotations to retrieve (and print) for each gene/transcript
--fast      use fast algorithm
            (you should use this if gene/transcript ID's are duplicated in your data file)
--outfile   name of output file containing original data
            plus annotation data appended on same line
--tempfile  file containing annotation data

HELP
exit();
}

#my $confFile = "PATH TO YOUR REGISTRY FILE UNDER biomart-perl/conf/. For Biomart Central Registry navigate to http://www.biomart.org/biomart/martservice?type=registry";
my $confFile = "/ircf/ircfapps/share/biomart-perl/conf/registry.xml";# <-- should point to the location of registry.xml
$db ||= 'phytozome';
$cuffdiff_file ||= 'gene_exp.diff';
$annots ||= 100;
$outfile ||= 'biomart_outfile.txt';
$tempfile ||= 'biomart_tempfile.txt';

# open cuffdiff output file
# should conform to format of *.diff files
open(my $DIFF, "<", $cuffdiff_file);
open(my $TMP, ">", $tempfile);
open(my $OUT, ">", $outfile);
#
# NB: change action to 'clean' if you wish to start a fresh configuration  
# and to 'cached' if you want to skip configuration step on subsequent runs from the same registry
#
my $action='cached';
#my $action='update';
my $initializer = BioMart::Initializer->new('registryFile'=>$confFile, 'action'=>$action, mode => 'LAZYLOAD');
my $registry = $initializer->getRegistry;

my $query = BioMart::Query->new('registry'=>$registry,'virtualSchemaName'=>'default');
#my $query_runner = BioMart::QueryRunner->new();

$query->setDataset($db);
$query->addAttribute("gene_name1");
$query->addAttribute("transcript_name1");
$query->addAttribute("pfam_id");
$query->addAttribute("pfam_desc");
$query->addAttribute("smart_id");
$query->addAttribute("smart_desc");
$query->addAttribute("panther_id");
$query->addAttribute("panther_desc");
#    $query->addAttribute("kog_id");
#    $query->addAttribute("kog_desc");
$query->addAttribute("kegg_enzyme_id");
$query->addAttribute("kegg_enzyme_desc");
#    $query->addAttribute("ko_id");
#    $query->addAttribute("keggorth_desc");
$query->addAttribute("go_id");
$query->addAttribute("go_desc");
$query->limitSize($annots);# this works
$query->formatter("TSV");

my $mquery_runner = BioMart::QueryRunner->new();
if (!$fast) {
    $mquery_runner->execute($query);
    $mquery_runner->printHeader($TMP);
    $mquery_runner->uniqueRowsOnly(1);
}

my %buff = ();
my @idlist = ();

my $cnt = 0;
while (<$DIFF>) {
    next if (++$cnt == 1);
#    say "cnt = $cnt" if ($verbose);
    last if ($debug && $cnt >= 10);
    my $line = $_;
    my @linevals = split /\t/, $line;
    my $id = $linevals[0];
    chomp($id);
    push(@idlist,$id);

    next unless (!exists($buff{$id}));
    my $query_runner;
    if ($fast) {
        $query_runner = BioMart::QueryRunner->new();
        $query_runner->uniqueRowsOnly(1);
    }
        
    $query->addFilter("gene_name_filter", [$id]);

    my $rtn = -99;
    if ($fast) {
        $rtn = $query_runner->execute($query);
        $buff{$id} = $query_runner;
    } else {
        $rtn = $mquery_runner->execute($query);
        $mquery_runner->printResults($TMP);
        if ($verbose) {
            print "$id\t";
            $mquery_runner->printCompletionStamp();
        }
    }

#    $query_runner->printHeader();
#    $query_runner->printResults();
#    $query_runner->printFooter();
    #####################################################################

    #print "\n\nquery count: ", $query->count(),  "\n";
    sleep(1);
}
$mquery_runner->printFooter() unless ($fast);

close($DIFF);

if ($fast) {
    my $fcnt;
    for my $id (@idlist) {
        $buff{$id}->printHeader($TMP) if (++$fcnt == 1);
        $buff{$id}->printResults($TMP);
        $buff{$id}->printFooter($TMP) if ($fcnt == 1);
    }
}

system("paste $cuffdiff_file $tempfile > $outfile");

