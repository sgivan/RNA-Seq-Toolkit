#!/bin/env perl

use 5.010;      # Require at least Perl version 5.8
use strict;     # Must declare all variables before using them
use warnings;   # Emit helpful warnings
use autodie;
use Getopt::Long; # use GetOptions function to for CL args
use lib '/ircf/ircfapps/share/biomart-perl/lib';
use BioMart::Initializer;
use BioMart::Query;
use BioMart::QueryRunner;

my ($debug,$verbose,$help);
my ($db,$cuffdiff_file,$annots);

my $result = GetOptions(
    "db=s"          =>  \$db,
    "cuff=s"        =>  \$cuffdiff_file,
    "annots=i"      =>  \$annots,
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
--cuff      cuffdiff output file

HELP
exit();
}

#my $confFile = "PATH TO YOUR REGISTRY FILE UNDER biomart-perl/conf/. For Biomart Central Registry navigate to http://www.biomart.org/biomart/martservice?type=registry";
my $confFile = "/ircf/ircfapps/share/biomart-perl/conf/registry.xml";
$db ||= 'phytozome';
$cuffdiff_file ||= 'gene_exp.diff';
$annots ||= 100;

# open cuffdiff output file
# should conform to format of *.diff files
open(DIFF, "<", $cuffdiff_file);
#
# NB: change action to 'clean' if you wish to start a fresh configuration  
# and to 'cached' if you want to skip configuration step on subsequent runs from the same registry
#
my $action='cached';
#my $action='update';
my $initializer = BioMart::Initializer->new('registryFile'=>$confFile, 'action'=>$action, mode => 'LAZYLOAD');
my $registry = $initializer->getRegistry;

my $query = BioMart::Query->new('registry'=>$registry,'virtualSchemaName'=>'default');
my $query_runner = BioMart::QueryRunner->new();

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
#$query_runner->printHeader();
$query_runner->uniqueRowsOnly(1);

my $cnt = 0;
while (<DIFF>) {
    next if (++$cnt == 1);
    say "cnt = $cnt" if ($verbose);
    last if ($debug && $cnt >= 10);
    my $line = $_;
    my @linevals = split /\t/, $line;
    my $id = $linevals[0];
        
    #$query->setDataset("phytozome");
    #$query->setDataset($db);
    #$query->addFilter("pac_transcript_id", ["23530342"]);
    #$query->addFilter("gene_name_filter", ["Phvulv091019806m.g"]);
    $query->addFilter("gene_name_filter", [$id]);
#    $query->addAttribute("gene_name1");
#    $query->addAttribute("transcript_name1");
#    $query->addAttribute("pfam_id");
#    $query->addAttribute("pfam_desc");
#    $query->addAttribute("smart_id");
#    $query->addAttribute("smart_desc");
#    $query->addAttribute("panther_id");
#    $query->addAttribute("panther_desc");
##    $query->addAttribute("kog_id");
##    $query->addAttribute("kog_desc");
#    $query->addAttribute("kegg_enzyme_id");
#    $query->addAttribute("kegg_enzyme_desc");
##    $query->addAttribute("ko_id");
##    $query->addAttribute("keggorth_desc");
#    $query->addAttribute("go_id");
#    $query->addAttribute("go_desc");
#
#    $query->formatter("TSV");
    #$query->limitSize(1);# this works

    #my $query_runner = BioMart::QueryRunner->new();
    ############################## GET COUNT ############################
    # $query->count(1);
    # $query_runner->execute($query);
    # print $query_runner->getCount();
    #####################################################################


    ############################## GET RESULTS ##########################
    # to obtain unique rows only
    #$query_runner->uniqueRowsOnly(1);

    $query_runner->execute($query);

    #print "\$query_runner isa '" . ref($query_runner) . "'\n";
    #print $query_runner->toString();
    #my $result_table = $query_runner->_getResultTable();
    #print ref($result_table) . "\n";

#    $query_runner->printHeader();
    $query_runner->printResults();
#    $query_runner->printFooter();
    #####################################################################

    #print "\n\nquery count: ", $query->count(),  "\n";

}
$query_runner->printFooter();

close(DIFF);

