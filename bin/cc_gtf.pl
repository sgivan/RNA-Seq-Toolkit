#!/usr/bin/perl
#
use strict;
use warnings qw/all/;
use Getopt::Long;

#
# This script takes several output files from cuffcompare and generates
# a GTF file that can subsequently be used as a reference file for
# other analyses
#

my ($trackingfile,$inputgtf,$verbose,$debug,$help,$mincoverage,$Rfile);

GetOptions (
    "tracking=s"        =>  \$trackingfile,
    "gtf=s"             =>  \$inputgtf,
    "mincov=f"          =>  \$mincoverage,
    "Rfile"             =>  \$Rfile,
    "verbose"           =>  \$verbose,
    "debug"             =>  \$debug,
    "help"              =>  \$help,
);

if ($help)
{

print <<HELP;

--tracking      input tracking file
--gtf           input GTF file
--mincov        minimum coverage value to accept
--Rfile         generate tab-delimited statistics per transcript
--verbose
--debug
--help

HELP
exit(0);
}

$trackingfile ||= 'all.tracking';
$inputgtf ||= 'all.combined.gtf';
$verbose = 1 if ($debug);

my %coverage = ();

print "tracking file = '$trackingfile'\ninput GTF = '$inputgtf'\n" if ($debug);

open (TRACK,$trackingfile) or die "can't open $trackingfile: $!";

while (<TRACK>) {
    #print "$_" if ($debug);
    #
    # split line into tab-delimited fields
    #
    # if no reference GTF file was provided, column number is different
    #
    my @vals = split /\s+/, $_;
    my ($tid,$lid,$rgid,$rtid,$cc,@sdata);
    if (@vals == 6) {
        ($tid,$lid,$rgid,$rtid,$cc,@sdata) = @vals;
    } else {
        ($tid,$lid,$rgid,$cc,@sdata) = @vals;
        $rtid = '';
    }

    #
    # now parse sdata
    #
    # my initial goal is to extract the maximum coverage value
    #

    my $maxcov = 0;
    foreach my $dstring (@sdata) {
        next if (!$dstring || $dstring eq '-');
        my @splitvals = split /\|/, $dstring;
        $maxcov = $splitvals[6] if ($splitvals[6] > $maxcov);
    }
    $coverage{$tid} = $maxcov;
}

close(TRACK) or warn "can't close $trackingfile properl: $!";;

open (GTF,$inputgtf) or die "can't open $inputgtf: $!";

open(OUT,">out.gtf") or die "can't open out.gtf: $!";
if ($Rfile) {
    open(R,">Rout.txt") or die "can't open Rout.txt: $!";
    print R "TCONS\texons\tlength\tcoverage\n";
}

my ($tid,$start,$stop,@buff,$cnt,$ltid,$lattrs);

while (<GTF>) {
    #print "line:\t", $_ if ($debug);
    my @lvals = split /\t/, $_;
    if ($lvals[8] =~ /transcript_id\s\"(.+?)\"/) {
        $tid = $1;
    } else {
        print "can't parse transcript_id from '$lvals[8]'\n";
    }

    if ($ltid && $tid ne $ltid) {
        pout(\@buff,$start,$stop,$ltid);
        @buff = ();
        $start = $lvals[3];
        $stop = $lvals[4];
    }
    push(@buff,[@lvals]);
    $ltid = $tid;
    $start = $lvals[3] if (!$start || ($lvals[3] < $start));
    $stop = $lvals[4] if (!$stop || ($lvals[4] > $stop));

#    exit if (++$cnt == 10);
} continue {
   if (eof) {
        pout(\@buff,$start,$stop,$ltid);
   } 
}

close(GTF) or warn "can't close $inputgtf properly: $!";
close(OUT) or warn "can't close out.gtf properly: $!";
close(R) or warn "can't close Rout.txt: $!" if ($Rfile);

sub pout {
    my ($buff,$start,$stop,$tid) = @_;

    my @transcript = @{$buff[0]};
    my $exoncnt = scalar(@$buff);
    my $length = $stop - $start + 1;
    #
    # implement a coverage filter
    #
    return if ($mincoverage && $coverage{$tid} < $mincoverage);

    chomp($transcript[8]);
    $transcript[8] =~ s/\sexon_number.+?;//;
    $transcript[2] = 'transcript';
    $transcript[3] = $start;
    $transcript[4] = $stop;
    $transcript[8] .= " cov \"$coverage{$tid}\"; exons \"$exoncnt\"; length \"$length\";";
    $transcript[8] .= "\n";


    if ($Rfile) {
        print R "$tid\t$exoncnt\t$length\t$coverage{$tid}\n";
    }

    print "\n" if ($debug);
    unshift(@buff,[@transcript]);
    foreach my $lref (@buff) {
        print join "\t", ("out:",@$lref) if ($debug);
        print OUT join "\t", (@$lref);
    }
    print "\n" if ($debug);

}

