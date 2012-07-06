#!/usr/bin/env perl
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
#
# Copyright 2012 Scott A. Givan
#
use strict;
use warnings qw/all/;
use Getopt::Long;

#
# This script takes several output files from cuffcompare and generates
# a GTF file that can subsequently be used as a reference file for
# other analyses. It also will filter by minimum coverage value.
#

my ($trackingfile,$inputgtf,$verbose,$debug,$help,$mincoverage,$Rfile,$clean,%clean);

GetOptions (
    "tracking=s"        =>  \$trackingfile,
    "gtf=s"             =>  \$inputgtf,
    "mincov=f"          =>  \$mincoverage,
    "Rfile"             =>  \$Rfile,
    "clean"             =>  \$clean,
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
--clean         remove transcripts with class codes i,p,x,s
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
    print "\@sdata = '@sdata'\n" if ($debug);

    #
    # we may have already parsed this and cleaned it out
    #if ($clean) {
    #    next if ($clean{$lid});
    #}
    #
    # now parse sdata
    #
    # my initial goal is to extract the maximum coverage value
    #

    my $maxcov = 0.0;
    foreach my $dstring (@sdata) {
        print "\$dstring = '$dstring'\n" if ($debug);
        next if (!$dstring || $dstring eq '-');
        my @splitvals = split /\|/, $dstring;
        print "\$splitvals[6] = '", $splitvals[6], "'\n" if ($debug);
        $maxcov = $splitvals[6] if ($splitvals[6] && $splitvals[6] > $maxcov);
    }
    $coverage{$tid} = $maxcov;

    if ($clean) {
        $clean{$lid} = 1 if ($cc eq 'i' || $cc eq 'p' || $cc eq 'x' || $cc eq 's');
    }
}

close(TRACK) or warn "can't close $trackingfile properl: $!";;

open (GTF,$inputgtf) or die "can't open $inputgtf: $!";

open(OUT,">out.gtf") or die "can't open out.gtf: $!";
if ($Rfile) {
    open(R,">Rout.txt") or die "can't open Rout.txt: $!";
    print R "TCONS\texons\tlength\tcoverage\n";
}

my ($lgid,$tid,$start,$stop,@buff,$cnt,$ltid,$lattrs);

while (<GTF>) {
    #print "line:\t", $_ if ($debug);
    my @lvals = split /\t/, $_;

    if ($lvals[8] =~ /gene_id\s\"(.+?)\"/) {
        $lgid = $1;
    } else {
        print "can't parse gene_id from '$lvals[8]'\n";
    }
    if ($clean) {
        next if ($clean{$lgid});
    }

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

