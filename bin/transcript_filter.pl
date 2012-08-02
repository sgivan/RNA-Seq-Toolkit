#!/usr/bin/env perl
#
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
use strict;
use warnings;
use Getopt::Long;

my ($debug,$infile,$min_length,$max_length,$outfails,$failsfile,$cumulative,$help);
#
# Perl script to filter a collection of transripts specificed in a GTF file
# by a minimum or maximum length value
#
# --cumulative tracks sum of lengths of each exon for a transcript, whereas
# if --cumulative is not specified, then length is based on start and stop coords
# of first and last exons, respectively.
#
GetOptions(

    "file=s"            =>  \$infile,
    "min_length=i"      =>  \$min_length,
    "max_length=i"      =>  \$max_length,
    "outfails"          =>  \$outfails,
    "cumulative"        =>  \$cumulative,
    "debug"             =>  \$debug,
    "help"              =>  \$help,
);

$infile = "-" unless ($infile);
$min_length ||= 500;
$max_length ||= 1e6;

if ($help) {
    print "usage: transcript_filter.pl --file infile_name\n";
    print "options:\n \
    \t--file infile_name or - \
    \t--min_length integer \
    \t--max_length integer \
    \t--outfails \
    \t--cumulative \
    \t--debug\n\n";
    exit;
}

open(IN,$infile) or die "can't open $infile: $!";
if ($outfails) {
    $failsfile = 'outfails.gtf';
    open(FAILS,">$failsfile") or die "can't open 'outfails.txt': $!";
}

my ($new_transcript_id,$transcript_id,$start,$stop,@buffer,$cnt,$cumlength);
($start,$stop,$cnt,$cumlength) = (0,0,0,0);
while (<IN>) {
    my $line = $_;
    my @vals = split /\t/, $line;
    if ($vals[8] =~ /\Wtranscript_id\s\"(.+?)\"\;/) {
        $new_transcript_id = $1;
    }
    $transcript_id = $new_transcript_id if (!$cnt);
    
    if ($new_transcript_id ne $transcript_id) {
#        print @buffer; 
        evalout($start,$stop,$cumlength,@buffer);
        @buffer = ();
        $transcript_id = $new_transcript_id;
        ($start,$stop,$cumlength) = ($vals[3],$vals[4],0);
    }

    push(@buffer,$line);
    $start = $vals[3] if ($vals[3] < $start || !$start);
    $stop = $vals[4] if ($vals[4] > $stop);
    $cumlength += $vals[4] - $vals[3] + 1;
#    print "loop '$cnt', start: '$start', vals[3]: '" . $vals[3] . "', stop: '$stop', vals[4]: '" . $vals[4] . "'\n";


} continue {
    ++$cnt;
    if (eof) {
#        print @buffer;
        evalout($start,$stop,$cumlength,@buffer);
    }
}

close(IN) or warn("can't close $infile properly: $!");
close(FAILS) or warn("can't close $failsfile properl: $!") if ($outfails);

sub evalout {
    my ($start,$stop,$cumlength,@bufferout) = @_;
    my $diff = $stop - $start;

#    if ($diff >= $min_length && $diff <= $max_length) {
    if (($cumulative && ($cumlength >= $min_length && $cumlength <= $max_length)) || (!$cumulative && ($diff >= $min_length && $diff <= $max_length))) {
        print "\nstart: $start; stop: $stop; cumulative: $cumlength\n" if ($debug);
        print @bufferout;
    } elsif ($outfails) {
        print FAILS @bufferout;
    }
}
