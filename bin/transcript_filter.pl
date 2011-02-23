#!/usr/bin/perl
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
