#!/usr/bin/perl
#
use 5.8.8;
use strict;
use Getopt::Long;
use Data::Dumper;
#use lib '/local/cluster/dev/lib/perl5/site_perl/5.8.8';
#use Bio::SeqIO;
#use Bio::Seq::Quality;

my ($read_1,$read_2,$debug,$help,$verbose,$read_1_out,$read_2_out,$maxN,$nomaxN,$smaller,$dump,$debug_max);

GetOptions(
            "read_1=s"        =>  \$read_1,
            "read_2=s"        =>  \$read_2,
            "read_1_out=s"    =>  \$read_1_out,
            "read_2_out=s"    =>  \$read_2_out,
            "debug"           =>  \$debug,
            "verbose"         =>  \$verbose,
            "help"            =>  \$help,
#            "maxN=i"            =>  \$maxN, # crude quality filter; max allowable number of N's in read
            "nomaxN"          =>    \$nomaxN, # don't use maxN filter
            "smaller=s"       =>  \$smaller, # to specify which file is shorter (not yet implemented)
            "dump"            =>  \$dump, # print Data::Dumper output
);

#
# determine which file has fewer reads:
#
_help() if ($help);
$read_1 = 'infile1' unless ($read_1);
$read_2 = 'infile2' unless ($read_2);
#$read_1_out = $read_1 . ".out.fq" unless ($read_1_out);
$read_1_out = $read_1 . ".matched.fq" unless ($read_1_out);
#$read_2_out = $read_2 . ".out.fq" unless ($read_2_out);
$read_2_out = $read_2 . ".matched.fq" unless ($read_2_out);
#$maxN = 4 unless ($maxN);
$verbose = 1 if ($debug);
$debug_max = 0;

if ($debug) {
#    print "read_1 = $read_1\nread_2 = $read_2\nread_1_out = $read_1_out\nread_2_out = $read_2_out\nmaxN = $maxN\n";
    print "read_1 = $read_1\nread_2 = $read_2\nread_1_out = $read_1_out\nread_2_out = $read_2_out\n\n";
    #$debug_max = 200000;
    $debug_max = 2e10;
}

my ($infile1_size,$infile2_size,$followsize) = (0,0,0);
print "determining relative size of input files\n" if ($verbose);
if (!$smaller) { # must determine fastq file contains fewer sequences if not explicitly passed via --smaller
  if (-e $read_1) {
    $infile1_size = `grep -c -E '@.+' $read_1`;
    chomp($infile1_size);
  } else {
    print "'$read_1' is missing\n";
    print "typical usage:\nillumina_pe_matchup.pl --read_1 <infile1> --read_2 <infile2>\n\n";
    exit();
  }
  
  if (-e $read_2) {
    $infile2_size = `grep -c -E '@.+' $read_2`;
    chomp($infile2_size);
  } else {
    print "'$read_2' is missing\n";
    exit();
  }
}


my ($masterfile,$followfile,$master_matched_filename,$follower_matched_filename) = ('','','','');

if ($smaller) {
  $masterfile = $smaller; # use the smaller fastq file to build master hash
  
  $read_1 eq $masterfile ? $followfile = $read_2 : $followfile = $read_1; # other file is the "follower"
  
  print "user-specified master file: '$masterfile'\n" if ($debug);

} elsif ($infile1_size > $infile2_size) {
  print "$read_1 has more sequences ($infile1_size > $infile2_size)\n" if ($debug);
  $masterfile = $read_2;
  $master_matched_filename = $read_2_out;
  $followfile = $read_1;
  $followsize = $infile1_size;
  $follower_matched_filename = $read_1_out;
} else {
  print "$read_2 has more sequences ($infile2_size > $infile1_size)\n" if ($debug);
  $masterfile = $read_1;
  $master_matched_filename = $read_1_out;
  $followfile = $read_2;
  $followsize = $infile2_size;
  $follower_matched_filename = $read_2_out;
}

print "using '$masterfile' as master (smaller) fastq file\n" if ($verbose);

#
# proceed with building content hash for smaller fastq file
#

my ($MASTER,%master,%follow,%temp,$lastID);
open($MASTER,$masterfile) or die "can't open '$masterfile': $!";

my ($linecnt,$seqcnt) = (0,0);

while (my $read1 = _fastq_in($MASTER)) {
  ++$linecnt;
  last if (++$seqcnt > $debug_max && $debug);
  my $id = $read1->[0];

  if ($id) {
    my ($type,$string);

        $master{$id} = {
                            seq     =>  $read1->[1],
                            qual    =>  $read1->[3],
                            };

  } else { # end of if ($line =~ /([\@\+])(\S+)/)
#    print "don't know what to do with this line: '$line'\n";
    last;
  }

} # end of while loop 
close($MASTER);
#
# now I should have a large %master of contents of smaller fastq file
#

if ($dump) {  # this could potentially be LOTS of text; ie gigabytes to terabytes

print "\n\n", "+" x 24, "DATA DUMP", "+" x 24, "\n\n";
$Data::Dumper::Indent = 3;
$Data::Dumper::Useqq = 1;
print Dumper(\%master);

print "-" x 50, "\n";

}

#
# open follower input fastq file
# and output files that will contain:
#   first mate pair fastq
#   second mate pair fastq
#   mate pairs in follower file missing from master file
#
print "opening second fastq file\n" if ($verbose);

my $FOLLOW;
open($FOLLOW,$followfile) or die "can't open $followfile: $!";

print "opening output fastq files\n" if ($verbose);
my $FOLLOWOUT;
open($FOLLOWOUT,">$follower_matched_filename") or die "can't open $follower_matched_filename: $!";

my $MASTEROUT;
open($MASTEROUT,">$master_matched_filename") or die "can't open $master_matched_filename: $!";

my $FOLLOW_NOMATE;
open($FOLLOW_NOMATE,">$followfile" . ".nomate.fq") or die "can't open $followfile" . ".nomate.fq: $!";

my ($out1_cnt,$out2_cnt,$loop_cnt) = (0,0,0);
#
# Now, loop through follower file and determine which reads have mates
# in master file and print them to a file.
# Also, print reads in follower file that don't have mates to their own file
#
print "stepping through follower file\n" if ($verbose);

while (my $read2 = _fastq_in($FOLLOW)) {
  last if (++$loop_cnt >= $debug_max && $debug);
  my $seqid = $read2->[0];
#  last if (!$seqid);
  if (!$seqid) {
    print "reached end of follower file\n" if ($debug);
    last;
  }
#  print "seqid: '$seqid'\n" if ($debug);
  my $seqid_tr = $seqid;
 

#
#   Sequence ID's are suffixed with either /1 or /2, depending on
#   which read they are. When querying %master, I need to convert
#   from one suffix to the other, because %master was created from
#   file containing the matching mate pair so it will have the
#   opposite suffix.
#

    my $digit = substr($seqid_tr,-1,1,"");
    if ($digit == 1) {
        $seqid_tr .= "2";
    } else {
        $seqid_tr .= "1";
    }


#  print "checking for '$seqid_tr'\n" if ($debug);

  if ($master{$seqid_tr}) { # if a mate is in %master, print it to mate file
    _fastq_out($read2->[0],$read2->[1],$read2->[3],$FOLLOWOUT);

#    if ($debug) {
#      my $quals = $read2->qual();
#      print "writing to '$followfile':\n\@$seqid\n", $read2->seq(), "\n\+$seqid\n@$quals\n";
#    }
    ++$master{$seqid_tr}->{mate}; # remember that this master sequence had a mate

    _fastq_out($seqid_tr,$master{$seqid_tr}->{seq},$master{$seqid_tr}->{qual},$MASTEROUT);


  } else {  # else print it to the nomate file
    _fastq_out($read2->[0],$read2->[1],$read2->[3],$FOLLOW_NOMATE);
  }
} # end of while (my $read2 = $follow->next_seq())

#
# now, write_seq() for read_1's with no mate ...
#
# order doesn't matter here
#
my $NOMATE;
open($NOMATE, ">$masterfile" . ".nomate.fq") or die "can't open $masterfile.nomate.fq: $!";
while (my($key,$value) = each(%master)) {
  if (!$value->{mate}) {
    _fastq_out($key,$value->{seq}, $value->{qual},$NOMATE);
  
  }
}

close($NOMATE) or warn("can't close $masterfile.nomate.fq properly: $!");
close($FOLLOW) or warn();
close($FOLLOWOUT) or warn();
close($FOLLOW_NOMATE) or warn();
close($MASTEROUT) or warn();

sub _fastq_in {
    my $fh = shift;
    my $id = <$fh>;
    my $seq = <$fh>;
    my $id2 = <$fh>;
    my $qual = <$fh>;

#   trim first character from ID's
    substr($id,0,1,"");
    substr($id2,0,1,"");

    my @seqdata = ($id,$seq,$id2,$qual);
#    print "@seqdata\n";
#    exit();
    chomp(@seqdata);
    return [@seqdata];
}

sub _fastq_out {
    my $id = shift;
#    my $hashref = shift;
    my $seqstring = shift;
    my $qualstring = shift;
    my $filehandle = shift;

    print $filehandle "\@$id\n$seqstring\n\+$id\n$qualstring\n";
    return 1;
}

sub _help {
print <<HELP;

            "read_1=s"        =>  \$read_1,
            "read_2=s"        =>  \$read_2,
            "read_1_out=s"    =>  \$read_1_out,# not yet implemented
            "read_2_out=s"    =>  \$read_2_out,# not yet implemented
            "debug"           =>  \$debug,
            "verbose"         =>  \$verbose,
            "help"            =>  \$help,
            "smaller=s"       =>  \$smaller, # to specify which file is shorter (not yet implemented)
            "dump"            =>  \$dump, # print Data::Dumper output

HELP
exit(0);
}
