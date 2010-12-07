#!/usr/bin/perl
#
use strict;
use Getopt::Long;
use Data::Dumper;
#use lib '/local/cluster/dev/lib/perl5/site_perl/5.8.8';
use Bio::SeqIO;
use Bio::Seq::Quality;

my ($read_1,$read_2,$debug,$help,$verbose,$read_1_out,$read_2_out,$maxN,$nomaxN,$smaller,$dump);

GetOptions(
            "read_1=s"        =>  \$read_1,
            "read_2=s"        =>  \$read_2,
            "read_1_out=s"    =>  \$read_1_out,
            "read_2_out=s"    =>  \$read_2_out,
            "debug"           =>  \$debug,
            "verbose"         =>  \$verbose,
            "help"            =>  \$help,
            "maxN=i"            =>  \$maxN, # crude quality filter; max allowable number of N's in read
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
$maxN = 4 unless ($maxN);
$verbose = 1 if ($debug);

if ($debug) {
    print "read_1 = $read_1\nread_2 = $read_2\nread_1_out = $read_1_out\nread_2_out = $read_2_out\nmaxN = $maxN\n";
}

my ($infile1_size,$infile2_size) = (0,0);
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
  $follower_matched_filename = $read_1_out;
} else {
  print "$read_2 has more sequences ($infile2_size > $infile1_size)\n" if ($debug);
  $masterfile = $read_1;
  $master_matched_filename = $read_1_out;
  $followfile = $read_2;
  $follower_matched_filename = $read_2_out;
}

print "using '$masterfile' as master (smaller) fastq file\n" if ($verbose);

#
# proceed with building content hash for smaller fastq file
#

my (%master,%follow,%temp,$lastID);
open(MASTER,$masterfile) or die "can't open '$masterfile': $!";

my ($linecnt,$seqcnt) = (0,0);
while (<MASTER>) {
  ++$linecnt;
  my $line;
  $line = $_;
  last if (++$seqcnt > 24 && $debug);
  chomp($line);
  print "$line\n" if ($debug);

#  if ($line =~ /([\@\+])(\S+)/) {
  if ($line =~ /([\@\+])(\S*)/) {
    my ($type,$id,$string);
    if ($2) {
        $id = $2
    } else {
        $id = '';
    }

    #
    # each sequence should have DNA string and quality string
    # sequence lines are preceded by lines starting with @
    # quality lines are preceded by lines starting with +
    #
    if ($1 eq '@') {
      $type = 'seq';
    } else {
      $type = 'qual';
    }

    #
    # if this is a new ID, save current %temp and start building new one
    #
    if ($id ne $lastID) {
    
      if ($debug) {
        print "\n\n\%temp:\n";
        print Dumper(\%temp), "\n\n";
      }

      # this should not save %temp unless both qual and seq are defined
      #
      $master{$lastID} = {%temp} if (defined($temp{qual}) && defined($temp{seq}));
      $lastID = $id;  # set $lastID to current ID
      %temp = ();     # clear %temp
    }

    #
    # collect next line of text from input fastq file
    # we should already know what to expect, based on $type determined above
    #
    $string = <MASTER>;
    chomp($string);
    print "string = '$string'\n" if ($debug);
    #
    # do some basic QC checks
    # in future, I should write seqs that fail to their own file
    #
    my $copy = $string;
    my $Ns = $copy =~ tr/N/N/;
    $maxN = (length($copy) + 1) if ($nomaxN);
    if ($type eq 'seq' && $Ns >= $maxN) { # only do this for DNA sequence lines
      print "$id has too many N's, skipping\n" if ($verbose);
      --$seqcnt;
      #
      # if we are skipping this sequence, must also skip the following quality lines
      #
      <MASTER>;
      <MASTER>;
      %temp = (); # and erase %temp to start new
      next;       # re-evaluate loop
    } elsif ($type eq 'seq') {
      print "$id is good, continuing\n" if ($debug);
    }
    #
    # save data to %temp
    #
    $temp{$type} = $string;
    
#    $master{$id}{$type} = $string;
  } else { # end of if ($line =~ /([\@\+])(\S+)/)
    print "don't know what to do with this line: '$line'\n";
  }
} # end of while (<MASTER>) loop
close(MASTER);
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

# print "'", $master{'HWI-EAS121:2:1:2:1808#0/2'}{'seq'}, "'\n";
# print "'", $master{'HWI-EAS121:2:1:2:1808#0/2'}{'qual'}, "'\n";
# 
# print "'", $master{'HWI-EAS121:2:1:11:363#0/2'}{'seq'}, "'\n";
# print "'", $master{'HWI-EAS121:2:1:11:363#0/2'}{'qual'}, "'\n";
# # 'CAGTGGTATCAACGCAGAGTACTTTTTTTTTTTTTTTTTTTTTTTTTNNTTTTTNTTTTTTTTTTTGGGGGGTGTTTGTT'
# # '\bb`aba`]_b``_PZba`G_Qaaaaab_Zabbbaabbbbba`aaa]DD]]QHFDMPKRPNTWZY[FFFNZLHLFMPMNH'
# 

#
# open follower input fastq file
# and output files that will contain:
#   first mate pair fastq
#   second mate pair fastq
#   mate pairs in follower file missing from master file
#
print "opening second fastq file\n" if ($verbose);

my $follow = Bio::SeqIO->new(
                                -file     =>  $followfile,
                                -format   =>  'fastq',
                              );

my $followout = Bio::SeqIO->new(
#                             -file      =>  ">$followfile" . ".out.fq",
#                             -file      =>  ">$followfile" . ".matched.fq",
                             -file      =>  ">$follower_matched_filename",
                             -format    =>  'fastq',
                          );

my $masterout = Bio::SeqIO->new(
#                              -file     =>  ">$masterfile" . ".out.fq",
#                              -file     =>  ">$masterfile" . ".matched.fq",
                              -file     =>  ">$master_matched_filename",
                              -format   =>  'fastq',
                            );

# my $master_nomate = Bio::SeqIO->new(
#                                     -file   =>  ">$masterfile" . ".nomate",
#                                     -format =>  'fastq',
#                                   );


my $follow_nomate = Bio::SeqIO->new(
                                -file     =>    ">$followfile" . ".nomate.fq",
                                -format   =>    'fastq',
                                );

my ($out1_cnt,$out2_cnt,$loop_cnt) = (0,0,0);
#
# Now, loop through follower file and determine which reads have mates
# in master file and print them to a file.
# Also, print reads in follower file that don't have mates to their own file
#
while (my $read2 = $follow->next_seq()) {
  last if (++$loop_cnt >= 2000 && $debug);
  my $seqid = $read2->id();
  my $seqid_tr = $seqid;
 
  my $seq = $read2->seq();
  #
  # Apply the same QC as for master seqs
  #
  my $Ns = $seq =~ tr/N/N/;
  $maxN = (length($seq) + 1) if ($nomaxN);
  next if ($Ns >= $maxN);

#
#   Sequence ID's are suffixed with either /1 or /2, depending on
#   which read they are. When querying %master, I need to convert
#   from one suffix to the other, because %master was created from
#   file containing the matching mate pair so it will have the
#   opposite suffix.
#
#  print "\$seqid_tr starts as '$seqid_tr'\n" if ($debug);
  if ($seqid_tr =~ /\/1/) {
    $seqid_tr =~ s/\/1/\/2/;
  } else {
    $seqid_tr =~ s/\/2/\/1/;
  }

#  print "checking for '$seqid_tr'\n" if ($debug);

  if ($master{$seqid_tr}) { # if a mate is in %master, print it to mate file
    $followout->write_fastq($read2);
    if ($debug) {
      my $quals = $read2->qual();
      print "writing to '$followfile':\n\@$seqid\n", $read2->seq(), "\n\+$seqid\n@$quals\n";
    }
    ++$master{$seqid_tr}->{mate}; # remember that this master sequence had a mate

    my @quals = split(//, $master{$seqid_tr}->{qual});
    print "quals : @quals\n" if ($debug);
    #
    # Illumina quality values are weird
    # translate them to standard values
    #
    @quals = map(ord() - 33, @quals);
    print "quals2: @quals\n" if ($debug);

    my $seq1 = Bio::Seq::Quality->new(
                                    -seq      =>  $master{$seqid_tr}->{seq},
 #                                   -qual     =>  $master{$seqid_tr}->{qual},
                                    -qual     =>  \@quals,
                                    -id       =>  $seqid_tr,
                                    );
    $masterout->write_fastq($seq1);

   last if (++$out2_cnt >= 5 && $debug);
  } else {  # else print it to the nomate file
    $follow_nomate->write_fastq($read2);
  }
} # end of while (my $read2 = $follow->next_seq())

#
# now, write_seq() for read_1's with no mate ...
#
# order doesn't matter here
#
open(NOMATE, ">$masterfile" . ".nomate.fq") or die "can't open $masterfile.nomate.fq: $!";
while (my($key,$value) = each(%master)) {
  if (!$value->{mate}) {
    print NOMATE "\@$key\n", $value->{seq}, "\n\+$key\n", $value->{qual}, "\n";
  
  }
}

close(NOMATE) or warn("can't close $masterfile.nomate.fq properly: $!");

sub _help {
print <<HELP;

            "read_1=s"        =>  \$read_1,
            "read_2=s"        =>  \$read_2,
            "read_1_out=s"    =>  \$read_1_out,# not yet implemented
            "read_2_out=s"    =>  \$read_2_out,# not yet implemented
            "debug"           =>  \$debug,
            "verbose"         =>  \$verbose,
            "help"            =>  \$help,
            "maxN"            =>  \$maxN, # crude quality filter; max allowable number of N's in read
            "smaller=s"       =>  \$smaller, # to specify which file is shorter (not yet implemented)
            "dump"            =>  \$dump, # print Data::Dumper output

HELP
exit(0);
}
