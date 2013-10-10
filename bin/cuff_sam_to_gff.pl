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
use strict;
use Getopt::Long;
use Bio::DB::Sam;
use Bio::DB::Sam::Constants; # this will import RFLAGS into this namespace
use File::Temp; # for the sort command

my ($infile,$outfile,$debug,$source,$type,$help,$integers,$mapped_id_file,$onlymapped,$sort_by_refmol,$temp_directory);
my $split = 0;


GetOptions(
            "infile=s"        =>  \$infile,
            "outfile=s"       =>  \$outfile,
            "source=s"        =>  \$source,
            "type=s"          =>  \$type,
            "debug"           =>  \$debug,
            "split"           =>  \$split,
            "integers"        =>  \$integers,
            "mapped"          =>  \$mapped_id_file,
            "onlymapped"      =>  \$onlymapped,
            "sort"            =>  \$sort_by_refmol,
            "tempdir=s"       =>  \$temp_directory,
            "help"            =>  \$help,
);
_help() if ($help);
$outfile = 'cuff_sam_to_gff.txt' unless ($outfile);
$infile = 'infile' unless ($infile);
$type = 'match' unless ($type);
#$source = 'cufflinks' unless ($source);
$source = 'QDread' unless ($source);
$temp_directory = '/tmp' unless ($temp_directory);

if ($debug) {
    print "outfile: '$outfile'\ninfile: '$infile'\ntype: '$type'\nsource: '$source'\ntempdir: '$temp_directory'\n";
    #exit();
}

my %GFFLOC = (

  refmol  =>  0,
  source  =>  1,
  type    =>  2,
  start   =>  3,
  stop    =>  4,
  score   =>  5,
  strand  =>  6,
  phase   =>  7,
  group   =>  8,
  
);

# my $RFLAGS = RFLAGS;
# my @flags = keys %$RFLAGS;
# print "flags: @flags\n";

# if (0x0010 & RFLAGS->{REVERSED}) {
#   print "match\n";
# } else {
#   print "don't match\n";
# }
# exit();

my $IN;
if ($infile eq '-') {
    $IN = *STDIN;
} else {
    open($IN,"<",$infile) or die "can't open $infile: $!";
}
open(OUT,">",$outfile) or die "can't open $outfile: $!";

if ($mapped_id_file) {
    open(IDFILE,">","idfile.txt") or die "can't open 'idfile.txt': $!";
}

my $cnt = 0;
#foreach my $line (<IN>) {
foreach my $line (<$IN>) {
	#last if (++$cnt > 10 && $debug);
    next if (substr($line,0,1) eq '@');
	last if (++$cnt > 10 && $debug);
	my @outvals = ();
	$outvals[$GFFLOC{source}] = $source;
	$outvals[$GFFLOC{type}] = $type;
	chomp($line);
#	print "line: '$line'\n";
	my @values = split /\t/, $line;
	my $flag = $values[1];
	$outvals[$GFFLOC{refmol}] = $values[2];
	$outvals[$GFFLOC{score}] = $values[4];
	$outvals[$GFFLOC{phase}] = '.';

    if ($onlymapped) {
        next if ($flag & 0x0004);
    }

	print "$values[0]\t" if ($debug);
	$outvals[$GFFLOC{start}] = $values[3];
	print "$outvals[$GFFLOC{start}]\t" if ($debug);

  if ($flag & 0x0010) {
#  if ($flag & RFLAGS->{REVERSED}) { # this doesn't work for some reason
    print "-\t" if ($debug);
    $outvals[$GFFLOC{strand}] = "-";
  } else {
    print "+\t" if ($debug);
    $outvals[$GFFLOC{strand}] = "+";
  }

	
	print "$values[5]\t" if ($debug);
	my $CIGAR = $values[5];

  if ($split) { # for discontinuous matches, explicitly split the match segments into HSP's
  
  } else {

    foreach my $pos (split/\D/, $CIGAR) {
      $values[3] += $pos;
    }
    $outvals[$GFFLOC{stop}] = $values[3] - 1;
    print "$outvals[$GFFLOC{stop}]\t" if ($debug);
    
#     if ($flag & 0x0010) {
#   #  if ($flag & RFLAGS->{REVERSED}) { # this doesn't work for some reason
#       print "-\t" if ($debug);
#       $outvals[$GFFLOC{strand}] = "-";
#     } else {
#       print "+\t" if ($debug);
#       $outvals[$GFFLOC{strand}] = "+";
#     }

  }
  
  #
  # for non-paired matches, we should have enough to print a GFF line
  #
  
  
	if ($flag & 0x0001) { # these are paired-end sequencing data
	
	  my $col9 = ucfirst($type) . " $values[0]";
	  
	  if ($flag & 0x0040) {
#	    $col9 .= "/a";
        if ($integers) {
            $col9 .= "/1";
        } else {
            $col9 .= "/a";
        }
	  } elsif ($flag & 0x0080) {
#	    $col9 .= "/b";
        if ($integers) {
            $col9 .= "/2";
        } else {
            $col9 .= "/b";
        }
	  }

#    if ($mapped_id_file) {
#        print IDFILE "$col9\n";
#    }    
	  
	  $outvals[$GFFLOC{group}] = $col9;
	  
	  if ($debug) {
          print "PAIRED $values[7]\t" if ($debug);
        
          if ($flag & 0x0020) {
        #  if ($flag & $RFLAGS->{M_REVERSED}) {
            print "-\t" if ($debug);
          } else {
            print "+\t" if ($debug);
          }
          
          print "$values[8]\t" if ($debug);
        
          print "[";
        
          print " read is mapped in a proper pair," if ($flag & 0x0002);
          print " query sequence is unmapped," if ($flag & 0x0004);
          print " mate is unmapped," if ($flag & 0x0008);
          print " read is first in a pair," if ($flag & 0x0040);
          print " read is second in a pair," if ($flag & 0x0080);
          print " alignment is not primary," if ($flag & 0x0100);
          print " read fails platform quality checks," if ($flag & 0x0200);
          print " read is PCR duplicate," if ($flag & 0x0400);
        
          print "]";
        }

	} else {
        my $col9 = ucfirst($type) . " $values[0]";
        if ($integers) {
            $col9 .= "/1";
        } else {
            $col9 .= "/a";
        }
        $outvals[$GFFLOC{group}] = $col9;

	}
	print "\n" if ($debug);
	
    if ($mapped_id_file) {
#        my $id = substr($outvals[$GFFLOC{group}],index($outvals[$GFFLOC{group}]," "));
        my $id = substr($outvals[$GFFLOC{group}],index($outvals[$GFFLOC{group}]," ") + 1);
#        print IDFILE $outvals[$GFFLOC{group}] . "\n";
        print IDFILE "$id\n";
    }
    gffline(@outvals);

}

close(IN);
close(OUT);
close(IDFILE) if ($mapped_id_file);

if ($sort_by_refmol) {
    #system("sort -k 1,4 $outfile > tempoutfile");
    print "sorting by refmol\n" if ($debug);
    my $tmp = File::Temp->new(
                                TEMPLATE => 'tempXXXXX',
                                DIR => '.',
                                SUFFIX => '.dat'
                            );
    my $tmpname = $tmp->filename();
    #system("sort -k 1,1 -k 4,4n $outfile > tempoutfile");# sort by col 1, then do numeric sort on col 4 (start coordinate)
    #system("sort --temporary-directory=$temp_directory -k 1,1 -k 4,4n $outfile > tempoutfile");# sort by col 1, then do numeric sort on col 4 (start coordinate)
    system("sort --temporary-directory=$temp_directory -k 1,1 -k 4,4n $outfile > $tmpname");# sort by col 1, then do numeric sort on col 4 (start coordinate)
    #system("cp tempoutfile $outfile");
    #unlink('tempoutfile');
    system("cp $tmpname $outfile");
    unlink($tmpname);
}

sub gffline {
  my @vals = @_;
#  print "\n" if ($debug);
    my $gffstring = join "\t", @vals;
    print $gffstring . ";\n" if ($debug);
    print OUT$gffstring . ";\n"; 
#	print join "\t", @vals, ";\n" if ($debug);
#    print OUT join "\t", @vals, ";\n";
#	print "\n";

}

sub _help {

print <<HELP;

            "infile=s"        =>  \$infile,
            "outfile=s"       =>  \$outfile,
            "source=s"        =>  \$source,
            "type=s"          =>  \$type,
            "debug"           =>  \$debug,
            "split"           =>  \$split,
            "integers"        =>  \$integers, # print /1 or /2 suffix to read ID instead of /a or /b
            "mapped"          =>  \$mapped_id_file,
            "onlymapped"      =>  only generate output for mapped reads
                                    some software outputs mapped & unmapped
                                    reads in the sam file (ie, bowtie)
            "sort"            =>  sort_by_refmol,# sort output file by refmol
                                    GFF files *must* be sorted like this if they will
                                    be used as input to QuantDisplay preload script
            "tempdir"         =>  temporary directory to use for sort command
            "help"            =>  help

HELP
exit();
}
