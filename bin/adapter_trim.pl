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
# copyright Scott Givan, University of Missouri, 2010
#
use strict;
use Getopt::Long;

my ($infile, $outfile, $overwrite, $adapterseq, $help, $verbose, $debug, $fastq,$idlist,$idlistfile,$printall,$notwoadapters,$id2adapters,$diversity);
my ($prefix,$internal,$suffix,$minlength,$lowercase);

my $options = GetOptions(
    "infile=s"      =>  \$infile,
    "outfile=s"     =>  \$outfile,
    "overwrite"     =>  \$overwrite,
    "prefix"        =>  \$prefix,# default
    "suffix"        =>  \$suffix,
    "minlength=i"   =>  \$minlength,
    "adapterseq=s"  =>  \$adapterseq,
    "lowercase"     =>  \$lowercase,
    "help"          =>  \$help,
    "verbose"       =>  \$verbose,
    "debug"         =>  \$debug,
    "fastq"         =>  \$fastq,
    "idlist"        =>  \$idlist,
    "idlistfile=s"  =>  \$idlistfile,
    "printall"      =>  \$printall,
    "notwoadapters" =>  \$notwoadapters,
    "id2adapters"   =>  \$id2adapters,
    "diversity=f"   =>  \$diversity,
);

if (!$options) {
  print "couldn't parse options\n";
  exit();
}

if ($help) {
print <<HELP;

    "infile=s"      =>  input file [default = infile]
    "outfile=s"     =>  output file [default = outfile]
    "overwrite"     =>  overwrite output file [default = no overwrite]
    "prefix"        =>  identify adapter seq at 5' end
    "suffix"        =>  identify adapter seq at 3' end
    "minlength"     =>  minimum acceptable length after trimming [default = 10]
    "adapterseq=s"   => the adapter sequence to identify and trim
    "lowercase"     =>  trim lowercase nucleotides & their quality values
    "help"          =>  print the help menu
    "verbose"       =>  verbose output to terminal
    "debug"         =>  debugging output to terminal
Options below only affect fastq input files:
    "fastq"         =>  input file is fastq [default = fastq]
    "idlist"        =>  generate a file called idlist that contains id's of sequences trimmed
    "idlistfile=s"  =>  input idlist file
    "printall"      =>  print all sequences, even if they don't contain adapter sequence
    "notwoadapters"  => skip sequences that contain > 1 adapter sequence, even if --printall
    "id2adapters"    => output sequences with >1 adapter to STDERR, requires --notwoadapters
    "diversity=f"    => diversity filter [default = 0 -- no diversity filter]

notes:

--diversity flag accepts a decimal value that specifies a nucleotide proportional threshold above which
a read is discarded. For example, if you specify 0.90, a read will be discarded if a any single nucleotide
composes more than 90% of the length of the read.
--notwoadapters filters out sequences with more than one adapter sequence. Typically, these reads are composed
almost exclusively of adapter sequence, so "two" usually means "at least two". You can see the sequences that
are filtered using --id2adapters, which directs those sequences to STDERR.

ie:  adapter_trim.pl --infile set1.fq --fastq --outfile trimmed1.fq --adapterseq AAGCAGTGGTATCAACGCAGAGTAC --notwoadapters --id2adapters >& 2adapters.txt

NOTE: Recent development has focused almost entirely on working with fastq files -- please do not use this script with fasta files.

HELP
exit;
}

$infile = 'infile' unless ($infile);
$outfile = 'outfile' unless ($outfile);
$idlistfile = 'idlist' unless ($idlistfile);
$diversity = 0.00 unless ($diversity);
#$adapterseq = '';
my $adapterseq_length = length($adapterseq);
$prefix = 1  unless ($suffix);
$minlength ||= 10;
$fastq ||= 1;

if (-e $outfile && !$overwrite) {
  print "$outfile already exists and you didn't specify to overwrite\n";
  exit();
}

#if (!$adapterseq) {
#  print "you must enter a adapter sequence using the --adapterseq argument\n";
#  exit();
#}

open(IN,$infile) or die "can't open '$infile': $!";
open(OUT,">$outfile") or die "can't open '$outfile': $!";

if ($idlist) {
  open(ID,">$idlistfile") or die "can't open 'idlist': $!";
}

my ($seqname, $seq, $parsed);

if (!$fastq) {
  while (<IN>) {
    print $_ if ($debug);
    
    if ($_ =~ /^>(.+)\n/) {
      $seqname = $1;
      next;
    }
    
    chomp($_);
    $seq = \$_;
  #  print "seq '$seqname' = '" . $$seq . "'\n" if ($debug);
    
    if ($adapterseq && $$seq =~ /^$adapterseq/) {
      print "parsed read:\n" if ($debug);
      $parsed = $$seq;
      $parsed =~ s/^$adapterseq//;
      print "\t>$seqname\n\t$parsed\n" if ($debug);
      print OUT ">$seqname\n$parsed\n";
    }
    
  }
} else { # this is a fastq file

  my ($initial,$gotseqname,$gotqual,$getqual,$getseq,$trim) = (0,0,0,0,0);
  my ($quality,$sequence,$qualname);
  #
  # loop through file, collect sequence ID, sequence and quality string
  # trim sequence if necessary
  # if sequence gets trimmed, trim same number of characters from quality string
  # print trimmed fastq sequence and quality string
  #
  
  my $loopcnt = 0;
  my $noread = 0;
  while (<IN>) {    
    if ($_ =~ /^\@(.+)\n/) {
      $seqname = $1;
       
    $sequence = <IN>;
    chomp($sequence);
    $qualname = <IN>;
    $qualname =~ s/[\+\n]//g;
    $quality = <IN>;
    chomp($quality);

    print "\nBEFORE:\nseqname:\t'$seqname'\nsequence:\t'$sequence'\nqualname:\t'$qualname'\nquality:\t'$quality'\n\n" if ($debug);
    #last if (++$loopcnt == 26 && $debug);

#        if ($debug) {
#            print "sequence: '$sequence'\n";
#            print "quality:  '$quality'\n";
#        }

    if ($diversity) {
        my $length = length($sequence);
        my $cA = $sequence =~ tr/A/A/;
        my $cT = $sequence =~ tr/T/T/;
        my $cG = $sequence =~ tr/G/G/;
        my $cC = $sequence =~ tr/C/C/;
        if ($cA/$length >= $diversity || $cT/$length >= $diversity || $cG/$length >= $diversity || $cC/$length >= $diversity) {
            if ($debug) {
                print "discarding\n$seqname\n$sequence\ndue to low diversity\n";
                print "\tA: $cA\n\tT: $cT\n\tG: $cG\n\tC: $cC\n";
            }
            next;
        }
    }# end of diversity section

    if ($lowercase) {
        if ($sequence =~ /([a-z]+)/) {
#            print "identified lowercase string '$1' in $seqname\n" if ($debug);
            $adapterseq = $1;
            $adapterseq_length = length($1);
        }
#        next unless ($adapterseq_length);
    }

      my $loc = index($sequence,$adapterseq,0);
    
    if ($debug) {
        print "adapter: '$adapterseq'\nadapterseq_length: '$adapterseq_length'\nloc = '$loc'\n";
    }

      #if ($adapterseq_length && $loc >= 0) {
      if ($adapterseq_length && (($prefix && $loc == 0) || ($suffix && $loc > 0))) {
      
        #if ($loc == 0 && $prefix) {
        if ($prefix) {
            print "prefix search\n" if ($debug);
          $sequence = substr($sequence,$adapterseq_length);
          $quality = substr($quality,$adapterseq_length);

          # notwoadapters only applies under prefix rules
          # suffix rules will trim everything downstream of first occurrence
          if ($notwoadapters) {
              #next if (!$printall && $sequence =~ /$adapterseq/);
              #next if (!$printall && (index($sequence,$adapterseq) >= 0));
              if (index($sequence,$adapterseq) >= 0) {
                  print STDERR "\@$seqname\n$sequence\n\+$seqname\n$quality\n" if ($id2adapters);
                  next;
              }
          }

          print "AFTER: seqname:\t'$seqname'\nsequence:\t'$sequence'\nqualname:\t'$qualname'\nquality:\t'$quality'\n\n" if ($debug);
          
          if (length($sequence) >= $minlength) {
              print OUT "\@$seqname\n$sequence\n\+$seqname\n$quality\n"; 
              print ID "$seqname\n" if ($idlist);
          }

        #} elsif ($loc > 0 && $suffix) {
        } elsif ($suffix) {

            print "suffix search\n" if ($debug);
          #print "BEFORE: seqname:\t'$seqname'\nsequence:\t'$sequence'\nqualname:\t'$qualname'\nquality:\t'$quality'\n\n" if ($debug);

          $sequence = substr($sequence,0,$loc);
          $quality = substr($quality,0,$loc);

          print "AFTER: seqname:\t'$seqname'\nsequence:\t'$sequence'\nqualname:\t'$qualname'\nquality:\t'$quality'\n\n" if ($debug);

          if (length($sequence) >= $minlength) {
              print OUT "\@$seqname\n$sequence\n\+$seqname\n$quality\n"; 
              print ID "$seqname\n" if ($idlist);
          }
        }
      
      } elsif ($printall) {

        if (length($sequence)) {
            print OUT "\@$seqname\n$sequence\n\+$seqname\n$quality\n"; 
        }
          
      } 
    } # end of if statement identifying adapter sequence in input sequence
  }
}


close(IN);
close(OUT);
close(ID) if ($idlist);
