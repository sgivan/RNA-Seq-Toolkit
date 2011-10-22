#!/usr/bin/env perl
#
use strict;
use warnings;
#
# script to generate a GFF file of intergenic regions from a GFF file of genic regions.
#
my $debug = 1;

my ($cnt,$chrname,$chromcnt,@lastvals) = (0,'',0);
while (<>) {
    my $line1 = $_;
    my $line2 = <>;
    chomp($line1);
    chomp($line2);
    #print "\$line1 = '$line1'\n\$line2 = '$line2'\n\n" if ($debug);
    ++$cnt;
    #++$cnt;
    ++$chromcnt;
    #exit if ($debug && $cnt == 10);

    my @vals1 = split /\t/, $line1;
    my @vals2 = split /\t/, $line2;

    $vals1[8] = join("","id=intergenic.","XX",";name=intergenic.","XX",";");
    $vals2[8] = join("","id=intergenic.","XX",";name=intergenic.","XX",";");
    #$vals2[8] = "id=intergenic.$cnt;name=intergenic.$cnt;";

    if ($cnt == 1) {
    #if ($cnt <= 2) {
        $chrname = $vals1[0];
        $vals1[8] =~ s/XX/$cnt/g;
        print "$vals1[0]\t$vals1[1]\t$vals1[2]\t1\t$vals1[3]\t$vals1[5]\t$vals1[6]\t$vals1[7]\t$vals1[8]\n";
        $vals1[8] =~ s/$cnt/XX/g;
        ++$cnt;
    } else {
        if ($vals1[0] ne $chrname) { 
            print "using lastvals\n";
            $lastvals[8] =~ s/XX/$cnt/g;
            print "$lastvals[0]\t$lastvals[1]\t$lastvals[2]\t$lastvals[4]\t.\t$lastvals[5]\t$lastvals[6]\t$lastvals[7]\t$lastvals[8]\n";
            ++$cnt;
            $vals1[8] =~ s/XX/$cnt/g;
            print "$vals1[0]\t$vals1[1]\t$vals1[2]\t1\t$vals1[3]\t$vals1[5]\t$vals1[6]\t$vals1[7]\t$vals1[8]\n";
            $vals1[8] =~ s/$cnt/XX/g;
            ++$cnt;
            $chrname = $vals1[0];
        } elsif ($vals2[0] ne $chrname) {
            $vals1[8] =~ s/XX/$cnt/g;
            print "$vals1[0]\t$vals1[1]\t$vals1[2]\t$vals1[4]\t.\t$vals1[5]\t$vals1[6]\t$vals1[7]\t$vals1[8]\n";
            $vals1[8] =~ s/$cnt/XX/g;
            ++$cnt;
#            @vals2 = @vals1;
#            $line2 = <>;
#            chomp($line2);
#            @vals2 = split /\t/, $line2;
#            $vals2[8] = join("","id=intergenic.","XX",";name=intergenic.","XX",";");
#            $vals1[8] =~ s/XX/$cnt/g;
#            print "$vals1[0]\t$vals1[1]\t$vals1[2]\t$vals1[4]\t$vals2[3]\t$vals1[5]\t$vals1[6]\t$vals1[7]\t$vals1[8]\n";
#            $vals1[8] =~ s/$cnt/XX/g;
#            ++$cnt;
            $vals2[8] =~ s/XX/$cnt/g;
            print "$vals2[0]\t$vals2[1]\t$vals2[2]\t1\t$vals2[3]\t$vals2[5]\t$vals2[6]\t$vals2[7]\t$vals2[8]\n";
            $vals2[8] =~ s/$cnt/XX/g;
            ++$cnt;
            $chrname = $vals2[0];
            next;
        }

    }
    $vals1[8] =~ s/XX/$cnt/g;
    print "$vals1[0]\t$vals1[1]\t$vals1[2]\t$vals1[4]\t$vals2[3]\t$vals1[5]\t$vals1[6]\t$vals1[7]\t$vals1[8]\n";
    $vals1[8] =~ s/$cnt/XX/g;

    if (eof) {
        print "eof reached\n" if ($debug);
        ++$cnt;
        $vals1[8] =~ s/XX/$cnt/g;
        print "$vals1[0]\t$vals1[1]\t$vals1[2]\t$vals2[4]\t.\t$vals1[5]\t$vals1[6]\t$vals1[7]\t$vals1[8]\n";
        $vals1[8] =~ s/$cnt/XX/g;
    }
    @lastvals = @vals2;
}

