#!/usr/bin/env perl
#
use strict;
use warnings;
#
# script to generate a GFF file of intergenic regions from a GFF file of genic regions.
#
my $debug = 0;

my ($cnt,$chrname,$prevline,@lastvals) = (0,'',0);
while (<>) {
    my $line1 = $_;
    chomp($line1);
    print "\n\$line1 = '$line1'\n" if ($debug);
    ++$cnt;

    my (@vals1,@vals2) = ();
    my ($coord1,$coord2);

    if ($cnt == 1) {
        @vals1 = split /\t/, $line1;
        $chrname = $vals1[0];
        $prevline = $line1;
        #$vals1[8] =~ s/XX/$cnt/g;
        print "$vals1[0]\t$vals1[1]\t$vals1[2]\t1\t$vals1[3]\t$vals1[5]\t$vals1[6]\t$vals1[7]\tid=intergenic.$cnt; name=intergenic.$cnt;\n";
        #$vals1[8] =~ s/$cnt/XX/g;
        #++$cnt;
        next;
    }
    
    @vals1 = split /\t/, $prevline;
    @vals2 = split /\t/, $line1;

    if ($vals1[0] ne $chrname) { 
        print "\$vals1[0] '$vals1[0]' ne '$chrname'\n" if ($debug);
        print "$vals1[0]\t$vals1[1]\t$vals1[2]\t1\t$vals1[4]\t$vals1[5]\t$vals1[6]\t$vals1[7]\tid=intergenic.$cnt; name=intergenic.$cnt;\n";
        ++$cnt;
        ($coord1,$coord2) = coordorder($vals1[4],$vals2[3]);
        #print "$vals1[0]\t$vals1[1]\t$vals1[2]\t$vals1[4]\t$vals2[3]\t$vals1[5]\t$vals1[6]\t$vals1[7]\ttid=intergenic.$cnt; name=intergenic.$cnt;\n";
        print "$vals1[0]\t$vals1[1]\t$vals1[2]\t$coord1\t$coord2\t$vals1[5]\t$vals1[6]\t$vals1[7]\ttid=intergenic.$cnt; name=intergenic.$cnt;\n";
        #++$cnt;
        $chrname = $vals1[0];
        $prevline = $line1;
        next;
    } elsif ($vals2[0] ne $chrname) {
        print "\$vals2[0] '$vals2[0]' ne '$chrname'\n" if ($debug);
        print "$vals1[0]\t$vals1[1]\t$vals1[2]\t$vals1[4]\t.\t$vals1[5]\t$vals1[6]\t$vals1[7]\tid=intergenic.$cnt; name=intergenic.$cnt;\n";
        ++$cnt;
        print "$vals2[0]\t$vals2[1]\t$vals2[2]\t1\t$vals2[3]\t$vals2[5]\t$vals2[6]\t$vals2[7]\tid=intergenic.$cnt; name=intergenic.$cnt;\n";
        #++$cnt;
        $chrname = $vals2[0];
        $prevline = $line1;
        next;
    }

    $vals1[8] = "id=intergenic.$cnt; name=intergenic.$cnt;";
    #print "$vals1[0]\t$vals1[1]\t$vals1[2]\t$vals1[4]\t$vals2[3]\t$vals1[5]\t$vals1[6]\t$vals1[7]\tid=intergenic.$cnt; name=intergenic.$cnt;\n";
    ($coord1,$coord2) = coordorder($vals1[4],$vals2[3]);
    print "$vals1[0]\t$vals1[1]\t$vals1[2]\t$coord1\t$coord2\t$vals1[5]\t$vals1[6]\t$vals1[7]\tid=intergenic.$cnt; name=intergenic.$cnt;\n";

    if (eof) {
        print "eof reached\n" if ($debug);
        ++$cnt;
        print "$vals1[0]\t$vals1[1]\t$vals1[2]\t$vals2[4]\t.\t$vals1[5]\t$vals1[6]\t$vals1[7]\tid=intergenic.$cnt; name=intergenic.$cnt;\n";
    }
    $prevline = $line1;
#    @lastvals = @vals2;
}

sub coordorder {
    my $coord1 = shift;
    my $coord2 = shift;

    if ($coord1 < $coord2) {
        return ($coord1,$coord2);
    } else {
        return ($coord2,$coord1);
    }
}

