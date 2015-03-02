#!/usr/bin/env perl

use 5.010;      # Require at least Perl version 5.8
use strict;
use autodie;
use Getopt::Long; # use GetOptions function to for CL args
use Webservice::InterMine;

# Set the output field separator as tab
$, = "\t";
# Print unicode to standard out
binmode(STDOUT, 'utf8');
# Silence warnings when printing null fields
no warnings ('uninitialized');

# The following import statement sets MouseMine as your default
# You must also supply your login details here to access this query
use Webservice::InterMine 0.9904 'http://www.mousemine.org/mousemine', 'M1R4W0S5Yfk6y7913eW1';

my $query = new_query(class => 'Gene');

# The view specifies the output columns
$query->add_view(qw/
    primaryIdentifier
    symbol
    name
    sequenceOntologyTerm.name
    chromosome.primaryIdentifier
    chromosomeLocation.id
    description
/);

$query->add_constraint(
    path  => 'Gene',
    op    => 'IN',
    value => '25genes',
    code  => 'A',
);

## Use an iterator to avoid having all rows in memory at once.
#my $it = $query->iterator();
#while (my $row = <$it>) {
#    print $row->{'primaryIdentifier'}, $row->{'symbol'}, $row->{'name'},
#        $row->{'sequenceOntologyTerm.name'}, $row->{'chromosome.primaryIdentifier'},
#        $row->{'chromosomeLocation.id'}, $row->{'description'}, "\n";
#}

my $results = $query->results(as => 'arrayrefs');

for my $row (@$results) {
    say $row->[0];
}

