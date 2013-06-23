#!perl -T

use Test::More tests => 1;

BEGIN {
    use_ok( 'Diversity' ) || print "Bail out!\n";
}

diag( "Testing Diversity $Diversity::VERSION, Perl $], $^X" );
