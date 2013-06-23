#!/usr/bin/perl
use warnings;
use strict;
use Carp;
use Data::Dumper;
use Test::Simple qw(no_plan);
use lib '../lib';

use Diversity;

=head1 NAME

 <app name> -- one line description of of application's purpose

=head1 VERSION

This documentation refers to <app name> version 0.0.1

=head1 SYNOPSIS
=head1 REQUIRED ARGUMENTS
=head1 OPTIONS
=head1 DESCRIPTION
=head1 DIAGNOSTICS
=head1 CONFIGURATION AND ENVIRONMENT
=head1 DEPENDENCIES
=head1 INCOMPATIBILITIES
=head1 BUGS AND LIMITATIONS

There are no known bugs in this module. 
Please report problems to Quinn Weaver <quinn@fairpath.com>
Patches are welcome.

=head1 AUTHOR

Fernando J. Pineda <fernando.pineda@jhu.edu>

=head1 LICENSE AND COPYRIGHT

Copyright (c) 2012 Fernando J. Pineda (<fernando.pineda@jhu.edu>). All rights reserved.

This module is free software; you can redistribute it and/or
modify it under the same terms as Perl itself. See L<perlartistic>.
 
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

=cut

=begin
# Delete this block if you don't need command line options

# USAGE: demo.pl --length=10 --file='data.txt' --verbose
# MORE:  http://perldoc.perl.org/Getopt/Long.html

use Getopt::Long;

my $string  = "blah,blah,blah"; # default string
my $length  = 0; 				 # default integer
my $verbose = undef;			 # default flag

my $the_options = GetOptions (
	"length=i" => \$length,		# numeric
	"file=s"   => \$string,		# string
	"verbose"  => \$verbose		# flag
	);
=cut


my $div_obj = Diversity->new(INDELS => 1);
# $div_obj->initialize("fastatmp.far");
$div_obj->initialize("test_data_small.far");

$div_obj->brute_force_diversity(); 
# die("debug exit");


# $div_obj->initialize("test_data.far");
my $d  = $div_obj->diversity();
my $s  = $div_obj->sigma();

#my $v = $div_obj->variance();

print("diversity = $d +/- $s\n");

ok($d eq '0.0105620784550071',"expected diversity = 0.0105620784550071, got diversity=$d");

my $c = $div_obj->consensus();
warn("$c\n");

exit;
