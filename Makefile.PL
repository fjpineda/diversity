use inc::Module::Install;

name     'Diversity';
all_from 'lib/Diversity.pm';
author   q{Fernando J. Pineda <fernando.pineda@jhu.edu>};
license  'perl';

build_requires 'Test::More';

auto_install;

WriteAll;

