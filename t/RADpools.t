#!/usr/bin/env perl

# Copyright 2010 John Davey, University of Edinburgh john.davey@ed.ac.uk

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

#############################################################################

# RADpools test suite

# History:
# 05/08/10 Initial version
# 04/11/10 1.1.1 51 tests written across all functions
# 12/09/11 1.2.3 Added test for P2 MIDs
use strict;
use warnings;

use Test::More tests => 52;
use Test::Trap
  qw(trap $trap :flow :stderr(systemsafe) :stdout(systemsafe) :warn :exit);
use File::Path;

# Load script
ok( require('RADpools'), 'loaded RADpools OK' ) or exit;

# Test run with no parameters
my @r = trap { main() };
is( $trap->exit, 7, 'Print usage when run with no options' );

# Test usage message
@ARGV = ('-u');
@r = trap { main() };
is( $trap->exit, 4, 'Print usage when run with -u' );
@ARGV = ('--usage');
@r = trap { main() };
is( $trap->exit, 4, 'Print usage when run with --usage' );

# Test help message
@ARGV = ('-h');
@r = trap { main() };
is( $trap->exit, 5, 'Print help when run with -h' );
@ARGV = ('--help');
@r = trap { main() };
is( $trap->exit, 5, 'Print help when run with --help' );

# Test man page
@ARGV = ('-m');
@r = trap { main() };
is( $trap->exit, 3, 'Query -m option alone (should have max_processes value)' );

@ARGV = ('--man');
@r = trap { main() };
is( $trap->exit, 6, 'Print man page when run with --man' );

# Test required input options
@ARGV = ('-d dummy_name');
@r = trap { main() };
is( $trap->exit, 7, 'Input read file missing' );

@ARGV = ('-i dummy_name');
@r = trap { main() };
is( $trap->exit, 8, 'Directory missing' );

chdir 'data';

# TEST POOLS FILE
@ARGV = ( '-iBeatles.fastq', '-ddummyname' );
@r = trap { main() };
like(
    $trap->die,
    qr/Couldn't open 'dummyname.pools'/,
    'Non-existent pools file'
);

rmtree('Beatles');    # Remove any previous testing output

@ARGV = ( '-iBeatles.fastq', '-dBeatles', '-v' );
@r = trap { main() };
like(
    $trap->stdout,
    qr/Loading pools from Beatles.pools\.\.\./,
    'Recognise existing pools file'
);
like(
    $trap->stdout,
    qr/Created \.\/Beatles\.\.\./,
    'Create directory when none exists'
);
like( $trap->stdout, qr/4 records loaded/, 'Reads with true MIDs loaded' );
like(
    $trap->stdout,
    qr/12 records not matching any listed P1 MIDs/,
    'Reads with any mismatches in MIDs rejected when -f not specified'
);

@ARGV = ( '-iBeatles.fastq', '-dBeatles.pools', '-v', '-f' );
@r = trap { main() };
like(
    $trap->stdout,
    qr/Loading pools from Beatles.pools\.\.\./,
    'Strip .pools from directory argument'
);
unlike(
    $trap->stdout,
    qr/Created \.\/Beatles\.\.\./,
    'Don\'t recreate existing directory'
);
like( $trap->stdout, qr/14 records loaded/, 'Reads with fuzzy MIDs loaded' );
like(
    $trap->stdout,
    qr/2 records not matching any listed P1 MIDs/,
    'Reads with more than one mismatch or N rejected when -f specified'
);

@ARGV = ( '-iBeatles_fuzzyressite.fastq', '-dBeatles', '-v' );
@r = trap { main() };
like(
    $trap->stdout,
    qr/12 records loaded/,
    'Reads with fuzzy restriction sites loaded'
);
like(
    $trap->stdout,
    qr/2 records not matching restriction enzyme site/,
    'Reads with more than one mismatch in restriction site rejected'
);

@ARGV = ( '-iBeatles.fastq', '-pBeatles2.fastq', '-dBeatles_p2mids', '-v' );
@r = trap { main() };
like( $trap->stdout, qr/4 records loaded/, 'Reads with P2 MIDs loaded' );

# Test faulty pools file going in

my %mid_pools;
@r = trap { load_pools_and_mids( \%mid_pools, 'Beatles_faulty', 0, 0 ) };

# Only one field
like(
    $trap->stderr,
qr/Faulty pool: Ringo,TATAC - must have at least a name and one MID separated by whitespace/,
    'Reject pool lines with only one field'
);

# MIDs contain more than AGCT
like(
    $trap->stderr,
    qr/Faulty MID : Mal Evans - MIDs must contain only ACGT/,
    'Reject MIDs containing more than ACGT'
);

# Two pools have the same name
like(
    $trap->stderr,
    qr/Faulty pool: George - pool names must be unique/,
    'Reject duplicated pool names'
);

# Empty lines are ignored
is_deeply(
    \%mid_pools,
    {
        'CGATT' => { '' => ['Brian'] },
        'AGTCA' => { '' => ['Brian'] },
        'GACTA' => { '' => ['George'] },
        'CTGAC' => { '' => ['Mal'] },
        'CAACT' => { '' => ['Paul'] },
        'ATATC' => { '' => ['John'] }
    },
    'Faulty MIDs processed correctly'
);

# Test pools w only one or two mismatches between them
@ARGV = ( '-iBeatles_mismatch.fastq', '-dBeatles_mismatch', '-v' );
@r = trap { main() };
like(
    $trap->stdout,
    qr/4 records loaded/,
    'Fuzzy MIDs not loaded without -f option'
);

@ARGV = ( '-iBeatles_mismatch.fastq', '-dBeatles_mismatch', '-v', '-f' );
@r = trap { main() };
like( $trap->stdout, qr/7 records loaded/, 'Fuzzy MIDs loaded with -f option' );
like( $trap->stdout, qr/John.Paul/,        'MIDs 1bp apart merged' );
like( $trap->stdout, qr/George.Ringo/,     'MIDs 2bp apart merged' );

# FASTQ TESTS
@ARGV = ( '-ifaulty_RAD_reads.fastq', '-dBeatles', '-v' );
@r = trap { main() };
like(
    $trap->stdout,
    qr/4 records not in FASTQ format/,
    'Bad FASTQ format rejected'
);
like(
    $trap->stdout,
    qr/1 records not matching restriction enzyme site/,
    'Bad restriction site rejected'
);
like(
    $trap->stdout,
    qr/1 records not matching any listed P1 MIDs/,
    'Unknown MIDs rejected'
);

# QUALITY
@ARGV = ( '-iquality.fastq', '-dBeatles', '-v' );
@r = trap { main() };
like(
    $trap->stdout,
    qr/2 records loaded/,
    'By default, Illumina qualities accepted'
);
like(
    $trap->stdout,
    qr/1 records not in FASTQ format/,
    'By default, Sanger qualities rejected'
);
@ARGV = ( '-iquality.fastq', '-dBeatles', '-v', '-s' );
@r = trap { main() };
like(
    $trap->stdout,
    qr/1 records loaded/,
    'With -s, Sanger qualities accepted'
);
like(
    $trap->stdout,
    qr/2 records not in FASTQ format/,
    'With -s, Illumina qualities rejected'
);

# SINGLE / PAIRED END READ COMBINATIONS
@ARGV = ( '-iBeatles.fastq', '-pBeatles2.fastq', '-dBeatles', '-v' );
@r = trap { main() };
like(
    $trap->stdout,
    qr/4 records loaded/,
    'Paired end matches loaded correctly'
);
like(
    $trap->stdout,
    qr/1 records not in FASTQ format/,
    'Paired end reads with mismatched headers rejected'
);

# ENZYME SITE OPTION
@ARGV = ( '-iBeatles.fastq', '-dBeatles', '-v', '-eTGCA' );
@r = trap { main() };
like(
    $trap->stdout,
    qr/4 records loaded/,
    'Partial restriction enzyme still matches'
);
@ARGV = ( '-iBeatles.fastq', '-dBeatles', '-v', '-eTGCATT' );
@r = trap { main() };
like(
    $trap->stdout,
    qr/0 records loaded/,
    'Incorrect restriction enzyme given'
);
@ARGV = ( '-iBeatles.fastq', '-dBeatles', '-v', '-eTGCAGGGG' );
@r = trap { main() };
like( $trap->stdout, qr/0 records loaded/,
    'Overlong restriction enzyme given' );

# Quality threshold is between 0 and 40
@ARGV = ( '-iread1.fastq', '-dBeatles', '-v', '-q-1' );
@r = trap { main() };
like(
    $trap->die,
    qr/Quality threshold must be between 0 and 40/,
    'Reject quality threshold below 0'
);

@ARGV = ( '-iread1.fastq', '-dBeatles', '-v', '-q41' );
@r = trap { main() };
like(
    $trap->die,
    qr/Quality threshold must be between 0 and 40/,
    'Reject quality threshold over 40'
);

# Trimming is between 0 and length of read
@ARGV = ( '-iread1.fastq', '-dBeatles', '-v', '-t-1' );
@r = trap { main() };
like(
    $trap->stdout,
    qr/12 records loaded/,
    'Negative trim values set to read length'
);

@ARGV = ( '-iread1.fastq', '-dBeatles', '-v', '-t52' );
@r = trap { main() };
like(
    $trap->stdout,
    qr/12 records loaded/,
    'Trim values longer than read length set to read length'
);

# Spot check quality and trim combinations
@ARGV = ( '-iread1.fastq', '-dBeatles', '-v', '-q10', '-t45' );
@r = trap { main() };
like( $trap->stdout, qr/9 records loaded/, 'Spot check q=10, t=45' );
@ARGV = ( '-iread1.fastq', '-dBeatles', '-v', '-q10', '-t40' );
@r = trap { main() };
like( $trap->stdout, qr/12 records loaded/, 'Spot check q=10, t=40' );
@ARGV = ( '-iread1.fastq', '-dBeatles', '-v', '-q20', '-t40' );
@r = trap { main() };
like( $trap->stdout, qr/10 records loaded/, 'Spot check q=20, t=40' );

# Check reads output
@ARGV = ( '-iread1.fastq', '-pread2.fastq', '-dBeatles' );
@r = trap { main() };
open my $john_fh, '<', "Beatles/John.reads"
  or die "Can't open Beatles/John.reads\n";
my @john = <$john_fh>;
close $john_fh;
is_deeply(
    \@john,
    [
"TGCAGGCACCAATGATGGATTTCGCTTGCATTACGTTCGTGGCAAA HHHHHHHHHHHHHHGHHHHHHHHHHHIHHHHHHHHIHHHIHAGGGG ATCAGGTGTCCGATACCCATATCACAGGCTCTTACTAGCTTGGGGTCGGAT HHHEHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH\@D-<A5\n",
"TGCAGGGATCTGGAAATTCCTCAGGAGCCTTGGCGTGGGAAAACCT IHHHHIHHHEHHHHGHHHHHHGGHGIFGFFHHGFHGHEHFDEABAB GGTTGCTCNCNTTCTGTGTCGCNTGTNNTCATACCTCGCGCATGCAGCACC ###################################################\n",
"TGCAGGGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTA HHHHHHHHHHHHHHHCEEE8DD\@DDGGGGGGGAFGGGGGE0<:\@:= AACTCTAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGG GGGGGGGGGGGGGEGGGGGGGGGGGEGGEGGGBGGGEGEGGFFGEGGGGGE\n"
    ],
    'Check sorted reads output'
);

# Output FASTQ
@ARGV = ( '-iread1.fastq', '-pread2.fastq', '-dBeatles', '-o' );
@r = trap { main() };
open my $john1_fh, '<', "Beatles/John_1.fastq"
  or die "Can't open Beatles/John_1.fastq\n";
my @john1 = <$john1_fh>;
close $john1_fh;
is_deeply(
    \@john1,
    [
        "\@HWUSI-EAS721_0001:2:4:978:315#0/1\n",
        "TGCAGGCACCAATGATGGATTTCGCTTGCATTACGTTCGTGGCAAA\n",
        "+HWUSI-EAS721_0001:2:4:978:315#0/1\n",
        "HHHHHHHHHHHHHHGHHHHHHHHHHHIHHHHHHHHIHHHIHAGGGG\n",
        "\@HWUSI-EAS721_0001:2:4:978:1113#0/1\n",
        "TGCAGGGATCTGGAAATTCCTCAGGAGCCTTGGCGTGGGAAAACCT\n",
        "+HWUSI-EAS721_0001:2:4:978:1113#0/1\n",
        "IHHHHIHHHEHHHHGHHHHHHGGHGIFGFFHHGFHGHEHFDEABAB\n",
        "\@HWUSI-EAS721_0001:2:4:978:105#0/1\n",
        "TGCAGGGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTA\n",
        "+HWUSI-EAS721_0001:2:4:978:105#0/1\n",
        "HHHHHHHHHHHHHHHCEEE8DD\@DDGGGGGGGAFGGGGGE0<:\@:=\n"
    ],
    'Check Read 1 sorted FASTQ output'
);

open my $john2_fh, '<', "Beatles/John_2.fastq"
  or die "Can't open Beatles/John_2.fastq\n";
my @john2 = <$john2_fh>;
close $john2_fh;
is_deeply(
    \@john2,
    [
        "\@HWUSI-EAS721_0001:2:4:978:315#0/2\n",
        "ATCAGGTGTCCGATACCCATATCACAGGCTCTTACTAGCTTGGGGTCGGAT\n",
        "+HWUSI-EAS721_0001:2:4:978:315#0/2\n",
        "HHHEHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH\@D-<A5\n",
        "\@HWUSI-EAS721_0001:2:4:978:1113#0/2\n",
        "GGTTGCTCNCNTTCTGTGTCGCNTGTNNTCATACCTCGCGCATGCAGCACC\n",
        "+HWUSI-EAS721_0001:2:4:978:1113#0/2\n",
        "###################################################\n",
        "\@HWUSI-EAS721_0001:2:4:978:105#0/2\n",
        "AACTCTAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGG\n",
        "+HWUSI-EAS721_0001:2:4:978:105#0/2\n",
        "GGGGGGGGGGGGGEGGGGGGGGGGGEGGEGGGBGGGEGEGGFFGEGGGGGE\n"
    ],
    'Check Read 2 sorted FASTQ output'
);

