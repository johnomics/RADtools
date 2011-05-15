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
# 02/11/10 Initial version
# 11/11/10 1.1.1 28 tests written to test input options
# 15/05/11 1.2.1 Replace -e option with -c

use strict;
use warnings;

use Test::More tests=>26;
use Test::Trap
  qw(trap $trap :flow :stderr(systemsafe) :stdout(systemsafe) :warn :exit);

# Load script
ok( require('RADtags'), 'loaded RADtags okay' ) or exit;

# Test basic options
my @r;

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
@ARGV = ( '-z' );
@r = trap { main() };
is( $trap->exit, 8, 'FASTQ specified with no site given' );

# Negative
@ARGV = ('-c-10');
@r = trap { main() };
like(
    $trap->die,
    qr/Please specify a positive value for cluster distance/,
    'Reject negative expected tags'
);

# Empty
@ARGV = ('-c');
@r = trap { main() };
is( $trap->exit, 3, 'Empty cluster distance value' );

# Test site option
# Empty
@ARGV = ( '-s' );
@r = trap { main() };
is( $trap->exit, 3, 'Empty restriction site value' );

# Non-ACGT
@ARGV = ( '-s1a!X' );
@r = trap { main() };
like(
    $trap->die,
    qr/should only contain AGCT/,
    'Reject restriction site containing non-ACGT characters'
);

chdir('data');

# Test directory option
@ARGV = ( '-dfakedir' );
@r = trap { main() };
like(
    $trap->die,
    qr/Can\'t open directory fakedir!/,
    'Fail on unknown directory'
);

@ARGV = ( '-dBeatles_tags', '-v' );
@r = trap { main() };
like( $trap->stdout, qr/Done/, 'Test successful run' );

# Longer than read length
@ARGV = (
    '-dBeatles_tags',
    '-sAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'
);
@r = trap { main() };
like(
    $trap->stdout,
    qr/is longer than read/,
    'Reject restriction site longer than read'
);

# Test file extension option

@ARGV = ( '-dBeatles_tags', '-v', '-ffakeext' );
@r = trap { main() };
like(
    $trap->die,
    qr/No read files found with extension 'fakeext'/,
    'Fail on unknown file extension'
);

chdir('Beatles_tags');
@ARGV = ( '-v' );
@r = trap { main() };
like( $trap->stdout, qr/Done/, 'Test current directory lookup' );

chdir('..');

# Test reads input with or without Illumina option
# If data is in reads format, this option makes no difference
# Assuming reads is only output by RADpools, and all RADpools data is Sanger

@ARGV = ( '-dBeatles_tags' );
my @output = (
"TGCAGGCACCAATGATGGATTTCGCTTGCATTACGTTCGTGGCAAA HHHHHHHHHHHHHHGHHHHHHHHHHHIHHHHHHHHIHHHIHAGGGG 3 1\n",
    "	ATCAGGTGTCCGATACCCATATCACAGGCTCTTACTAGCTTGGGGTCGGAT 3\n",
    "\n",
"TGCAGGGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTA HHHHHHHHHHHHHHHCEEE8DD\@DDGGGGGGGAFGGGGGE0<:\@:= 3 1\n",
    "	AACTCTAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGG 3\n",
    "\n"
);
match_file( \@ARGV, "Beatles_tags/John.tags", \@output,
    'Check Sanger reads output with no Illumina option' );

@ARGV = ( '-dBeatles_tags' );
@output = (
"TGCAGGGTTCACGGCAACTCCTTCGCTGAGCGTCTTACGAAGATTA hgggghgggdggggfggggggffgfhefeeggfegfgdgecd`a`a 3 1\n",
    "	ATTTTGGCGGAAGGTTGAAAGGTATTTTTTAACAAATTTTACGTTTACGAC 3\n",
    "\n",
"TGCAGGTGGCTGGATGAGGAAGGCCAATTGGGAATGTAATTACGGT fgffffggegggfggggggfgfffffeffedffff^bfffeb`]eB 3 1\n",
    "	AACAAAATAAAAATAGAGTAAGTACATTTTTTTGGTGATAACAGTCATGAT 3\n",
    "\n"
);
match_file( \@ARGV, "Beatles_tags/Paul.tags", \@output,
    'Check Illumina reads output with no Illumina option' );

@ARGV = ( '-dBeatles_tags', '-i' );
@output = (
"TGCAGGCACCAATGATGGATTTCGCTTGCATTACGTTCGTGGCAAA HHHHHHHHHHHHHHGHHHHHHHHHHHIHHHHHHHHIHHHIHAGGGG 3 1\n",
    "	ATCAGGTGTCCGATACCCATATCACAGGCTCTTACTAGCTTGGGGTCGGAT 3\n",
    "\n",
"TGCAGGGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTA HHHHHHHHHHHHHHHCEEE8DD\@DDGGGGGGGAFGGGGGE0<:\@:= 3 1\n",
    "	AACTCTAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGG 3\n",
    "\n"
);
match_file( \@ARGV, "Beatles_tags/John.tags", \@output,
    'Check Sanger reads output with Illumina option' );

@ARGV = ( '-dBeatles_tags', '-i' );
@output = (
"TGCAGGGTTCACGGCAACTCCTTCGCTGAGCGTCTTACGAAGATTA hgggghgggdggggfggggggffgfhefeeggfegfgdgecd`a`a 3 1\n",
    "	ATTTTGGCGGAAGGTTGAAAGGTATTTTTTAACAAATTTTACGTTTACGAC 3\n",
    "\n",
"TGCAGGTGGCTGGATGAGGAAGGCCAATTGGGAATGTAATTACGGT fgffffggegggfggggggfgfffffeffedffff^bfffeb`]eB 3 1\n",
    "	AACAAAATAAAAATAGAGTAAGTACATTTTTTTGGTGATAACAGTCATGAT 3\n",
    "\n"
);
match_file( \@ARGV, "Beatles_tags/Paul.tags", \@output,
    'Check Illumina reads output with Illumina option' );

# Test FASTQ->reads conversion
@ARGV = ( '-dBeatles_tags_fastq', '-z', '-f.fastq', '-sTGCAGG', );
@output = (
"TGCAGGCACCAATGATGGATTTCGCTTGCATTACGTTCGTGGCAAA ggggggggggggggfggggggggggghgggggggghggghg`ffff ATCAGGTGTCCGATACCCATATCACAGGCTCTTACTAGCTTGGGGTCGGAT gggdggggggggggggggggggggggggggggggggggggggggg_cL[`T\n",
"TGCAGGCCCATCAACGAAAATCGCATTGTTATTCGTTGAATGATAG Rc`aadggbgfa_f]V__[`TJTTQVXY`[`YYOccaS__BBBBBB CCTAACTTAAAATGTCTAATATTTCACATGACGGAATTAGGACATAAAAAT eee^eeff`fcffffeffffceeeefffffffeffffffffefffdfddff\n",
"TGCAGGCCCATCAACGAAAATCGCATTGTTATTCGTTGAATGATAG fgffggfggeggggggggggfgefffbfggggfffgefffcdYbfB AACAAAATAAAAATAGAGTAAGTACATTTTTTTGGTGATAACAGTCATGAT fgggggffgfgggggggggggfghgggggggggggcggfggfggggggggg\n",
"TGCAGGCCCATCAACGAAAATCGCATTGTTATTCGTTGAATGATAG ggggecccc^ccccbe`ce\\__^]Udcggedgeeedcgagagdgef CTCCTGCGGCTCTTTTGTGGGAAACATCGAAAAATCTTATTTGTGAAGATT fgegggfeegace`deZ^aegggeggdfgffdgaeg`dfgddgfgf]gdaf\n",
"TGCAGGGATCTGGAAATTCCTCAGGAGCCTTGGCGTGGGAAAACCT `e`cV^\\WYRagfccd]c[cZRVMS[KV__bfd`XXead^`BBBBB TGTGCACATTAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTC efbffgggggggfdgghggge`fffggggggfgggX]V][gegegeafffg\n",
"TGCAGGGATCTGGAAATTCCTCAGGAGCCTTGGCGTGGGAAAACCT fgffffggegggfggggggfgfffffeffedffff^bfffeb`]eB ATCGAAGTTAATAATTTTCCTGGCCTTTTTCAGGGCAAGAGATGCTCACGA eeeeeecffafggggghgggggggggggggggfggggeggdgggeeggggh\n",
"TGCAGGGATCTGGAAATTCCTCAGGAGCCTTGGCGTGGGAAAACCT hgggghgggdggggfggggggffgfhefeeggfegfgdgecd`a`a GGTTGCTCNCNTTCTGTGTCGCNTGTNNTCATACCTCGCGCATGCAGCACC BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB\n",
"TGCAGGGCCTCCTTCTGTTTAAGGATACAATATGTGATGGCAGCAT efffcgggedggefgegadge_dgcb\\ebaeb_e[egegegffa\\_ AAAGCGGAAGAGATCTCGACGACATCGTCGTAAGACATNNNNNNNNNNNNN ^[`__``W`^ZZUT]Yea\\c]W_V_ccYaBBBBBBBBBBBBBBBBBBBBBB\n",
"TGCAGGGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTA gggggggggggggggbdddWcc_ccfffffff`efffffdO[Y_Y\\ AACTCTAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGG fffffffffffffdfffffffffffdffdfffafffdfdffeefdfffffd\n",
"TGCAGGGTTCACGGCAACTCCTTCGCTGAGCGTCTTACGAAGATTA ggggggfecgXW`X_ffeffeb``efffdb_fYccYTZ_XOcXRVX ATTTTGGCGGAAGGTTGAAAGGTATTTTTTAACAAATTTTACGTTTACGAC TTTSOQTTSSL\\\\Q[NTTTSLTKTS[\\[QLW_VOXcW][c[Z`XcBBBBBB\n",
"TGCAGGTCTAGTTATTTCAGGGAATTCTCGAAGAAGTTCATGGAAT ggegggfggegggggggghggg_gggggggcbeeaggggggggggh AACTGACCACGCATCAAAAATGCAGTTCCTCGTGGACACTGGATCTGACTT ffghggggggggfggeeegfggggggggggggggggggggggafggg_ggd\n",
"TGCAGGTGATGGTGGAGCAAGCGTGGCAGCGGGGGTTCGATGTACA ggggggfgggggfggggggegggfghfhghggggfTceg^Y_BBBB CTTCCTGCTTTCTTTGTCTGGAAGGGATTATGTTTAGTAGTAAGGCTGCAA gggggggggggggggggggggggggggggggggggggggghggggggggeg\n",
);
match_file(
    \@ARGV,   "Beatles_tags_fastq/illumina.reads",
    \@output, 'Check FASTQ->reads conversion'
);

# Check processing Illumina/Sanger FASTQ data with/without Illumina option

@ARGV = ( '-dBeatles_tags_fastq', '-z', '-f.fastq', '-sTGCAGG' );
@output = (
"CCCATCAACGAAAATCGCATTGTTATTCGTTGAATGATAG fgcefccfbe`ce`__^]UdbfgedgeeedcfafadYbeB 3 3\n",
    "	AACAAAATAAAAATAGAGTAAGTACATTTTTTTGGTGATAACAGTCATGAT 1\n",
    "	CCTAACTTAAAATGTCTAATATTTCACATGACGGAATTAGGACATAAAAAT 1\n",
    "	CTCCTGCGGCTCTTTTGTGGGAAACATCGAAAAATCTTATTTGTGAAGATT 1\n",
    "\n",
"GATCTGGAAATTCCTCAGGAGCCTTGGCGTGGGAAAACCT ggedggfgfggggfgfffffefeedffef^edfecb`]`B 3 2\n",
    "	ATCGAAGTTAATAATTTTCCTGGCCTTTTTCAGGGCAAGAGATGCTCACGA 1\n",
    "	TGTGCACATTAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTC 1\n",
    "\n"
);
match_file( \@ARGV, "Beatles_tags_fastq/illumina.tags", \@output,
'Check Illumina FASTQ processing without Illumina option (outputs Illumina format)'
);

@ARGV = ( '-dBeatles_tags_fastq', '-z', '-f.fastq', '-sTGCAGG' );
@output = (
"CCCATCAACGAAAATCGCATTGTTATTCGTTGAATGATAG CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCBCCC 3 3\n",
    "	AACAAAATAAAAATAGAGTAAGTACATTTTTTTGGTGATAACAGTCATGAT 1\n",
    "	CCTAACTTAAAATGTCTAATATTTCACATGACGGAATTAGGACATAAAAAT 1\n",
    "	CTCCTGCGGCTCTTTTGTGGGAAACATCGAAAAATCTTATTTGTGAAGATT 1\n",
    "\n",
"GATCTGGAAATTCCTCAGGAGCCTTGGCGTGGGAAAACCT CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 3 2\n",
    "	ATCGAAGTTAATAATTTTCCTGGCCTTTTTCAGGGCAAGAGATGCTCACGA 1\n",
    "	TGTGCACATTAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTC 1\n",
    "\n",
);
match_file( \@ARGV, "Beatles_tags_fastq/sanger.tags", \@output,
'Check Sanger FASTQ processing without Illumina option (outputs Sanger format)'
);

# Check processing Illumina/Sanger FASTQ data with Illumina option
@ARGV =
  ( '-dBeatles_tags_fastq', '-z', '-f.fastq', '-sTGCAGG', '-i' );
@output = (
"CCCATCAACGAAAATCGCATTGTTATTCGTTGAATGATAG GHDFGDDGCFADFA@@?>6ECGHFEHFFFEDGBGBE:CF! 3 3\n",
    "	AACAAAATAAAAATAGAGTAAGTACATTTTTTTGGTGATAACAGTCATGAT 1\n",
    "	CCTAACTTAAAATGTCTAATATTTCACATGACGGAATTAGGACATAAAAAT 1\n",
    "	CTCCTGCGGCTCTTTTGTGGGAAACATCGAAAAATCTTATTTGTGAAGATT 1\n",
    "\n",
"GATCTGGAAATTCCTCAGGAGCCTTGGCGTGGGAAAACCT HHFEHHGHGHHHHGHGGGGGFGFFEGGFG?FEGFDCA>A! 3 2\n",
    "	ATCGAAGTTAATAATTTTCCTGGCCTTTTTCAGGGCAAGAGATGCTCACGA 1\n",
    "	TGTGCACATTAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTC 1\n",
    "\n"
);
match_file( \@ARGV, "Beatles_tags_fastq/illumina.tags", \@output,
'Check Illumina FASTQ processing with Illumina option (outputs Sanger format)'
);

@ARGV =
  ( '-c7', '-dBeatles_tags_fastq', '-z', '-f.fastq', '-sTGCAGG', '-i' );
@output = (
"CCCATCAACGAAAATCGCATTGTTATTCGTTGAATGATAG \$\$\#\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\!\$\$\$ 4 4\n",
    "	AACAAAATAAAAATAGAGTAAGTACATTTTTTTGGTGATAACAGTCATGAT 1\n",
    "	ATCAGGTGTCCGATACCCATATCACAGGCTCTTACTAGCTTGGGGTCGGAT 1\n",
    "	CCTAACTTAAAATGTCTAATATTTCACATGACGGAATTAGGACATAAAAAT 1\n",
    "	CTCCTGCGGCTCTTTTGTGGGAAACATCGAAAAATCTTATTTGTGAAGATT 1\n",
    "\n",
"GATCTGGAAATTCCTCAGGAGCCTTGGCGTGGGAAAACCT \$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$ 3 2\n",
    "	ATCGAAGTTAATAATTTTCCTGGCCTTTTTCAGGGCAAGAGATGCTCACGA 1\n",
    "	TGTGCACATTAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTC 1\n",
    "\n",
"TGATGGTGGAGCAAGCGTGGCAGCGGGGGTTCGATGTACA \$\$\$\$\$\"\"\$\"\$\$\$\$\$\$\$\$\"\"\$\$\$\"\"\$\$\$\$\$\"\$\$\$\"\$\$\$\"\$\" 2 2\n",
    "	AACTGACCACGCATCAAAAATGCAGTTCCTCGTGGACACTGGATCTGACTT 1\n",
    "	CTTCCTGCTTTCTTTGTCTGGAAGGGATTATGTTTAGTAGTAAGGCTGCAA 1\n",
    "\n"

);
match_file( \@ARGV, "Beatles_tags_fastq/sanger.tags", \@output,
'Check Sanger FASTQ processing with Illumina option (outputs incorrect qualities)'
);

sub match_file {
    my ( $ARGV_ref, $filename, $output_ref, $name ) = @_;
    @ARGV = @{$ARGV_ref};
    my @r = trap { main() };
    open my $fh, '<', $filename or die " Can't open $filename \n ";
    my @lines = <$fh>;
    close $fh;

    is_deeply( \@lines, $output_ref, $name );

}
