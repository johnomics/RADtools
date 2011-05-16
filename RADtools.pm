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

# RADtools.pm

# History:
# 11/11/10 1.1.1 Initial version
# 16/05/11 1.2.1 Minor fixes

package RADtools;

use strict;
use warnings;
use Exporter;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK);

use Cwd;
use File::Basename;
use File::Copy;

$VERSION   = 1.2.1;
@ISA       = qw(Exporter);
@EXPORT    = ();
@EXPORT_OK = qw(get_pools_filename sort_reads_file load_fastq_pair QUAL_OFFSET);

use constant QUAL_OFFSET => 33;

sub get_pools_filename {
    my ($directory) = @_;

    my $current_dir  = getcwd;
    my $base_dirname = basename($directory);

    my $pools_directory =
      $directory eq $current_dir ? "$current_dir/.." : $current_dir;

    my @pools_files;
    opendir( my $poolsdir, $pools_directory );
    @pools_files = grep { /^$base_dirname\.pools$/ } readdir $poolsdir;

    if ( @pools_files == 0 ) {
        die
"No pools file found! Please check current directory or specify -d option. Pools file should be in same directory as the directory specified by -d (which defaults to current directory)\n";
    }

    if ( @pools_files > 1 ) {
        die
"More than one pools file found! Please rename pools files or use a different -d option\n";
    }
    my $pools_filename = $pools_files[0];

    open my $pools_file, '<', "$pools_directory/$pools_filename"
      or die "Can't open pools file $pools_directory/$pools_filename\n";

    return $pools_file;
}

sub sort_reads_file {
    my ($arg_ref) = @_;

    my $directory     = $arg_ref->{directory};
    my $read_filename = $arg_ref->{read_filename};
    my $sort_start    = $arg_ref->{sort_start};
    my $read1_field   = $arg_ref->{read1_field};

    system(
"sort -t ' ' -k$read1_field\.$sort_start -T ./$directory --output=$read_filename\.sort $read_filename"
    );

    move( "$read_filename\.sort", $read_filename );

    return;
}

sub load_fastq_pair {
    my ($arg_ref) = @_;

    my $valid          = 1;
    my $read1_handle   = $arg_ref->{read1_handle};
    my $read2_handle   = $arg_ref->{read2_handle};
    my $qual_low_char  = $arg_ref->{qual_low_char};
    my $qual_high_char = $arg_ref->{qual_high_char};
    my $illumina       = $arg_ref->{illumina};

    my $r1_ref = load_fastq_record(
        {
            fastq_handle   => $read1_handle,
            qual_low_char  => $qual_low_char,
            qual_high_char => $qual_high_char,
            illumina       => $illumina
        }
    );
    return {
        r1name => '',
        r1seq  => '',
        r1qual => '',
        r2name => '',
        r2seq  => '',
        r2qual => '',
        valid  => -1
      }
      if ( $r1_ref->{valid} == -1 );

    # Default read 2 to empty, valid record, so single end data will load
    my $r2_ref = { name => '', seq => '', qual => '', valid => 1, };
    if ($read2_handle) {
        $r2_ref = load_fastq_record(
            {
                fastq_handle   => $read2_handle,
                qual_low_char  => $qual_low_char,
                qual_high_char => $qual_high_char,
                illumina       => $illumina
            }
        );
        return {
            r1name => '',
            r1seq  => '',
            r1qual => '',
            r2name => '',
            r2seq  => '',
            r2qual => '',
            valid  => -1
          }
          if ( $r2_ref->{valid} == -1 );
    }

    # Reject pair if either read is invalid
    $valid = 0 if ( ( $r1_ref->{valid} == 0 ) || ( $r2_ref->{valid} == 0 ) );

    # Check that r1 and r2 headers match; strip off final 1/2
    # If no read2 file, just make r2 stub eq to r1 stub so it always matches
    my $r1name_stub = substr $r1_ref->{name}, 0, -1;
    my $r2name_stub =
      $read2_handle
      ? substr $r2_ref->{name}, 0, -1
      : $r1name_stub;
    $valid = 0 if ( $r1name_stub ne $r2name_stub );

    return {
        r1name => $r1_ref->{name},
        r1seq  => $r1_ref->{seq},
        r1qual => $r1_ref->{qual},
        r2name => $r2_ref->{name},
        r2seq  => $r2_ref->{seq},
        r2qual => $r2_ref->{qual},
        valid  => $valid,
    };
}

sub load_fastq_record {
    my ($arg_ref) = @_;

    my $valid = 1;

    my $fastq_handle   = $arg_ref->{fastq_handle};
    my $qual_low_char  = $arg_ref->{qual_low_char};
    my $qual_high_char = $arg_ref->{qual_high_char};
    my $illumina       = $arg_ref->{illumina};

    my $seqhead  = <$fastq_handle>;
    my $seq      = <$fastq_handle>;
    my $qualhead = <$fastq_handle>;
    my $qual     = <$fastq_handle>;

    return { name => '', seq => '', qual => '', valid => -1, }
      if ( ( !$seqhead )
        || ( !$seq )
        || ( !$qualhead )
        || ( !$qual ) );

    chomp $seqhead;
    chomp $seq;
    chomp $qualhead;
    chomp $qual;

    my $name;
    if ( $seqhead =~ /^@(\S+)/xms ) {
        $name = $1;
    }
    else { $name = ''; $valid = 0; }

    $valid = 0 if ( $seq !~ /^[ACGTN]+$/xms );
    $valid = 0 if ( $qualhead !~ /^\+/xms );
    $valid = 0 if ( $qual !~ /^[$qual_low_char-$qual_high_char]+$/xms );
    if ($illumina) {
        $qual =~ tr/\x40-\xff\x00-\x3f/\x21-\xe0\x21/;
    }

    return {
        name  => $name,
        seq   => $seq,
        qual  => $qual,
        valid => $valid,
    };
}
