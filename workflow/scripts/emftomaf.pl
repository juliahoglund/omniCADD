#!/usr/bin/env perl
# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


#####################################################################
##
## PROGRAM emf2maf.pl
##
## DESCRIPTION
##   This parser converts an EMF (Ensembl Multi Format) into an
##   MAF file.
##
#####################################################################


use strict;
use warnings;
use autodie;  # Automatically handle errors for file operations

my $VERSION = "1.0.1";

if (!@ARGV) {
  print STDERR qq"
  emf2maf.pl v$VERSION - EMF to MAF file converter.
  Use: emf2maf.pl file1.emf [file2.emf ...]
  MAF files will be named file1.maf, file2.maf...
";
  exit 1;
}

foreach my $emf_file (@ARGV) {
  print STDERR "Parsing file $emf_file...\n";
  my $maf_file = $emf_file;
  $maf_file =~ s/(\.emf)?$//;
  $maf_file .= ".maf";
  open(my $emf_fh, '<', $emf_file);
  open(my $maf_fh, '>', $maf_file);
  print $maf_fh "##maf version=1\n";
  print $maf_fh "# emf2maf.pl v$VERSION from file $emf_file\n";
  print $maf_fh "# Here is the header from the original file:\n";
  my $data = [];
  my $pattern = "";
  my $mode = "header";
  while ($_ = <$emf_fh>) {
    if ($_ =~ /^##/) {
     if ($_ =~ /^## ?DATE (.+)/) {
       print $maf_fh "# original dump date: $1\n";
     } elsif ($_ =~ /^## ?RELEASE (.+)/) {
       print $maf_fh "# ensembl release: $1\n";
     } elsif ($_ =~ /^### (.+)/) {
       print $maf_fh "# epo2x composite sequence: $1\n";
     }
    } elsif ($_ =~ /^# *(.+)/) {
      print $maf_fh "# emf comment: $1\n";
    } elsif ($_ =~ /^SEQ (.+)/) {
      my $info = $1;
      my ($species, $chromosome, $start, $end, $strand) =  $info =~ /(\S+)\s(\S+)\s(\S+)\s(\S+)\s(\S+)/;
      my ($chr_length) = $info =~ /chr_length=(\d+)/;
      $pattern .= " ?(\\S)";
      push(@$data, {type => "SEQ", species => $species, seq_region => $chromosome,
          start => $start, end => $end, strand => $strand, seq => "",
	  chr_length => $chr_length});
    } elsif ($_ =~ /^SCORE/) {
      $pattern .= " (\-?[\\d\.]+)";
      push(@$data, {type => "SCORE", values => []});
    } elsif ($_ =~ /^TREE (.+)/) {
      print $maf_fh "# tree: $1\n"
    } elsif ($_ =~ /^ID (.+)/) {
      print $maf_fh "# id: $1\n"
    } elsif ($_ =~ /^DATA/) {
      if ($mode eq "header") {
        $mode = "data";
      } else {
        die "Error while parsing line $.: Unexpected DATA line\n";
      }
    } elsif ($_ =~ /^\/\//) {
      if ($mode eq "data") {
        write_maf($maf_fh, $data);
	$data = [];
	$pattern = "";
        $mode = "header";
      } else {
        die "Error while parsing line $.: Unexpected end of block\n";
      }
    } elsif ($_ !~ /^\s*$/) {
      if ($mode eq "data") {
        my @this_line = $_ =~ /$pattern/;
	for (my $i=0; $i<@this_line; $i++) {
	  my $this_data = $data->[$i];
	  if ($this_data->{type} eq "SEQ") {
	    $this_data->{seq} .= $this_line[$i];
	  } elsif ($this_data->{type} eq "SCORE") {
	    push(@{$this_data->{values}}, $this_line[$i]);
	  }
	}
      }
    }
  }
  close($emf_fh);
  close($maf_fh);
}

sub write_maf {
  my ($maf_fh, $data) = @_;
  print $maf_fh "a\n";
  foreach my $this_data (@$data) {
    if ($this_data->{type} eq "SEQ") {
      if (!defined($this_data->{chr_length})) {
        die "Cannot write maf in reverse strand because there is no length info on EMF original file\n";
      }
      my ($maf_start, $maf_length, $maf_strand);
      if ($this_data->{strand}==1 or $this_data->{strand} eq "+") {
        $maf_start = $this_data->{start} - 1;
        $maf_strand = "+";
      } elsif ($this_data->{strand}==-1 or $this_data->{strand} eq "-") {
        $maf_start = $this_data->{chr_length} - $this_data->{end};
        $maf_strand = "-";
      } else {
        die "Cannot understand strand <".$this_data->{strand}.">";
      }
      $maf_length = $this_data->{end} - $this_data->{start} + 1;
      printf $maf_fh "s %-30s %10d %7d %s %10d ",
          ($this_data->{species}.".".$this_data->{seq_region}),
	  $maf_start, $maf_length, $maf_strand,
	  ($this_data->{chr_length} or 0);
      print $maf_fh $this_data->{seq}, "\n"
    } elsif ($this_data->{type} eq "SCORE") {
      print $maf_fh "# gerp scores: ", join(" ", @{$this_data->{values}}), "\n";
    }
  }
  print $maf_fh "\n";
}
