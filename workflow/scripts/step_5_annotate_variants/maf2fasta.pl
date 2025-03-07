#!/usr/bin/perl 
# from mugsy; multiple whole genome alignment tool
# from fork https://github.com/kloetzl/mugsy of 
# Angiuoli et al 2010; https://sourceforge.net/projects/mugsy/
# Convert MAF to FASTA 
# Optionally only convert blocks that contain label 
# ./maf2fasta.pl [label] < maf > fasta
# Modified by: Julia Höglund 2025-03-05

use strict; 
use warnings;  # Add warnings for better debugging

# Initialize variables
my $currscore; 
my $currlabel; 
my $currcoord; 
my $currorient; 
my $saveblock = 0; 
my @matches; 
my @blocks; 

# Read input line by line
while (my $line = <STDIN>) {
    if ($line =~ /a\s+score=([\d\.\-]+)/) {
        if ($saveblock > 0) {
            my @nmatches = @matches;
            push @blocks, [$currscore, $currlabel, $currorient, $currcoord, \@nmatches];
        }
        ($currscore) = ($line =~ /a\s+score=([\d\.\-]+)/);
        ($currlabel) = ($line =~ /label=(\d+)/);
        @matches = ();
    }
    elsif ($line =~ /s\s+(\S+)\s+(\d+)\s+(\d+)\s+([+-])\s+(\d+)\s+(\S+)/) {
        my $accession = $1;
        my $start = $2;
        my $len = $3;
        my $orientation = $4;
        my $seqlength = $5;
        my $seq = $6;
        if ($accession =~ /([^\.]+)\.(\S+)/) {
            # Add code to handle this case if needed
        }
        else {
            die "Invalid accession format: $accession" if ($accession =~ /\./);
        }
        push @matches, [$accession, $start, $len, $orientation, $seqlength, $seq];
        $saveblock = 1;
    }
}

# Save the last block if needed
if ($saveblock > 0) {
    my @nmatches = @matches;
    push @blocks, [$currscore, $currlabel, $currorient, $currcoord, \@nmatches];
}

# Sort and print blocks
foreach my $block (sort { $a->[3] <=> $b->[3] } @blocks) {
    if ($ARGV[0]) {
        if ($block->[1] && $block->[1] eq $ARGV[0]) {  # Ensure $block->[1] is defined
            printFASTA(@$block);
        }
    }
    else {
        printFASTA(@$block);
    }
}

# Subroutine to print FASTA format
sub printFASTA {
    my ($score, $label, $orient, $coord, $matches) = @_;
    foreach my $m (@$matches) {
        # print ">$m->[0].$label score=$score $m->[1] $m->[2] $m->[3] $m->[4]\n";
        print ">$m->[0] $m->[1] $m->[2] $m->[3] $m->[4]\n";
        for (my $i = 0; $i < length($m->[5]); $i += 60) {
            print substr($m->[5], $i, 60), "\n";
        }
    }
    print "=\n";
}
