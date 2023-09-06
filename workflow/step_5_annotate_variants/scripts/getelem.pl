#!/usr/bin/env perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
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

# :Edited: 18-11-2022 and 24-08-2023

# note!
# not used in current version of pipeline

use strict;
use warnings;

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Utils::Exception qw(throw);

# Imports Gepopt package which is doing the command line handling
use Getopt::Long;

# Save arguments in the variable names utilized by this script
GetOptions( "chr=s" => \my $seq_region
          , "start=i" => \my $seq_region_start
          , "end=i" => \my $seq_region_end
          , "sp=s" => \my $species
          );

# @ARGV: [ qw<arg1 arg2 arg3> ]


#
# Simple example to show how to get conservation scores for a slice.
#

my $reg = "Bio::EnsEMBL::Registry";

$reg->load_registry_from_db(
      -host => "ensembldb.ensembl.org",
      -user => "anonymous",
);

# Get method_link_species_set adaptor
my $mlss_adaptor = $reg->get_adaptor("Multi", "compara", "MethodLinkSpeciesSet");

# Get method_link_species_set object for gerp conservation scores for mammals
my $mlss = $mlss_adaptor->fetch_by_method_link_type_species_set_name("GERP_CONSTRAINED_ELEMENT", "mammals");
throw("Unable to find method_link_species_set") if (!defined($mlss));

# Get slice adaptor for $species
my $slice_adaptor = $reg->get_adaptor($species, 'core', 'Slice');
throw("Registry configuration file has no data for connecting to <$species>") if (!$slice_adaptor);

# Create slice 
my $slice = $slice_adaptor->fetch_by_region('toplevel', $seq_region, $seq_region_start, $seq_region_end);
throw("No Slice can be created with coordinates $seq_region:$seq_region_start-$seq_region_end") if (!$slice);

# Get conservation score adaptor
my $cs_adaptor = $reg->get_adaptor("Multi", 'compara', 'ConstrainedElement');   
                 
# To get one score per base in the slice, must set display_size to the size of the slice.
my $scores = $cs_adaptor->fetch_all_by_MethodLinkSpeciesSet_Slice($mlss, $slice);

print "number of scores " . @$scores . "\n";

# Print out the position, observed, expected and difference scores.
foreach my $score (@$scores) {
        printf("%d\t%d\t%d\t%.4f\t%.4e\n",  $seq_region, $score->start, $score->end, $score->score, $score->p_value);
}




