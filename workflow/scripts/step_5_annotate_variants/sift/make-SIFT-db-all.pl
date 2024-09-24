#!/usr/bin/perl -w

# This script's original distribution:
# https://github.com/pauline-ng/SIFT4G_Create_Genomic_DB/blob/master/make-SIFT-db-all.pl
# Modifications: Julia Höglund
# Date: 2024-09-23
# Changes:
# - add script directory to getoptions
# - add redirection to scripts when lookin for metadocs

# create the SIFT databases, must be run in scripts directory because it's lookin for metadocs


use strict;
#require 'common-utils.pl'; 
#use Getopt::Std;
use File::Basename;
use Cwd qw(abs_path);
my $directory_of_script = dirname(abs_path(__FILE__));
require $directory_of_script . '/common-utils.pl';

use Getopt::Long qw (GetOptions);

my %options  = ();
#getopts ("mhf:", \%options);

my $help = "";
my $ensembl_download = "";
my $metafile = "";
my $script_dir = "";
GetOptions ('help|h' => \$help, 
	    'ensembl_download' => \$ensembl_download,
	    'config=s' => \$metafile,
	    'script_directory=s' => \$script_dir);

if ($help || !$metafile || !$metafile || !$script_dir) {
	print "make-SIFT-db-all.pl\n";
	print " -h  help\n";
	print " -config  [required config file]\n";
	print " -ensembl_download \n";
	print " -script_directory  [required path to sift scripts]\n"
}

my $meta_href = readMeta($metafile);
my %meta_hash = %{$meta_href};

make_dir ($meta_hash{"PARENT_DIR"});
make_dir ($meta_hash{"PARENT_DIR"} . "/" . $meta_hash{"ORG_VERSION"});

for my $key (keys %meta_hash) {
	my $value = $meta_hash{$key};
#	print $value . "\n";
	if ($key =~ /_DIR$/ && $key ne "PARENT_DIR") {
#		print "making $value\n";
		my $dir = $meta_hash{"PARENT_DIR"} . "/" . $value;
		if (! -d $dir) { mkdir($dir, 0775) or die $dir . $!; }
	}
}
my $siftscore_dir = $meta_hash{"PARENT_DIR"} . "/" . $meta_hash{"SINGLE_REC_WITH_SIFTSCORE_DIR"};
if (! -d $siftscore_dir) { mkdir ($siftscore_dir, 0775); } 		


# changed from orinal; 
# store paths so this script does not need to be run in the script directory
my $ensembl_annotation = $script_dir . '/download-annotation-srcs.pl';
my $ensembl_fasta = $script_dir . '/download-fasta.pl';
my $ensembl_dbSNP = $script_dir . '/download-dbSNP-files.pl';
my $split_dbSNP = $script_dir . '/split-dbSNP-by-chr.pl';

if ($ensembl_download) {
	# not manual, this is from Ensembl and needs to download
	print "downloading gene annotation\n";
	`perl $ensembl_annotation $metafile`; 
	print "done downloading gene annotation\n"; 

	print "downloading fasta files\n"; 
	`perl $ensembl_fasta $metafile`; 
	print "done downloading DNA fasta sequences"; 

	print "download dbSNP files\n";
	`perl $ensembl_dbSNP $metafile`; 
	if (exists ($meta_hash{"DBSNP_VCF_FILE"}) && $meta_hash{"DBSNP_VCF_FILE"} ne "") 
	{
		print "splitting dbSNP file\n";
		`perl $split_dbSNP $metafile`;
	}
} # end if downloading files from Ensembl

my $convert_ucsc = $script_dir . '/gff_gene_format_to_ucsc.pl';

print "converting gene format to use-able input\n";
`perl $convert_ucsc $metafile`; 
#`perl ensembl_gene_format_to_ucsc.pl $metafile`;
print "done converting gene format\n";

# decompress chromosome files so that Bio::DB can process them
my $fasta_dir =  $meta_hash{"PARENT_DIR"} . "/" . $meta_hash{"CHR_DOWNLOAD_DEST"} ;
if (glob ("$fasta_dir/*.gz")) {
        $fasta_dir .= "/*.gz";
        `gunzip $fasta_dir`;
        # check that all DNA files finished downloading, otherwise do it again.
        # because all future programs won't work properly otherwise.
        if ($?) {
                die "DNA files do not exist or did not unzip properly\n";
        }
} elsif (glob ("$fasta_dir/*.fa") || glob ("$fasta_dir/*.fasta") ||  glob ("$fasta_dir/*.fna")) {
        # *.fa or *.fasta files already exist
} else {
        die "no DNA fasta files in $fasta_dir\n";
}

my $records = $script_dir . '/make-single-records-BIOPERL.pl';
my $noncoding = $script_dir . '/make-single-records-noncoding.pl';
my $fasta_seqs = $script_dir . '/generate-fasta-subst-files-BIOPERL.pl';

print "making single records file\n"; 
`perl $records $metafile`; 
print "done making single records template\n"; 

print "making noncoding records file\n";
`perl $noncoding $metafile`;
print "done making noncoding records\n";

print "make the fasta sequences\n"; 
`perl $fasta_seqs $metafile`; 
print "done making the fasta sequences\n";

print "start siftsharp, getting the alignments\n"; 

# remove all_prot.fasta if it already exists
my $all_prot_fasta  =  $meta_hash{"PARENT_DIR"} . "/all_prot.fasta";
if ( -e  $all_prot_fasta ) {
	unlink $all_prot_fasta;
}

my $combine_prot_fasta_command = "for file in " .  $meta_hash{"PARENT_DIR"} . "/" . $meta_hash{"FASTA_DIR"} . "/*.fasta ; do cat \"\$file\" >> " . $all_prot_fasta . "; done";

`$combine_prot_fasta_command`;


my $sift4g_command = $meta_hash{"SIFT4G_PATH"} .  " -d " . $meta_hash{"PROTEIN_DB"} . " -q " . $all_prot_fasta . " --subst " .  $meta_hash{"PARENT_DIR"} . "/" . $meta_hash{"SUBST_DIR"} . " --out " .  $meta_hash{"PARENT_DIR"} . "/" . $meta_hash{"SIFT_SCORE_DIR"} . " --sub-results " ;

print $sift4g_command . "\n";

#`$siftsharp_command`;
`$sift4g_command`;

if ($? != 0) {
        exit (-1);
}

print "done getting all the scores\n";

print "populating databases\n";
# check quality
# Pauline to do , by chromosame in script
my $populate_db = $script_dir . '/make-sift-scores-db-batch.pl';
my $check_genes = $script_dir . '/check_genes.pl'; 

`perl $populate_db $metafile`;

print "checking the databases\n";
system ("perl $check_genes $metafile");

my @chromosomes =  getChr ($meta_hash{"PARENT_DIR"} . "/" . $meta_hash{"GENE_DOWNLOAD_DEST"});
my $check_siftdb = $script_dir . '/check_SIFTDB.py';
my $check_pep = $script_dir . '/checkENSPep.pl';


foreach my $chr (@chromosomes) {

        my $db_file = $meta_hash{"PARENT_DIR"} . "/" . $meta_hash{"SINGLE_REC_WITH_SIFTSCORE_DIR"} . "/" . $chr .  "_scores.Srecords.with_dbSNPid.sorted";
        my $check_outfile = $meta_hash{"PARENT_DIR"}  .  "/" . $meta_hash{"ORG_VERSION"} . "/" . $chr .  "_SIFTDB_stats.txt";
        `python $check_siftdb $db_file $check_outfile`;
}

my $pep_check_file = $meta_hash{"PARENT_DIR"} . "/ENS_peptide_check.txt";
`perl $check_pep $metafile ENS $pep_check_file`; 

my $final_outfolder = $meta_hash{"PARENT_DIR"} ."/". $meta_hash{"ORG_VERSION"};
system ("cp  $metafile $final_outfolder"); 

# make some space
my $zip_dir = $meta_hash{"PARENT_DIR"} . "/" . $meta_hash{"CHR_DOWNLOAD_DEST"} . "/*" ;
print "zipping up $zip_dir\n";
system ("gzip $zip_dir");
my $rm_dir = $meta_hash{"PARENT_DIR"} . "/" . $meta_hash{"SINGLE_REC_WITH_SIFTSCORE_DIR"} . "/*";
system ("rm $rm_dir");

print "All done!\n";
exit (0);

