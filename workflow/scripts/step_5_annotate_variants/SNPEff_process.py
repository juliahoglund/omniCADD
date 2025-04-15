#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Updated script to process VCF input format with `ANN` annotations while retaining all functionality, including GC% and CpG% calculations.
"""

import csv
from argparse import ArgumentParser
import pysam
import sys

# OptionParser for input files.
parser = ArgumentParser(description=__doc__)

parser.add_argument("-i", "--input",
    help="Path to the VCF input file",
    type=str,
    required=True)
parser.add_argument("-r", "--reference",
    help="Path to the reference genome (FASTA file)",
    type=str,
    required=True)
parser.add_argument("-o", "--output",
    help="Output file (default: vep_annotated.tsv)",
    type=str,
    default="vep_annotated.tsv")

args = parser.parse_args()

##########################################
############### VARIABLES ################
##########################################

# Define headers for output file.
ELIST = ['#Chrom', 'Pos', 'Ref', 'Alt', 'isTv', 'Consequence', 'GC', 'CpG',
         'motifECount', 'motifEHIPos', 'motifEScoreChng', 'Domain', 'oAA',
         'nAA', 'Grantham', 'SIFTcat', 'SIFTval', 'cDNApos', 'relcDNApos',
         'CDSpos', 'relCDSpos', 'protPos', 'relprotPos']

# List for transversions and transitions.
TRANSVERSIONS = {('A', 'C'), ('C', 'A'), ('T', 'A'), ('A', 'T'), ('C', 'G'),
                 ('G', 'C'), ('G', 'T'), ('T', 'G')}

TRANSITIONS = {('C', 'T'), ('T', 'C'), ('G', 'A'), ('A', 'G')}

# List of hierarchy of the consequences.
HIERARCHY = [
    "splice_acceptor_variant", "splice_donor_variant", "stop_gained",
    "frameshift_variant", "stop_lost", "start_lost", "missense_variant",
    "inframe_insertion", "inframe_deletion", "synonymous_variant",
    "stop_retained_variant", "coding_sequence_variant", "mature_miRNA_variant",
    "5_prime_UTR_variant", "3_prime_UTR_variant", "intron_variant",
    "splice_region_variant", "downstream_gene_variant", "upstream_gene_variant",
    "intergenic_variant"
]

##########################################
############### FUNCTIONS ################
##########################################

def parse_annotation(info_field):
    """
    Parse the `ANN` field in the INFO column to extract the most severe consequence.
    """
    annotations = []
    for entry in info_field.split(";"):
        if entry.startswith("ANN="):
            annotations = entry[4:].split(",")  # Extract annotations after "ANN="
            break

    # Parse each annotation and find the most severe consequence.
    for annotation in annotations:
        fields = annotation.split("|")
        if len(fields) > 1:
            consequence = fields[1]  # The second field is the consequence.
            # Match the consequence against the HIERARCHY
            for severity in HIERARCHY:
                if severity in consequence:
                    return consequence  # Return the full consequence description
    return "Unknown"  # Default to "Unknown" if no known consequence is found

def is_transversion(ref, alt):
    """
    Check if the variant is a transversion.
    """
    return (ref, alt) in TRANSVERSIONS

def count_GC_CpG(chrom, start, end, window, seq_tabix):
    """
    Count GC and CpG sites in a window of 75 bases (75 bases before and after the variant).
    Returns percentage GC and CpG counts inside this window.
    """
    # Normalize chromosome name to match reference genome
    if chrom not in seq_tabix.references:
        if f"chr{chrom}" in seq_tabix.references:
            chrom = f"chr{chrom}"  # Prepend 'chr'
        elif f"chr_{chrom}" in seq_tabix.references:
            chrom = f"chr_{chrom}"  # Prepend 'chr_'
        elif chrom.startswith("chr") and chrom[3:] in seq_tabix.references:
            chrom = chrom[3:]  # Remove 'chr' prefix
        else:
            print(f"Chromosome {chrom} not found in reference genome.")
            return '-', '-'

    try:
        # Fetch sequence from the reference genome
        sequence = seq_tabix.fetch(chrom, max(0, start - window), end + window)
        CpG, GC = 0, 0
        count = 0
        lbase = ''
        for pos in range(len(sequence)):
            base = sequence[pos]
            count += 1
            if base in 'GC':
                GC += 1
            if lbase == 'C' and base == 'G':
                CpG += 1
            lbase = base
        if count > 0:
            return GC / float(count), CpG / (count * 0.5)
        else:
            return '-', '-'
    except Exception as e:
        print(f"Error fetching sequence for {chrom}:{start}-{end}: {e}")
        return '-', '-'

def extract_transcript_coding_prot_feature(output_dict, annotation, position, label1, label2):
    """
    Extract transcript-related features (e.g., cDNApos, relcDNApos, etc.) from an annotation.
    """
    try:
        # Check if the position exists in the annotation
        if position >= len(annotation) or annotation[position] == "":
            output_dict[label1], output_dict[label2] = ("-", "-")
            return output_dict

        # Check if the value is in the expected "value/length" format
        if "/" not in annotation[position]:
            output_dict[label1], output_dict[label2] = ("-", "-")
            return output_dict

        helper = []
        elength = None
        if (annotation[position] != "-") and (annotation[position][0] != "-"):
            helper = [x.strip() for x in annotation[position].split("/")]
            if len(helper) == 2:
                elength = int(helper[-1])
                annotation[position] = helper[0]
            else:
                sys.exit(
                    'Unexpected format in annotation. Chrom:%s Pos:%s :%s' % (
                        output_dict['#Chrom'], output_dict['Pos'], annotation[position]))
            output_dict[label1] = annotation[position].replace('?-', '').replace('-?', '').split('-')[0]
            if elength is not None:
                output_dict[label2] = "%.2f" % (min(1.0, float(
                    annotation[position].replace('?-', '').replace('-?', '').split('-')[0]) / elength))
        else:
            output_dict[label1], output_dict[label2] = ("-", "-")

        return output_dict
    except ValueError:
        sys.exit("Error processing annotation. Location: %s %s %s" % (
            output_dict['#Chrom'], output_dict['Pos'], annotation[position])) 

def extract_consequences(output_dict, vepfields, fVconseq):
    """
    Extract the consequences from the VEP annotated VCF file and translate them into abbreviations.
    Appends the consequence abbreviation to the given dict and returns it.
    """
    # Define a mapping of consequences to abbreviations
    consequence_map = {
        "stop_gained": "SG",
        "start_lost": "SG",
        "missense_variant": "NS",
        "initiator_codon_variant": "NS",
        "protein_altering_variant": "NS",
        "inframe_insertion": "IF",
        "inframe_deletion": "IF",
        "frameshift_variant": "FS",
        "stop_lost": "SL",
        "incomplete_terminal_codon_variant": "SL",
        "splice_donor_variant": "CS",
        "splice_acceptor_variant": "CS",
        "splice_region_variant": "S",
        "non_coding_exon_variant": "NC",
        "mature_miRNA_variant": "NC",
        "non_coding_transcript_exon_variant": "NC",
        "synonymous_variant": "SN",
        "stop_retained_variant": "SN",
        "intergenic_variant": "IG",
        "intergenic_region": "IG",
        "downstream_gene_variant": "DN",
        "upstream_gene_variant": "UP",
        "feature_truncation": "O",
        "feature_elongation": "O",
        "regulatory_region_variant": "R",
        "TF_binding_site_variant": "R",
        "regulatory_region_amplification": "R",
        "5_prime_UTR_variant": "U5",
        "3_prime_UTR_variant": "U3",
        "intron_variant": "I",
        "coding_sequence_variant": "O",
        "non_coding_transcript_variant": "NC",
        "regulatory_region_ablation": "R",
        "sequence_feature": "SF",
        "conserved_intron_variant": "CI",
        "5_prime_UTR_premature_start_codon_gain_variant": "U5",  # Added this
        "3_prime_UTR_truncation": "U3",  # Added this
        "bidirectional_gene_fusion": "BF",  # Added this
        "non_coding_transcript_exon_variant": "NC",  # Added this
        "splice_polypyrimidine_tract_variant": "SP",  # Added this
    }

    # Extract the consequences from the VEP fields
    consequences = set([x.strip() for x in vepfields[fVconseq].split(",")])

    # Handle combined consequences (e.g., "splice_region_variant&intron_variant")
    for consequence in consequences:
        terms = consequence.split("&")  # Split combined consequences
        for term in terms:
            if term in consequence_map:
                output_dict["Consequence"] = consequence_map[term]
                return output_dict

    # Handle unrecognized consequences
    if len(consequences) == 1:
        sys.stderr.write(f"Unrecognized Consequence: {consequences}\n")
        output_dict["Consequence"] = "O"  # Default to 'UNKNOWN'
    else:
        sys.stderr.write(f"Need simplification: {consequences}\n")
        output_dict["Consequence"] = "O"  # Default to 'UNKNOWN'

    return output_dict

# Function for extracting the changed amino acid from VEP. Appends the AA
# before and after the mutation to the given dict and returns it.
def extract_Aminoacids(output_dict, vepfields, fVAA):
    """
    Extract the original and resulting amino acids from the VEP annotation.
    Save the original amino acid as oAA and the resulting amino acid as nAA.
    Use '*' for stop codons.
    """
    # Mapping of three-letter amino acid codes to one-letter codes
    aa_three_to_one = {
        "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
        "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
        "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
        "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
        "Ter": "*", "Stop": "*"
    }

    # Debug: Print the entire vepfields and the HGVS.p field
    print(f"DEBUG: vepfields = {vepfields}")
    if fVAA >= len(vepfields) or not vepfields[fVAA]:  # Check if HGVS.p field is missing or empty
        print(f"DEBUG: Missing or empty HGVS.p field: {vepfields}")
        output_dict['oAA'] = "-"  # Assign "-" for missing or irrelevant HGVS.p field
        output_dict['nAA'] = "-"
        return output_dict

    try:
        # Extract the HGVS.p field (e.g., p.Leu245Arg)
        hgvs_p = vepfields[fVAA]
        print(f"DEBUG: vepfields[fVAA] = {hgvs_p}")  # Debug: Print the HGVS.p field

        # Extract the original and resulting amino acids
        if "p." in hgvs_p:
            aa_change = hgvs_p.split("p.")[-1]  # Get the part after "p."
            original_aa = ''.join(filter(str.isalpha, aa_change[:3]))  # First three letters are the original AA
            new_aa = ''.join(filter(str.isalpha, aa_change[3:]))  # Remaining letters are the new AA

            # Map the original and new amino acids to one-letter codes
            output_dict['oAA'] = aa_three_to_one.get(original_aa, "-")  # Map original AA
            output_dict['nAA'] = aa_three_to_one.get(new_aa, "-")  # Map new AA
        else:
            print(f"DEBUG: HGVS.p field does not contain 'p.': {hgvs_p}")
            output_dict['oAA'] = "-"
            output_dict['nAA'] = "-"
    except IndexError:
        print(f"DEBUG: IndexError for vepfields[fVAA] = {vepfields[fVAA]}")  # Debug: Handle malformed input
        output_dict['oAA'] = "-"
        output_dict['nAA'] = "-"  # Handle malformed input

    return output_dict

##########################################
############### MAIN SCRIPT ##############
##########################################

# Open the reference genome (FASTA file).
ref_fasta = pysam.FastaFile(args.reference)

# Open the input and output files.
with open(args.input, "r") as infile, open(args.output, "w") as outfile:
    reader = (line for line in infile if not line.startswith("##"))  # Skip metadata lines.
    header = next(reader).strip().split("\t")  # Read the header line.
    writer = csv.DictWriter(outfile, fieldnames=ELIST, delimiter="\t")
    
    # Write the header to the output file.
    writer.writeheader()
    
   # Process each row in the input file.
    for line in reader:
        row = dict(zip(header, line.strip().split("\t")))
        output_dict = {}
        output_dict['#Chrom'] = row['#CHROM']
        output_dict['Pos'] = int(row['POS'])
        output_dict['Ref'] = row['REF']
        output_dict['Alt'] = row['ALT']
        output_dict['isTv'] = is_transversion(row['REF'], row['ALT'])
        
        # Calculate GC% and CpG% using the `count_GC_CpG` function.
        output_dict['GC'], output_dict['CpG'] = count_GC_CpG(
            output_dict['#Chrom'], output_dict['Pos'], output_dict['Pos'], 75, ref_fasta
        )
        
        # Extract transcript-related features and consequences from the `ANN` field.
        for annotation in row['INFO'].split(";"):
            if annotation.startswith("ANN="):
                annotations = annotation[4:].split(",")  # Extract annotations after "ANN="
                for ann in annotations:
                    vepfields = ann.split("|")  # Split annotation into fields
                    
                    # Extract transcript-related features
                    output_dict = extract_transcript_coding_prot_feature(output_dict, vepfields, 10, 'cDNApos', 'relcDNApos')
                    output_dict = extract_transcript_coding_prot_feature(output_dict, vepfields, 11, 'CDSpos', 'relCDSpos')
                    output_dict = extract_transcript_coding_prot_feature(output_dict, vepfields, 12, 'protPos', 'relprotPos')
                    
                    # Extract consequences
                    output_dict = extract_consequences(output_dict, vepfields, 1)  # Use the correct index for consequences
                    
                    # Extract amino acids
                    output_dict = extract_Aminoacids(output_dict, vepfields, 10)  # Assuming index 10 for HGVS.p
                    
                    break  # Process only the first annotation for simplicity
        
        # Fill in placeholders for other fields (e.g., motifECount, Domain, etc.).
        output_dict.update({key: "-" for key in ELIST if key not in output_dict})
        
        # Write the formatted row to the output file.
        writer.writerow(output_dict)