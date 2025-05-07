#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Updated script to process VCF input format with `ANN` annotations while retaining all functionality, including GC% and CpG% calculations.
"""

import csv
from argparse import ArgumentParser
import pysam
import sys, os

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
parser.add_argument("-g", "--grantham",
    help="Path to the Grantham scores file",
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
# Function for reading the Grantham file and returning it as a dictionary.
def read_grantham(filename):
    """
    Read the Grantham scores file and return it as a dictionary.
    """
    grantham = {}
    if os.path.exists(filename):
        with open(filename) as in_h:
            for g_line in in_h:
                fields = g_line.split()
                if len(fields) == 2:
                    aminos = tuple(fields[0].upper().split("-"))
                    grantham[aminos] = fields[1]
                else:
                    sys.stderr.write("Grantham scores, unexpected line, skipping: %s\n" % g_line.strip())
    else:
        sys.stderr.write("Grantham scores input file does not exist: %s\n" % filename)
    return grantham

def extract_alleles_locs(output_dict, fVCoord, fVallele, fVName, vepfields):
    """
    Extract chromosome, position, Ref, and Alt alleles from the new input format.
    """
    output_dict['#Chrom'] = vepfields[fVCoord].split(':')[0]
    output_dict['Pos'] = int(vepfields[fVCoord].split(':')[1])
    output_dict['Alt'] = vepfields[fVallele].upper()
    output_dict['Ref'] = vepfields[fVName].split('/')[0][-1]  # Assuming Ref is in this format
    return output_dict

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

# Function for returning the most deleterious annotation for the same variant,
# when there are two annotations given for a single variant.
def indexing(previous, current):
    global HIERACHY1

    index_current = HIERACHY1.index(current)
    index_previous = HIERACHY1.index(previous)

    if index_previous > index_current:
        return current
    else:
        return previous
    
def extract_Aminoacids(output_dict, vepfields, fVAA, grantham_scores):
    """
    Extract the original and resulting amino acids from the VEP annotation.
    Save the original amino acid as oAA, the resulting amino acid as nAA,
    and the Grantham score for the amino acid change.
    """
    # Mapping of three-letter amino acid codes to one-letter codes
    aa_three_to_one = {
        "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
        "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
        "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
        "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
        "Ter": "*", "Stop": "*"
    }

    # Check if the HGVS.p field is valid
    if fVAA >= len(vepfields) or not vepfields[fVAA]:
        output_dict['oAA'] = "-"
        output_dict['nAA'] = "-"
        output_dict['Grantham'] = "-"
        return output_dict

    # Extract the HGVS.p field (e.g., p.Leu245Arg)
    hgvs_p = vepfields[fVAA]

    # Extract the original and resulting amino acids
    if "p." in hgvs_p:
        aa_change = hgvs_p.split("p.")[-1]
        original_aa = ''.join(filter(str.isalpha, aa_change[:3])).capitalize()
        new_aa = ''.join(filter(str.isalpha, aa_change[3:])).capitalize()

        # Map the original and new amino acids to one-letter codes
        oAA = aa_three_to_one.get(original_aa, "-")
        nAA = aa_three_to_one.get(new_aa, "-")

        output_dict['oAA'] = oAA
        output_dict['nAA'] = nAA

        # Handle synonymous changes (oAA == nAA)
        if oAA == nAA:
            output_dict['Grantham'] = "-"  # No Grantham score for synonymous changes
        elif "*" in (oAA, nAA):  # Handle stop codon changes; as of now also set to missing
            output_dict['Grantham'] = "-"  # Assign a default score for stop codon changes
            
        else:
            # Look up the Grantham score
            grantham_key = (oAA, nAA)
            reverse_grantham_key = (nAA, oAA)
            if grantham_key in grantham_scores:
                output_dict['Grantham'] = grantham_scores[grantham_key]
            elif reverse_grantham_key in grantham_scores:
                output_dict['Grantham'] = grantham_scores[reverse_grantham_key]
            else:
                output_dict['Grantham'] = "-"
    else:
        output_dict['oAA'] = "-"
        output_dict['nAA'] = "-"
        output_dict['Grantham'] = "-"

    return output_dict

def extract_extra(output_dict, vepfields, fVExtra):
    """
    Extract extra information from the VEP annotation and append it to the output dictionary.
    Handles SIFT, DOMAINS, HIGH_INF_POS, and MOTIF_SCORE_CHANGE fields.
    """
    # Check if the fVExtra field exists and is not empty
    if fVExtra >= len(vepfields) or not vepfields[fVExtra]:
        # Assign default values if the field is missing
        output_dict['SIFTcat'] = "-"
        output_dict['SIFTval'] = "-"
        output_dict['Domain'] = "-"
        output_dict['motifEHIPos'] = "-"
        output_dict['motifEScoreChng'] = "0.0"
        output_dict['motifECount'] = "0"
        return output_dict

    # Process the extra field
    for elem in [x.strip() for x in vepfields[fVExtra].split(";")]:
        hfields = [x.strip() for x in elem.split("=")]
        if len(hfields) < 2:
            continue  # Skip malformed entries

        if hfields[0] == "SIFT":
            # Extract SIFT category and value
            hfields2 = [x.strip() for x in hfields[1].rstrip(")").split("(")]
            output_dict['SIFTcat'] = hfields2[0] if len(hfields2) > 0 else "-"
            output_dict['SIFTval'] = hfields2[1] if len(hfields2) > 1 else "-"
        elif hfields[0] == "DOMAINS":
            # Extract domain-related information
            ncoils, tmhmm, sigp = False, False, False
            ndomain, lcompl = False, False
            for dfields in [x.strip() for x in hfields[1].split(',')]:
                if len([x.strip() for x in dfields.split(":")]) < 2:
                    continue

                category, name = [x.strip() for x in dfields.split(":")][0:2]
                if ("_domain" in category) or ("_profile" in category):
                    ndomain = True
                else:
                    name = name.lower()
                    if "coil" in name:
                        ncoils = True
                    elif "tmhelix" in name:
                        tmhmm = True
                    elif "signalp" in name:
                        sigp = True
                    elif name == "seg":
                        lcompl = True

            # Implement simple hierarchy of domain annotations:
            if ncoils:
                output_dict['Domain'] = "ncoils"
            elif tmhmm:
                output_dict['Domain'] = "tmhmm"
            elif sigp:
                output_dict['Domain'] = "sigp"
            elif ndomain:
                output_dict['Domain'] = "ndomain"
            elif lcompl:
                output_dict['Domain'] = "lcompl"
        elif hfields[0] == "HIGH_INF_POS":
            # Extract high information position
            output_dict['motifEHIPos'] = "True" if hfields[1] == "Y" else "False"
        elif hfields[0] == "MOTIF_SCORE_CHANGE":
            # Extract motif score change
            output_dict['motifEScoreChng'] = hfields[1]
            output_dict['motifECount'] = "1"

    # Assign default values for missing fields
    if 'motifEScoreChng' not in output_dict:
        output_dict['motifEScoreChng'] = "0.0"
        output_dict['motifECount'] = "0"
    if 'motifEHIPos' not in output_dict:
        output_dict['motifEHIPos'] = "-"
    if 'Domain' not in output_dict:
        output_dict['Domain'] = "-"
    if 'SIFTval' not in output_dict:
        output_dict['SIFTval'] = "-"
    if 'SIFTcat' not in output_dict:
        output_dict['SIFTcat'] = "-"

    return output_dict


##########################################
############### MAIN SCRIPT ##############
##########################################

# Open the reference genome (FASTA file).
ref_fasta = pysam.FastaFile(args.reference)
# Load the Grantham scores
grantham_scores = read_grantham(args.grantham)

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
                    output_dict = extract_Aminoacids(output_dict, vepfields, 10, grantham_scores)  # Assuming index 10 for HGVS.p
                    
                    break  # Process only the first annotation for simplicity
        
        # Fill in placeholders for other fields (e.g., motifECount, Domain, etc.).
        output_dict.update({key: "-" for key in ELIST if key not in output_dict})
        
        # Write the formatted row to the output file.
        writer.writerow(output_dict)