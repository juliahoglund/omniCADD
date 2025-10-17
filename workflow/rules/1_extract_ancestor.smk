# -*- snakemake -*-

'''
 The snakemake file goes through part 1 of extracting an ancestral sequence from a multiple alignment. 
 This file can/will be called from the snakefile, but can be run separately as shown below.
 The scripts directory contains all the used scripts by the snakemake file.
 This pipeline takes in maf files. 
 If the user has another file format, these should be converted before, 
 either with the emf2maf.pl script of with the msaconverter. 

 :Author: Job van Schipstal
 :Date: 21-9-2023

 :Exttension and modification: Julia Höglund
 :Date: 1-4-2023
'''

import sys
from snakemake.io import expand, glob_wildcards

"""
 Parse MAF file and removes ambiguous nucleotides from the alignment.
 All 11 ambiguous symbols are converted to N.
 Only needs to be used if directly processing the .maf files in maftools results in errors.
 This is because MAF duplicate finder only supports [actgACTG-Nn].
"""
rule clean_ambiguous:
    input:
        "resources/alignment/{part}.maf.gz"
    params:
        name_ancestor=config["mark_ancestor"]["name_ancestor"],
        name_reference=config["alignment"]["name_species_interest"],
        p_n=config["mark_ancestor"].get("p_n", 0.8),
        p_gap=config["mark_ancestor"].get("p_gap", 0.8)
    conda:
        get_conda_env("alignment")
    resources:
        runtime=lambda wildcards, attempt: min(480, 60 * attempt)
    threads: 2
    output:
        temp("results/alignment/cleaned_maf/{part}.maf.gz")
    script:
        "../scripts/step_1_extract_ancestor/clean_maf.py"


"""
 Identifies the most recent common ancestor between two given species and marks it with an identifier.
 Config input:
    "ancestor", 	what to name the ancestral node of interest (example: Mouse_Rat)
    "sp1_ab",		name of sp1 in the tree 
    				(how it is named in the alignment file tree section)
    "sp2_ab", 		name of sp2 in the tree (the ancestor of sp1 and 2 will be selected)
    "name_sp1", 	name/label of the species of interest 
    				(how it is named in the alignment file alignment section)
"""
rule mark_ancestor:
    input:
        script = workflow.source_path("../scripts/step_1_extract_ancestor/mark_ancestor.py"),  # Fixed path
        maf = lambda wildcards: get_df_input_maf()  # Move maf to input section
    params:
        sp1_ab = config['mark_ancestor']['sp1_tree_ab'],
        sp2_ab = config['mark_ancestor']['sp2_tree_ab'],
        name_sp1 = config['alignment']['name_species_interest'],
        ancestor = config['mark_ancestor']['name_ancestor']
    conda:
        get_conda_env("ancestor")
    log:
        "results/logs/{part}_mark_ancestor_log.txt"
    output:
        temp("results/alignment/marked_ancestor/{part}.maf.gz")
    shell:
        "python3 {input.script}"
        " -i {input.maf}"
        " -o {output}"
        " -a {params.ancestor}"
        " -l {log}"
        " --sp1-label {params.name_sp1}"
        " --sp1-ab {params.sp1_ab}"
        " --sp2-ab {params.sp2_ab}"


def get_df_input_maf():
    """
    Input based on configuration. If ancestor must be marked that rule is input, if not and also no cleaning is needed,
    the source maf file is taken as input instead. If cleaning is needed that rule is added instead.
    Otherwise the input MAF file is required directly, skipping the other two steps and saving some time.
    :return: str, input file
    """

    # Since we only have one alignment now, check if marking is needed
    if 'name_ancestor' in config["mark_ancestor"]:
        return "results/alignment/marked_ancestor/{part}.maf.gz"

    if config["alignment"]["clean_maf"] == "True":
        return "results/alignment/cleaned_maf/{part}.maf.gz"

    return f"{config['alignment']['path']}{{part}}.maf.gz"

"""
 Removes all duplicate sequences and keeps only the one sequence that is the most similar to the block consensus.
 Can be run in a container, as mafTools is python2.7 dependent and can cause version issues.
"""
rule maf_df:
    input:
        lambda wildcards: get_df_input_maf()
    container:
        "docker://juliahoglund/maftools:latest"
    conda:
        get_conda_env("ancestor")
    threads: 2
    output:
        temp("results/alignment/dedup/{part}.maf.lz4")
    shell:
        ensure_dir("{output}") +
        "lz4 -dc {input} | "
        "mafDuplicateFilter --maf /dev/stdin | lz4 -f stdin {output}"

"""
 Reorders species within any alignment block, so that the wanted species are in front.
 (it also removes sequences that are not from species given in the order)
"""
rule maf_ro:
    input:
        "results/alignment/dedup/{part}.maf.lz4"
    params:
        order = config["alignment"]["filter_order"]
    conda:
        get_conda_env("ancestor")
    container:
        "docker://juliahoglund/maftools:latest"
    threads: 2
    output:
        temp("results/alignment/row_ordered/{part}.maf.lz4")
    shell:
        "lz4 -dc {input} | mafRowOrderer --maf /dev/stdin --order {params.order} | lz4 -f stdin {output}"

"""
 Helper function to gather alignment part files so they can be merged for each chromosome.
 Output: list of str, all part files for that prefix
 Exits the program if no files are found since creating a rule with no inputs would break the workflow.
 
 Amount of parts can be dynamic, so gather parts by looking how many are present.
 We are looking at the original input files so we can build the DAG ahead of time,
 otherwise checkpoints would be needed to reevaluate the DAG.
 both emf and maf input files can be checked, based on the config.
"""
def gather_part_files():
    alignment_config = config['alignment']
    input_pattern = f"{alignment_config['path']}{{part}}.{alignment_config['type']}"
    parts = glob_wildcards(input_pattern).part
    parts_filtered = []
    for part in parts:
        if not any(pattern in part for pattern in alignment_config["exclude_patterns"]):
            parts_filtered.append(part)

    # Formulate filenames as output from the previous step
    infiles = expand(
        f"results/alignment/row_ordered/{{part}}.maf.lz4", part=parts_filtered)

    # Handle the case when no files are found
    if len(infiles) == 0:
        # For linting, return a placeholder
        import sys
        if "--lint" in sys.argv:
            print(f"Warning: No alignment parts found in the form {input_pattern} (linting mode)")
            return ["results/alignment/row_ordered/placeholder.maf.lz4"]
        else:
            # For actual runs, exit with error
            sys.exit(f"No alignment parts found in the form {input_pattern}")

    return infiles

"""
 Go through all MAF alignment files and sort the blocks by the chromosome of the species of interest
 lz4 compression is fast, 500Mb/s compression and multi-GB/s decompression for a single modern cpu core.
"""
rule sort_by_chr:
    input:
        maf=gather_part_files(),
        script=workflow.source_path("../scripts/step_1_extract_ancestor/sort_by_chromosome.py")
    params:
        species_name=config["alignment"]["name_species_interest"],
        chromosomes=config["chromosomes"]["karyotype"],
        ancestor=config["mark_ancestor"]["name_ancestor"],
        directory="results/alignment/merged/"
    conda:
        get_conda_env("alignment")
    resources:
        mem_mb=lambda wildcards, attempt: min(64000, 8000 * attempt),
        runtime=lambda wildcards, attempt: min(720, 120 * attempt),
        tmpdir="results/tmp/sort_chr"
    threads: 8
    output:
        out_chr = expand("results/alignment/merged/chr{chr}.maf",
            chr=config["chromosomes"]["karyotype"])
    shell:
        ensure_dirs(*[f"results/alignment/merged/chr{chr}.maf" for chr in config["chromosomes"]["karyotype"]]) +
        "python3 {input.script} "
        "-s {params.species_name} "
        "-i {input.maf} "
        "-c {params.chromosomes} "
        "-a {params.ancestor} "
        "-o {params.directory}"

"""	
 Flips all alignment blocks in which the species of interest and its ancestors have been on the negative strand. 
"""
rule maf_str:
    input:
        "results/alignment/merged/chr{chr}.maf.gz"
    params:
        species_label = config['alignment']['name_species_interest']
    conda:
        get_conda_env("ancestor")
    container:
        "docker://juliahoglund/maftools:latest"
    threads: 6
    output:
        temp("results/alignment/stranded/chr{chr}.maf.gz")
    shell:
        "gzip -dc {input} | mafStrander --maf /dev/stdin --seq {params.species_label}. --strand + | gzip > {output} && gzip -9 {input}"


"""
 Sorts alignment blocks with respect to coordinates of the first species of interest using its genome.
 Takes input as the fast .lz4 but saves as the more compressed lz4 
 since this final alignment is not marked as temporary.
 If the file was defined to be presorted in the config we skip sorting for a speed benefit.
"""
rule maf_sorter:
    input:
        "results/alignment/stranded/chr{chr}.maf.gz"
    params:
        species_label=config['alignment']['name_species_interest'],
    conda:
        get_conda_env("ancestor")
    container:
        "docker://juliahoglund/maftools:latest"
    threads: 6
    output:
        "results/alignment/sorted/chr{chr}.maf.gz"
    shell:
        "gzip -dc {input} | mafSorter --maf /dev/stdin --seq {params.species_label}. > {output}"

"""
 Reconstructs the marked ancestor sequences in the preprocessed maf files using the identifiers 
 and outputs per chromosome a fasta file of the ancestral sequence. 
"""
rule gen_ancestor_seq:
    input:
        maf=f"results/alignment/sorted/chr{{chr}}.maf.gz",
        script=workflow.source_path("../scripts/step_1_extract_ancestor/extract_ancestor.py")
    params:
        species_name=config["alignment"]["name_species_interest"],
        ancestor=config['mark_ancestor']['name_ancestor'],
        reference=config['mark_ancestor']['reference_genome']
    conda:
        get_conda_env("ancestor")
    output:
        "results/ancestral_seq/{params.ancestor}/chr{chr}.fa"
    shell:
        '''
        faidx -v {params.reference}
        
        python3 {input.script} \
         -i {input.maf} \
         -o {output} \
         -a {params.ancestor} \
         -n {params.species_name} \
         -r {params.reference}.fai \
        '''