"""
 Module to generate whole genome PHRED-like CADD scores.
 First all possible variants are generated, these are then annotated and scored
 using the annotate_variants and train_test_model modules of the workflow.
 Afterwards the raw model scores are sorted in descending order and are
 assigned the final PHRED-like CADD score. These are then again sorted by
 chromosome and position. Finally the scores for each position,
 are summarised by a mean, min and max value per position (not yet implemnted).

 :Author: Job van Schipstal
 :Date: 01-11-2023
 :Extension and modification: Julia Höglund
 :Datwe: 13-08-2024

 The scripts and workflow have been adopted from the work of Christian Gross.
"""


wildcard_constraints:
    part="[0-9]+",


CHROMSOME_LIST = config["chromosomes"]["score"]


"""
Generates all possible variants for a chromosome in blocks.
The files are also directly bgzipped.
"""


checkpoint generate_all_variants:
    input:
        reference=config["generate_variants"]["reference_genome_wildcard"],
        script=workflow.source_path(f"{SCRIPTS_8}create_variants.py"),
    params:
        blocksize=config["parallelization"]["whole_genome_positions_per_file"],
    conda:
        get_conda_env("score")
    resources:
        mem_mb=8000,
        runtime=240,  # 4 hours for variant generation
    log:
        "results/whole_genome_variants/chr{chr}/generation.log",
    output:
        directory("results/whole_genome_variants/chr{chr}/"),  # Directory output
    benchmark:
        "logs/benchmarks/generate_all_variants_chr{chr}.tsv"
    shell:
        """
         mkdir -p {output} logs/benchmarks
         
         python3 {input.script} -o {output} \
         -s {params.blocksize} -r {input.reference} \
         -c {wildcards.chr} > {log} && \
         for file in {output}/*.vcf
         do
           bgzip "$file" && tabix "$file.gz"
         done
         """


# hardcoded to All columns as of now
def gather_scores(wildcards):
    """Gather all score files from generate_all_variants checkpoint"""
    checkpoint_output = checkpoints.generate_all_variants.get(**wildcards).output[0]
    parts = glob_wildcards(os.path.join(checkpoint_output, "{part}.vcf.gz")).part
    parts_sorted = natsorted(parts, alg=ns.INT)  # Natural sort
    return expand(
        "results/whole_genome_scores/raw_parts/{cols}/chr{chr}/{part}.csv",
        cols=wildcards.cols,
        chr=wildcards.chr,
        part=parts_sorted,
    )


"""
Annotate a vcf file using Ensembl-VEP.
The VEP cache can automatically be downloaded if should_install is True in the config, 
otherwise a path to an existing cache should be given.
An indexed cache is faster than the standard one, so that is what the vep_cache rule provides.
This rule expects SIFT scores to be available but this is not the case for many species,
"""


# TODO make sift a config option
# TODO: move output of problem_out; see why problem and if they can be incorporated later
rule run_genome_vep:
    input:
        vcf="results/whole_genome_variants/chr{chr}/{part}.vcf.gz",
        script=workflow.source_path(f"{SCRIPTS_5}vep.sh"),
        cache=(
            rules.vep_cache.output
            if config["annotation"]["vep"]["cache"]["should_install"] == "True"
            else []
        ),
    params:
        cache_dir=config["annotation"]["vep"]["cache"]["directory"],
        species_name=config["species_name"],
    log:
        "results/logs/score_variants/run_genome_vep/chr{chr}_{part}.log",
    conda:
        get_conda_env("annotation")
    threads: 2
    resources:
        mem_mb=4000,
        runtime=180,
    priority: 1
    benchmark:
        "logs/benchmarks/run_genome_vep_chr{chr}_{part}.tsv"
    output:
        temp("results/whole_genome_annotations/chr{chr}/{part}_vep_output.tsv"),
    shell:
        """
         mkdir -p $(dirname {output}) logs/benchmarks
         
         chmod +x {input.script} 2>> {log} && \\
         {input.script} {input.vcf} {output} \\
         {params.cache_dir} {params.species_name} {threads} 2>> {log} && \\
         [[ -s {output} ]] 2>> {log}
         """


"""
Processes VEP output into the tsv format used by the later steps.
The VEP consequences are summarised and basic annotations are calculated here as well.
"""


rule process_genome_vep:
    input:
        vcf="results/whole_genome_variants/chr{chr}/{part}.vcf.gz",
        vep="results/whole_genome_annotations/chr{chr}/{part}_vep_output.tsv",
        genome=config["generate_variants"]["reference_genome_wildcard"],
        grantham=workflow.source_path(
            "../../resources/grantham_matrix/grantham_table.tsv"
        ),
        script=workflow.source_path(f"{SCRIPTS_5}VEP_process.py"),
    log:
        "results/logs/score_variants/process_genome_vep/chr{chr}_{part}.log",
    conda:
        get_conda_env("common")
    resources:
        mem_mb=6000,
        runtime=120,
    priority: 1
    benchmark:
        "logs/benchmarks/process_genome_vep_chr{chr}_{part}.tsv"
    output:
        temp("results/whole_genome_annotations/chr{chr}/{part}.vep.tsv"),
    shell:
        """
         mkdir -p $(dirname {output}) 2>> {log}
         mkdir -p logs/benchmarks 2>> {log}
         
         python3 {input.script} -v {input.vep} -s {input.vcf} \
         -r {input.genome} -g {input.grantham} -o {output} --multiple 2>> {log}
         """


rule intersect_genomewide:
    input:
        vep="results/whole_genome_annotations/chr{chr}/{part}.vep.tsv",
        bed="results/annotation/constraint/constraint_chr{chr}.bed",
        script="workflow/scripts/combine_annotations/merge_annotations.py",
    log:
        "results/logs/score_variants/intersect_genomewide/chr{chr}_{part}.log",
    conda:
        get_conda_env("annotation")
    threads: 8
    resources:
        mem_mb=12000,
        runtime=90,
    benchmark:
        "logs/benchmarks/intersect_bed_chr{chr}_{part}.tsv"
    output:
        temp("results/whole_genome_annotations/chr{chr}/{part}_annotated.tsv"),
    shell:
        """
        mkdir -p $(dirname {output}) 2>> {log}
        mkdir -p logs/benchmarks 2>> {log}
        
        python3 {input.script} \
        -v {input.vep} \
        -b {input.bed} \
        -o {output} 2>> {log}
        """


rule prepare_whole_genome:
    input:
        data="results/whole_genome_annotations/chr{chr}/{part}_annotated.tsv",
        imputation="results/dataset/imputation_dict.txt",
        processing=config["annotation_config"]["processing"],
        interactions=config["annotation_config"]["interactions"],
        script=workflow.source_path(f"{SCRIPTS_6}prepare_annotated_data.py"),
    params:
        derived_variants="",  # Fixed empty params
        y="",
    threads: 5
    resources:
        mem_mb=8000,
        runtime=60,
    conda:
        get_conda_env("annotation")
    benchmark:
        "logs/benchmarks/prepare_whole_genome_chr{chr}_{part}.tsv"
    output:
        npz="results/dataset/whole_genome_snps/chr{chr}/{part}.npz",
        meta="results/dataset/whole_genome_snps/chr{chr}/{part}.npz.meta.csv.gz",
        cols="results/dataset/whole_genome_snps/chr{chr}/{part}.npz.columns.csv",
    log:
        "results/logs/data_preparation/whole_genome/chr{chr}/{part}.log",
    shell:
        """
        mkdir -p $(dirname {output.npz})
        mkdir -p $(dirname {log})
        mkdir -p logs/benchmarks
        
        python3 {input.script} -i {input.data} --npz {output.npz} \
        --processing-config {input.processing} \
        --interaction-config {input.interactions} \
        --imputation-dict {input.imputation} \
        {params.derived_variants} {params.y} > {log}
        """


"""
Scores the predicted probability for all possible variants to be of class 1,
(proxy) deleterious. Saved as an csv with chr, pos, ref, alt and raw score.
"""


rule score_variants:
    input:
        data="results/dataset/whole_genome_snps/chr{chr}/{part}.npz",
        data_m="results/dataset/whole_genome_snps/chr{chr}/{part}.npz.meta.csv.gz",
        data_c="results/dataset/whole_genome_snps/chr{chr}/{part}.npz.columns.csv",
        scaler="results/model/{cols}/full.scaler.pickle",
        model="results/model/{cols}/full.mod.pickle",
        script=workflow.source_path(f"{SCRIPTS_8}model_predict.py"),
    log:
        "results/logs/score_variants/score_variants/{cols}_chr{chr}_{part}.log",
    conda:
        get_conda_env("score")
    threads: 4
    resources:
        mem_mb=6000,
        runtime=120,
    benchmark:
        "logs/benchmarks/score_variants_{cols}_chr{chr}_{part}.tsv"
    output:
        temp("results/whole_genome_scores/raw_parts/{cols}/chr{chr}/{part}.csv"),
    shell:
        """
        mkdir -p $(dirname {output}) logs/benchmarks
        
        python3 {input.script} \\
        -i {input.data} \\
        --model {input.model} \\
        --scaler {input.scaler} \\
        -o {output} \\
        --sort \\
        --no-header 2> {log}
        """


"""
sorts all raw scores from highest to lowest, i.e. ranking them later for the
PHRED score generation
"""


rule sort_raw_scores:
    input:
        gather_scores,
    log:
        "results/logs/score_variants/sort_raw_scores/chr{chr}.log",
    threads: 8
    resources:
        mem_mb=lambda wildcards, attempt: min(
            128000, config["memory"]["dataset_mb"] * attempt
        ),
        runtime=lambda wildcards, attempt: min(720, 90 * attempt),
        tmpdir="results/tmp/sort_raw_{chr}",
    benchmark:
        "logs/benchmarks/sort_raw_scores_chr{chr}.tsv"
    output:
        "results/whole_genome_scores/RAW_scores_chr{chr}.csv",
    conda:
        get_conda_env("common")
    shell:
        """
        mkdir -p $(dirname {output}) {resources.tmpdir} logs/benchmarks
        
        LC_ALL=C sort \
        --merge \
        -t "," \
        -k5gr \
        -S {resources.mem_mb}M \
        --parallel={threads} \
        --temporary-directory={resources.tmpdir} \
         {input} > {output} 2> {log}
        """


rule count_positions:
    input:
        "results/whole_genome_scores/RAW_scores_chr{chr}.csv",
    log:
        "results/logs/score_variants/count_positions/chr{chr}.log",
    resources:
        mem_mb=1000,
        runtime=15,
    benchmark:
        "logs/benchmarks/count_positions_chr{chr}.tsv"
    output:
        "results/whole_genome_scores/counts/chr{chr}.txt",
    conda:
        get_conda_env("common")
    shell:
        """
        mkdir -p $(dirname {output}) 2>> {log}
        mkdir -p logs/benchmarks 2>> {log}
        
        wc -l {input} > {output} 2>> {log}
        """


"""
Assigns phred scores to all variants, in addition to the raw scores.
The phred scores are based on a genome-wide level; or as many chromosomes
included in the analysis, rather than a chromosome wide ranking
"""


rule assign_phred_scores:
    input:
        data="results/whole_genome_scores/full_RAW_scores.csv",
        counts=expand(
            "results/whole_genome_scores/counts/chr{chr}.txt",
            chr=config["chromosomes"]["score"],
        ),
        script=workflow.source_path(f"{SCRIPTS_8}assign_phred_scores.py"),
    params:
        outmask="results/whole_genome_scores/phred/chrCHROM.tsv",
        chromosomes=config["chromosomes"]["score"],
    log:
        "results/logs/score_variants/assign_phred_scores.log",
    resources:
        mem_mb=lambda wildcards, attempt: min(192000, 24000 * attempt),  # PHRED assignment is very memory intensive
        runtime=lambda wildcards, attempt: min(960, 120 * attempt),
    benchmark:
        "logs/benchmarks/assign_phred_scores.tsv"
    output:
        expand(
            "results/whole_genome_scores/phred/chr{chr}.tsv",
            chr=config["chromosomes"]["score"],
        ),
    conda:
        get_conda_env("score")
    shell:
        """
        mkdir -p results/whole_genome_scores/phred logs/benchmarks
        
        python3 {input.script} \\
        -i {input.data} \\
        -o {params.outmask} \\
        --chroms {params.chromosomes} \\
        --count-file {input.counts} 2> {log}
        """


rule sort_and_merge_scores:
    input:
        gather_scores,
    log:
        "results/logs/score_variants/sort_and_merge_scores.log",
    threads: 16
    resources:
        mem_mb=config["memory"]["dataset_mb"],
        runtime=240,
        tmpdir="results/tmp/sort_merge",
    output:
        "results/whole_genome_scores/full_RAW_scores.csv.gz",  # Direct to final compressed
    conda:
        get_conda_env("common")
    shell:
        """
        mkdir -p {resources.tmpdir} 2>> {log}
        
        # Pipe directly from sort to compression
        LC_ALL=C sort \
        --merge \
        -t "," \
        -k5gr \
        -S {resources.mem_mb}M \
        --parallel={threads} \
        --temporary-directory={resources.tmpdir} \
         {input} 2>> {log} | gzip > {output}
        """


rule index_cadd_scores:
    input:
        "results/cadd_scores/chr{chr}.tsv.gz",
    log:
        "results/logs/score_variants/index_cadd_scores/chr{chr}.log",
    conda:
        get_conda_env("common")
    resources:
        mem_mb=2000,
        runtime=30,
    output:
        "results/cadd_scores/chr{chr}.tsv.gz.tbi",
    shell:
        "tabix -s1 -b2 -e2 {input} 2> {log}"


rule summarize_cadd_scores:
    input:
        scores=expand(
            "results/cadd_scores/chr{chr}.tsv.gz", chr=config["chromosomes"]["score"]
        ),
        indices=expand(
            "results/cadd_scores/chr{chr}.tsv.gz.tbi",
            chr=config["chromosomes"]["score"],
        ),
    log:
        "results/logs/score_variants/summarize_cadd_scores.log",
    conda:
        get_conda_env("common")
    resources:
        mem_mb=4000,
        runtime=60,
    output:
        summary=report("results/cadd_scores/scoring_summary.txt", category="Logs"),
    shell:
        """
        mkdir -p $(dirname {output.summary})
        echo "CADD scoring completed successfully" > {output.summary} 2>> {log}
        echo "Files created:" >> {output.summary} 2>> {log}
        ls -lh {input.scores} >> {output.summary} 2>> {log}
        echo "Total variants scored:" >> {output.summary} 2>> {log}
        zcat {input.scores} | wc -l >> {output.summary} 2>> {log}
        """
