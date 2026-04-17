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
        script="workflow/scripts/score_variants/create_variants.py",
    output:
        directory("results/whole_genome_variants/chr{chr}/"),  # Directory output
    log:
        "results/logs/generate_all_variants/chr{chr}.log",
    benchmark:
        "logs/benchmarks/generate_all_variants_chr{chr}.tsv"
    conda:
        get_conda_env("score")
    threads: get_resource("generate_all_variants", "threads")
    resources:
        mem_mb=get_resource("generate_all_variants", "mem_mb"),
        runtime=get_resource("generate_all_variants", "runtime"),
        time=get_resource("generate_all_variants", "time"),
        partition=get_resource("generate_all_variants", "partition"),
    params:
        blocksize=config["parallelization"]["whole_genome_positions_per_file"],
    shell:
        """
         mkdir -p {output} $(dirname {log}) logs/benchmarks
         
         python3 {input.script} -o {output} \
         -s {params.blocksize} -r {input.reference} \
         -c {wildcards.chr} > {log} && \
         for file in {output}/*.vcf
         do
           bgzip "$file" && tabix "$file.gz"
         done
         """


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
        script="workflow/scripts/vep_annotation/vep.sh",
        cache=(rules.vep_cache.output if config["annotation"]["vep"]["cache"]["should_install"] == "True" else []),
    output:
        temp("results/whole_genome_annotations/chr{chr}/{part}_vep_output.tsv"),
    log:
        "results/logs/run_genome_vep/chr{chr}/{part}.log",
    benchmark:
        "logs/benchmarks/run_genome_vep_chr{chr}_{part}.tsv"
    priority: 1
    conda:
        get_conda_env("annotation")
    threads: get_resource("run_genome_vep", "threads")
    resources:
        mem_mb=get_resource("run_genome_vep", "mem_mb"),
        runtime=get_resource("run_genome_vep", "runtime"),
        time=get_resource("run_genome_vep", "time"),
        partition=get_resource("run_genome_vep", "partition"),
    params:
        cache_dir=config["annotation"]["vep"]["cache"]["directory"],
        species_name=config["species_name"],
    shell:
        """
        mkdir -p $(dirname {output}) logs/benchmarks
        mkdir -p $(dirname {log})
        touch {log}
        chmod +x {input.script} >> {log} 2>&1
        {input.script} {input.vcf} {output} {params.cache_dir} {params.species_name} {threads} >> {log} 2>&1
        [[ -s {output} ]] >> {log} 2>&1
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
        grantham=config["annotation"]["grantham_matrix"],
        script="workflow/scripts/vep_annotation/VEP_process.py",
    output:
        temp("results/whole_genome_annotations/chr{chr}/{part}.vep.tsv"),
    log:
        "results/logs/process_genome_vep/chr{chr}/{part}.log",
    benchmark:
        "logs/benchmarks/process_genome_vep_chr{chr}_{part}.tsv"
    priority: 1
    conda:
        get_conda_env("common")
    threads: get_resource("process_genome_vep", "threads")
    resources:
        mem_mb=get_resource("process_genome_vep", "mem_mb"),
        runtime=get_resource("process_genome_vep", "runtime"),
        time=get_resource("process_genome_vep", "time"),
        partition=get_resource("process_genome_vep", "partition"),
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
    output:
        temp("results/whole_genome_annotations/chr{chr}/{part}_annotated.tsv"),
    log:
        "results/logs/intersect_genomewide/chr{chr}/{part}.log",
    benchmark:
        "logs/benchmarks/intersect_bed_chr{chr}_{part}.tsv"
    conda:
        get_conda_env("annotation")
    threads: get_resource("intersect_genomewide", "threads")
    resources:
        mem_mb=get_resource("intersect_genomewide", "mem_mb"),
        runtime=get_resource("intersect_genomewide", "runtime"),
        time=get_resource("intersect_genomewide", "time"),
        partition=get_resource("intersect_genomewide", "partition"),
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
        script="workflow/scripts/combine_annotations/prepare_annotated_data.py",
    output:
        npz="results/dataset/whole_genome_snps/chr{chr}/{part}.npz",
        meta="results/dataset/whole_genome_snps/chr{chr}/{part}.npz.meta.csv.gz",
        cols="results/dataset/whole_genome_snps/chr{chr}/{part}.npz.columns.csv",
    log:
        "results/logs/prepare_whole_genome/chr{chr}/{part}.log",
    benchmark:
        "logs/benchmarks/prepare_whole_genome_chr{chr}_{part}.tsv"
    conda:
        get_conda_env("annotation")
    threads: get_resource("prepare_whole_genome", "threads")
    resources:
        mem_mb=get_resource("prepare_whole_genome", "mem_mb"),
        runtime=get_resource("prepare_whole_genome", "runtime"),
        time=get_resource("prepare_whole_genome", "time"),
        partition=get_resource("prepare_whole_genome", "partition"),
    params:
        derived_variants="",  # Fixed empty params
        y="",
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
        script="workflow/scripts/score_variants/model_predict.py",
    output:
        temp("results/whole_genome_scores/raw_parts/{cols}/chr{chr}/{part}.csv"),
    log:
        "results/logs/score_variants/{cols}/chr{chr}/{part}.log",
    benchmark:
        "logs/benchmarks/score_variants_{cols}_chr{chr}_{part}.tsv"
    priority: 1
    conda:
        get_conda_env("score")
    threads: get_resource("score_variants", "threads")
    resources:
        mem_mb=get_resource("score_variants", "mem_mb"),
        runtime=get_resource("score_variants", "runtime"),
        time=get_resource("score_variants", "time"),
        partition=get_resource("score_variants", "partition"),
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
    output:
        "results/whole_genome_scores/RAW_scores_chr{chr}.csv",
    log:
        "results/logs/sort_raw_scores/chr{chr}.log",
    benchmark:
        "logs/benchmarks/sort_raw_scores_chr{chr}.tsv"
    priority: 1
    conda:
        get_conda_env("common")
    threads: get_resource("sort_raw_scores", "threads")
    resources:
        mem_mb=lambda wildcards, attempt: min(128000, get_resource("sort_raw_scores", "mem_mb") * attempt),
        runtime=lambda wildcards, attempt: get_resource("sort_raw_scores", "runtime") * attempt,
        time=lambda wildcards, attempt: get_resource("sort_raw_scores", "time") * attempt,
        partition=get_resource("sort_raw_scores", "partition"),
        tmpdir="results/tmp/sort_raw_{chr}",
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
    output:
        "results/whole_genome_scores/counts/chr{chr}.txt",
    log:
        "results/logs/count_positions/chr{chr}.log",
    benchmark:
        "logs/benchmarks/count_positions_chr{chr}.tsv"
    priority: 1
    conda:
        get_conda_env("common")
    threads: get_resource("count_positions", "threads")
    resources:
        mem_mb=get_resource("count_positions", "mem_mb"),
        runtime=get_resource("count_positions", "runtime"),
        time=get_resource("count_positions", "time"),
        partition=get_resource("count_positions", "partition"),
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
        script="workflow/scripts/score_variants/assign_phred_scores.py",
    output:
        expand(
            "results/whole_genome_scores/phred/chr{chr}.tsv",
            chr=config["chromosomes"]["score"],
        ),
    log:
        "results/logs/assign_phred_scores.log",
    benchmark:
        "logs/benchmarks/assign_phred_scores.tsv"
    priority: 1
    conda:
        get_conda_env("score")
    threads: get_resource("assign_phred_scores", "threads")
    resources:
        mem_mb=lambda wildcards, attempt: get_resource("assign_phred_scores", "mem_mb") * attempt,
        runtime=lambda wildcards, attempt: get_resource("assign_phred_scores", "runtime") * attempt,
        time=lambda wildcards, attempt: get_resource("assign_phred_scores", "time") * attempt,
        partition=get_resource("assign_phred_scores", "partition"),
    params:
        outmask="results/whole_genome_scores/phred/chrCHROM.tsv",
        chromosomes=config["chromosomes"]["score"],
    shell:
        """
        mkdir -p results/whole_genome_scores/phred logs/benchmarks
        
        python3 {input.script} \\
        -i {input.data} \\
        -o {params.outmask} \\
        --chroms {params.chromosomes} \\
        --count-file {input.counts} 2> {log}
        """


rule merge_all_raw_scores:
    """Merge all per-chromosome sorted score files into a single genome-wide file."""
    input:
        expand(
            "results/whole_genome_scores/RAW_scores_chr{chr}.csv",
            chr=config["chromosomes"]["score"],
        ),
    output:
        "results/whole_genome_scores/full_RAW_scores.csv",
    log:
        "results/logs/merge_all_raw_scores.log",
    benchmark:
        "logs/benchmarks/merge_all_raw_scores.tsv"
    priority: 1
    conda:
        get_conda_env("common")
    threads: get_resource("merge_all_raw_scores", "threads")
    resources:
        mem_mb=lambda wildcards, attempt: get_resource("merge_all_raw_scores", "mem_mb") * attempt,
        runtime=get_resource("merge_all_raw_scores", "runtime"),
        time=get_resource("merge_all_raw_scores", "time"),
        partition=get_resource("merge_all_raw_scores", "partition"),
        tmpdir="results/tmp/sort_merge",
    shell:
        """
        mkdir -p $(dirname {output}) {resources.tmpdir} logs/benchmarks
        LC_ALL=C sort --merge \
        -t "," \
        -k5gr \
        -S {resources.mem_mb}M \
        --temporary-directory={resources.tmpdir} \
        {input} > {output} 2> {log}
        """


rule format_cadd_scores:
    """Sort PHRED scores by position and bgzip for tabix indexing."""
    input:
        "results/whole_genome_scores/phred/chr{chr}.tsv",
    output:
        "results/cadd_scores/chr{chr}.tsv.gz",
    log:
        "results/logs/format_cadd_scores/chr{chr}.log",
    benchmark:
        "logs/benchmarks/format_cadd_scores_chr{chr}.tsv"
    priority: 1
    conda:
        get_conda_env("common")
    threads: get_resource("format_cadd_scores", "threads")
    resources:
        mem_mb=get_resource("format_cadd_scores", "mem_mb"),
        runtime=get_resource("format_cadd_scores", "runtime"),
        time=get_resource("format_cadd_scores", "time"),
        partition=get_resource("format_cadd_scores", "partition"),
    shell:
        """
        mkdir -p $(dirname {output}) logs/benchmarks
        (grep "^#" {input}; grep -v "^#" {input} | sort -k2,2n) | bgzip > {output} 2> {log}
        """


rule index_cadd_scores:
    input:
        "results/cadd_scores/chr{chr}.tsv.gz",
    output:
        "results/cadd_scores/chr{chr}.tsv.gz.tbi",
    log:
        "results/logs/index_cadd_scores/chr{chr}.log",
    priority: 1
    conda:
        get_conda_env("common")
    threads: get_resource("index_cadd_scores", "threads")
    resources:
        mem_mb=get_resource("index_cadd_scores", "mem_mb"),
        runtime=get_resource("index_cadd_scores", "runtime"),
        time=get_resource("index_cadd_scores", "time"),
        partition=get_resource("index_cadd_scores", "partition"),
    shell:
        "tabix -s1 -b2 -e2 {input} 2> {log}"


rule summarize_cadd_scores:
    input:
        scores=expand("results/cadd_scores/chr{chr}.tsv.gz", chr=config["chromosomes"]["score"]),
        indices=expand(
            "results/cadd_scores/chr{chr}.tsv.gz.tbi",
            chr=config["chromosomes"]["score"],
        ),
    output:
        summary=report("results/cadd_scores/scoring_summary.txt", category="Logs"),
    log:
        "results/logs/summarize_cadd_scores.log",
    priority: 1
    conda:
        get_conda_env("common")
    threads: get_resource("summarize_cadd_scores", "threads")
    resources:
        mem_mb=get_resource("summarize_cadd_scores", "mem_mb"),
        runtime=get_resource("summarize_cadd_scores", "runtime"),
        time=get_resource("summarize_cadd_scores", "time"),
        partition=get_resource("summarize_cadd_scores", "partition"),
    shell:
        """
        mkdir -p $(dirname {output.summary})
        echo "CADD scoring completed successfully" > {output.summary} 2>> {log}
        echo "Files created:" >> {output.summary} 2>> {log}
        ls -lh {input.scores} >> {output.summary} 2>> {log}
        echo "Total variants scored:" >> {output.summary} 2>> {log}
        zcat {input.scores} | wc -l >> {output.summary} 2>> {log}
        """
