"""
Basic helper rules for the other parts of the workflow.
Configuration loading helper functions are also found here.

:Author: Job van Schipstal
:Date: 21-9-2023
:extension and modification: Julia Höglund
:Date: 13-8-2025
"""

import subprocess
import os
import sys
import glob


# Global conda environment resolver
def get_conda_env(env_name):
    """Return path to conda environment file"""
    import os

    return os.path.join(workflow.basedir, "workflow", "envs", f"{env_name}.yml")


def get_alignment_config():
    """Get the correct alignment configuration based on the configured name (mark_ancestor: ancestral_alignment)."""
    # Try to get the name from config, falling back to any available or 43_mammals.epo
    alignment_name = config.get("mark_ancestor", {}).get("ancestral_alignment")

    if not alignment_name:
        alignments = config.get("ancestral_sequence", {}).get("alignment", {}).get("alignments", {})
        if alignments:
            alignment_name = next(iter(alignments))
        else:
            alignment_name = "43_mammals.epo"

    return config["ancestral_sequence"]["alignment"]["alignments"][alignment_name]


def should_skip_ancestral_path():
    """
    Check if ancestral and conservation paths use same species order.
    If true, the workflow can skip ancestral-specific MAF processing and
    reuse conservation outputs, saving ~50% of alignment processing time.
    """
    filter_order = get_alignment_config().get("filter_order", "")
    conservation_order = get_alignment_config().get("conservation_species_order", "")
    return filter_order == conservation_order and conservation_order != ""


def get_maf_ro_output(wildcards):
    """
    Return appropriate maf_ro output depending on whether paths are merged.
    When species orders match, returns conservation path to avoid duplication.
    """
    if should_skip_ancestral_path():
        return f"results/alignment/row_ordered_conservation/{wildcards.part}.maf.lz4"
    else:
        return f"results/alignment/row_ordered/{wildcards.part}.maf.lz4"


def get_sort_by_chr_output(wildcards):
    """
    Return appropriate sorted MAF output depending on whether paths are merged.
    When species orders match, returns conservation path to avoid duplication.
    """
    if should_skip_ancestral_path():
        return f"results/alignment/merged_conservation/chr{wildcards.chr}.maf"
    else:
        return f"results/alignment/merged/chr{wildcards.chr}.maf.gz"


def get_gene_annotation_file():
    """
    Return appropriate gene annotation file based on config.
    Used by: gene_prediction.smk
    """
    source = config["gene_annotation"]["source"]

    if source == "gff":
        return config["gene_annotation"]["gff"]
    elif source == "gtf":
        return config["gene_annotation"]["gtf"]
    elif source == "augustus":
        return "results/gene_prediction/genes_validated.gff3"
    elif source == "none":
        return None
    else:
        raise ValueError(f"Unknown gene annotation source: {source}")


def get_df_input_maf(wildcards):
    """
    Input for maf_df rule (duplicate filter step), which follows after marking.
    Pipeline: clean_ambiguous → mark_ancestor → maf_df → maf_ro
    If mark_ancestor is configured, maf_df takes marked_ancestor output.
    If only cleaning is configured, maf_df takes cleaned_maf output.
    Otherwise maf_df takes raw alignment files directly.
    Used by: ancestral_generation.smk (maf_df input)
    :param wildcards: Snakemake wildcards object
    :return: str, input file path with wildcard resolved
    """
    alignment_config = get_alignment_config()

    if "name_ancestor" in config["mark_ancestor"]:
        return f"results/alignment/marked_ancestor/{wildcards.part}.maf.gz"

    if alignment_config["clean_maf"]:
        return f"results/alignment/cleaned_maf/{wildcards.part}.maf.gz"

    return f"{alignment_config['path']}{wildcards.part}.maf.gz"


def get_mark_ancestor_input_maf(wildcards):
    """
    Input for mark_ancestor rule (the step before maf_df).
    Pipeline: clean_ambiguous → mark_ancestor
    If cleaning is configured, mark_ancestor takes cleaned_maf as input.
    Otherwise mark_ancestor takes raw alignment files directly.
    Used by: ancestral_generation.smk (mark_ancestor input)
    :param wildcards: Snakemake wildcards object
    :return: str, input file path with wildcard resolved
    """
    alignment_config = get_alignment_config()

    if alignment_config["clean_maf"]:
        return f"results/alignment/cleaned_maf/{wildcards.part}.maf.gz"

    return f"{alignment_config['path']}{wildcards.part}.maf.gz"


def get_conservation_species_order():
    """
    Get species order for conservation files, ensuring ancestor is included.
    Conservation files need the ancestor for ancestral extraction when paths are shared.
    Returns: Space-separated string of species names including the ancestor.
    """
    alignment_config = get_alignment_config()
    ancestor_name = config["mark_ancestor"].get("name_ancestor", "")
    
    # Get the conservation order, defaulting to filter_order
    order = alignment_config.get("conservation_species_order", alignment_config["filter_order"])
    
    # If ancestor is configured and not already in the order, add it at the beginning
    if ancestor_name and ancestor_name not in order:
        # Assume order is space-separated string
        order_parts = order.split()
        # Insert ancestor after the reference species (first in list)
        if len(order_parts) > 0:
            order_parts.insert(1, ancestor_name)
            order = " ".join(order_parts)
    
    return order


def gather_part_files():
    """
    Gather all alignment part files based on configuration.
    Used by: ancestral_generation.smk
    
    If species orders are identical (should_skip_ancestral_path returns True),
    returns conservation files directly to avoid duplicate processing.
    """
    alignment_config = get_alignment_config()
    input_pattern = f"{alignment_config['path']}{{part}}.{alignment_config['type']}"
    parts = glob_wildcards(input_pattern).part

    parts_filtered = []
    for part in parts:
        if not any(pattern in part for pattern in alignment_config["exclude_patterns"]):
            parts_filtered.append(part)

    # If species orders are identical, use conservation outputs directly
    if should_skip_ancestral_path():
        infiles = expand(f"results/alignment/row_ordered_conservation/{{part}}.maf.lz4", part=parts_filtered)
    else:
        # Normal path: use row_ordered outputs
        infiles = expand(f"results/alignment/row_ordered/{{part}}.maf.lz4", part=parts_filtered)

    # Handle the case when no files are found
    if len(infiles) == 0:
        # For linting or dry-run, return a placeholder
        import sys

        if "--lint" in sys.argv or "--dry-run" in sys.argv or "-n" in sys.argv:
            mode = "lint/dry-run" if any(arg in sys.argv for arg in ["--lint", "--dry-run", "-n"]) else "lint"
            print(f"Warning: No alignment parts found in the form {input_pattern} ({mode} mode)")
            return ["results/alignment/row_ordered/placeholder.maf.lz4"]
        else:
            # For actual runs, exit with error
            sys.exit(f"No alignment parts found in the form {input_pattern}")

    return infiles


def gather_part_files_conservation():
    """
    Gather all alignment part files for conservation analysis (preserves all species).
    Used by: ancestral_generation.smk (sort_by_chr_conservation rule)
    """
    alignment_config = get_alignment_config()
    alignment_dir = alignment_config["path"]
    parts = glob_wildcards(os.path.join(alignment_dir, "{part}.maf.gz")).part

    # Handle the case when no files are found
    if len(parts) == 0:
        import sys

        if "--lint" in sys.argv or "--dry-run" in sys.argv or "-n" in sys.argv:
            print(f"Warning: No conservation alignment parts found in {alignment_dir} (dry-run mode)")
            return ["results/alignment/row_ordered_conservation/placeholder.maf.lz4"]
        else:
            sys.exit(f"No conservation alignment parts found in {alignment_dir}")

    return expand("results/alignment/row_ordered_conservation/{part}.maf.lz4", part=parts)


def get_folds(excluding=None) -> list:
    """
    Get list of numbers, one for each fold that is to be taken as input.
    Used by: train_test_model.smk
    :param excluding: optional int(-like), fold to exclude (def None)
    :return: List of numbers, useful for snakemake.expand
    """
    folds = list(range(config["model"]["n_folds"]))
    if excluding:
        excluding = int(excluding)
        if excluding in folds:
            folds.remove(excluding)
    return folds


def get_train_folds(fold):
    """
    Get all folds except the test fold.
    Used by: train_test_model.smk
    """
    all_folds = list(range(config["model"]["n_folds"]))
    return [f for f in all_folds if f != int(fold)]


def gather_scores(wildcards):
    """
    Gather all score files from generate_all_variants checkpoint.
    Used by: score_variants.smk
    cols defaults to "All" when not present in wildcards (e.g. sort_raw_scores).
    Parts are always numeric (wildcard_constraints part="[0-9]+"), so sorted(key=int) suffices.
    """
    cols = getattr(wildcards, "cols", "All")
    checkpoint_output = checkpoints.generate_all_variants.get(chr=wildcards.chr).output[0]
    parts = glob_wildcards(os.path.join(checkpoint_output, "{part}.vcf.gz")).part
    parts_sorted = sorted(parts, key=int)
    return expand(
        "results/whole_genome_scores/raw_parts/{cols}/chr{chr}/{part}.csv",
        cols=cols,
        chr=wildcards.chr,
        part=parts_sorted,
    )


def get_alignment_parts(wildcards):
    """Resolve chunk part IDs for a chromosome using the split_alignment checkpoint."""
    cp = checkpoints.split_alignment.get(chr=wildcards.chr)
    (parts,) = glob_wildcards(os.path.join(cp.output[0], f"chr{wildcards.chr}-{{part}}.maf"))
    return parts


"""
Global wildcard constraints, ease matching of wildcards in rules.
Chr is constrained to only be numbers or letters.
Name, label and file may not contain /, they may not be sub-folders.
"""


# Global wildcard constraints for use across modules
wildcard_constraints:
    chr="[0-9XY]+",
    part="[a-zA-Z0-9._-]+",  # Allow dots, underscores, hyphens for MAF part names
    fold="[0-9]+",
    type="(simulated|derived|validation)",
    cols="(All|[A-Za-z0-9_]+)",
    ancestor="[A-Za-z0-9_]+",
    c="[0-9.]+",
    iter="[0-9]+",
    tool="(phastCons|phyloP)",


"""
Bgzip_tabix combines bgzip and the tabix rule to reduce overhead.
This prioritises the combined rule over just tabix,
desired since they produce the same output.
Bgzip_validation_variants has the highest priority since
it also handles moving the variants into the results folder.
"""


ruleorder: bgzip_tabix > tabix


"""
Unzip MAF files since the tool needs a seekable file,
it is more effective to not have to decompress multiple times.
Marked as temp so automatically deleted when longer needed, leaving only the compressed original.
"""


rule unzip_maf:
    input:
        "{folder}/{file}.maf.gz",
    output:
        "{folder}/{file}.maf",  # TODO Mark TEMP after testing
    log:
        "results/logs/unzip_maf/{folder}/{file}.log",
    conda:
        get_conda_env("common")
    threads: get_resource("unzip_maf", "threads")
    resources:
        mem_mb=get_resource("unzip_maf", "mem_mb"),
        runtime=get_resource("unzip_maf", "runtime"),
        time=get_resource("unzip_maf", "time"),
        partition=get_resource("unzip_maf", "partition"),
    shell:
        "gzip -dc {input} > {output} 2> {log}"


"""
Counts variants in a VCF file, by counting all lines not starting with a comment or whitespace
"""


rule count_variants:
    input:
        "{folder}/{file}.vcf",
    output:
        "{folder}/{file}.vcf.count",
    log:
        "results/logs/count_variants/{folder}/{file}.log",
    conda:
        get_conda_env("common")
    threads: get_resource("count_variants", "threads")
    resources:
        mem_mb=get_resource("count_variants", "mem_mb"),
        runtime=get_resource("count_variants", "runtime"),
        time=get_resource("count_variants", "time"),
        partition=get_resource("count_variants", "partition"),
    shell:
        r"grep -c '^[^#\S]' {input} > {output} 2> {log}"


"""
Sums counts of all per-chromosome VCF files in a folder.
Since this is only used for training we sum for the training chromosomes.
"""


rule add_counts:
    input:
        expand("{{folder}}/chr{chr}.vcf.count", chr=config["chromosomes"]["karyotype"]),
    output:
        report("{folder}/total.count", category="Logs"),
    log:
        "results/logs/add_counts/{folder}/total.log",
    conda:
        get_conda_env("common")
    threads: get_resource("add_counts", "threads")
    resources:
        mem_mb=get_resource("add_counts", "mem_mb"),
        runtime=get_resource("add_counts", "runtime"),
        time=get_resource("add_counts", "time"),
        partition=get_resource("add_counts", "partition"),
    shell:
        "cat {input} | awk '{{s+=$1}} END {{print s}}' > {output} 2> {log}"


"""
Compress VCF file using bgzip and index using tabix
"""


rule bgzip_tabix:
    input:
        "{folder}/{file}.vcf",
    output:
        vcf=temp("{folder}/{file}.vcf.gz"),
        index=temp("{folder}/{file}.vcf.gz.tbi"),
    log:
        "results/logs/bgzip_tabix/{folder}/{file}.log",
    conda:
        get_conda_env("common")
    threads: get_resource("bgzip_tabix", "threads")
    resources:
        mem_mb=get_resource("bgzip_tabix", "mem_mb"),
        runtime=get_resource("bgzip_tabix", "runtime"),
        time=get_resource("bgzip_tabix", "time"),
        partition=get_resource("bgzip_tabix", "partition"),
    shell:
        "bgzip -c {input} > {output.vcf} 2> {log} && " "tabix -p vcf -f {output.vcf} 2>> {log}"


"""
Index VCF using tabix
"""


rule tabix:
    input:
        "{folder}/{file}.vcf.gz",
    output:
        temp("{folder}/{file}.vcf.gz.tbi"),
    log:
        "results/logs/tabix/{folder}/{file}.log",
    conda:
        get_conda_env("common")
    threads: get_resource("tabix", "threads")
    resources:
        mem_mb=get_resource("tabix", "mem_mb"),
        runtime=get_resource("tabix", "runtime"),
        time=get_resource("tabix", "time"),
        partition=get_resource("tabix", "partition"),
    shell:
        "tabix -p vcf -f {input} 2> {log}"


def load_tsv_configuration(file: str) -> dict:
    """
    Loads configuration from tsv table file.
    The double commented (##) line is read as column names

    :param file: str, filename to load configuration from
    :return: dict key is entry label and value is dict,
    with key column label and the value of that column
    """
    file_h = open(file, "r")
    elements = []
    samples = {}
    for line in file_h:
        line = line.strip()
        parts = line.split("\t")
        if line.startswith("##"):
            elements = parts[1:]
        if line.startswith("#") or len(line) == 0:
            continue
        samples[parts[0]] = dict([(label, value) for label, value in zip(elements, parts[1:])])

    return samples


def ensure_dir(path):
    """Helper function to ensure directory exists"""
    return f"mkdir -p $(dirname {path}) && "


def ensure_dirs(*paths):
    """Helper function to ensure multiple directories exist"""
    dirs = " ".join([f"$(dirname {path})" for path in paths])
    return f"mkdir -p {dirs} && "


# Add standard directory creation rule
rule create_base_directories:
    output:
        touch("results/.directories_created"),
    log:
        "results/logs/create_base_directories.log",
    conda:
        get_conda_env("common")
    threads: get_resource("create_base_directories", "threads")
    resources:
        mem_mb=get_resource("create_base_directories", "mem_mb"),
        runtime=get_resource("create_base_directories", "runtime"),
        time=get_resource("create_base_directories", "time"),
        partition=get_resource("create_base_directories", "partition"),
    shell:
        """
        mkdir -p results/{{annotation,dataset,logs,temp,model,scores,figures}}/ 2> {log}
        mkdir -p results/annotation/{{vep,constraint,gerp,phast}}/
        mkdir -p results/dataset/{{simulated,derived,validation,whole_genome_snps}}/
        mkdir -p results/whole_genome_{{variants,annotations,scores}}/
        mkdir -p results/cadd_scores/
        mkdir -p logs/benchmarks/
        """


# Add cleanup checkpoints:


checkpoint cleanup_temp_files:
    input:
        # After annotation is complete
        expand(
            "results/annotation/constraint/constraint_chr{chr}.bed",
            chr=config["chromosomes"]["karyotype"],
        ),
    output:
        touch("results/.cleanup_annotation"),
    log:
        "results/logs/cleanup_temp_files.log",
    conda:
        get_conda_env("common")
    threads: get_resource("cleanup_temp_files", "threads")
    resources:
        mem_mb=get_resource("cleanup_temp_files", "mem_mb"),
        runtime=get_resource("cleanup_temp_files", "runtime"),
        time=get_resource("cleanup_temp_files", "time"),
        partition=get_resource("cleanup_temp_files", "partition"),
    shell:
        """
        # Clean up large temporary alignment files
        rm -rf results/alignment/splitted/
        rm -rf results/alignment/fasta/
        rm -rf results/alignment/pruned/
        echo "Cleaned up annotation intermediate files"
        """


def get_folds():
    """Get all fold numbers."""
    return list(range(1, config["model"]["n_folds"] + 1))


def get_train_folds(test_fold):
    """Get all folds except the test fold."""
    all_folds = get_folds()
    return [f for f in all_folds if f != int(test_fold)]


def get_model_columns():
    """Get all model column subsets."""
    return list(config["model"]["column_subsets"].keys()) + ["All"]


# Add config validation
def validate_config():
    required_keys = [
        "chromosomes.karyotype",
        "chromosomes.train",
        "chromosomes.score",
        "species_name",
        "model.n_folds",
        "annotation_config.processing",
        "annotation_config.interactions",
    ]

    for key in required_keys:
        keys = key.split(".")
        value = config
        try:
            for k in keys:
                value = value[k]
        except KeyError:
            raise ValueError(f"Missing required config key: {key}")


# Call validation
validate_config()
