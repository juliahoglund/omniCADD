from snakemake.utils import min_version
from glob import glob
import os

##### pipeline version ######
omniCADD_version = "0.1.0"

##### set minimum snakemake version #####
min_version("8.0.0")


##### Load config files #####
configfile: "config/config.yaml"


# Load references config
import yaml

with open(config.get("references", "config/config_refs.yaml"), "r") as f:
    references = yaml.safe_load(f)

# Merge references into config for backward compatibility
# but also make them available as separate 'references' dict
config["references"] = references

# Process OMNICADD_REF_DATA variable in references
import re


def expand_ref_paths(d, base_path):
    """Recursively expand {{OMNICADD_REF_DATA}} placeholders"""
    if isinstance(d, dict):
        return {k: expand_ref_paths(v, base_path) for k, v in d.items()}
    elif isinstance(d, list):
        return [expand_ref_paths(item, base_path) for item in d]
    elif isinstance(d, str):
        return re.sub(r"\{\{OMNICADD_REF_DATA\}\}", base_path, d)
    else:
        return d


# Get base path from environment or config
import os

base_ref_path = os.environ.get("OMNICADD_REF_DATA", references.get("OMNICADD_REF_DATA", "resources"))
references = expand_ref_paths(references, base_ref_path)
config["references"] = references


##### setup report #####
report: "report/workflow.rst"


##### PARAMS #####

REFERENCE = "resources/genome/"
SCRIPTS_3 = "workflow/scripts/variant_simulation/"


##### Load configuration logic and helper functions #####
include: "workflow/rules/config_helpers.smk"
##### Load common rules #####
include: "workflow/rules/common.smk"


##### Conditionally load preprocessing rules #####
if should_include_preprocessing():

    include: "workflow/rules/data_preparation.smk"

    print("Including preprocessing rules")


##### Conditionally load gene prediction rules #####
if should_include_gene_prediction():

    include: "workflow/rules/gene_prediction.smk"

    print("Including gene prediction rules")


##### Load core workflow modules #####
# Only include alignment-based ancestor extraction if using pre-existing MAF files
ancestral_source = config.get("ancestral_sequence", {}).get("source", "alignment")
if ancestral_source == "alignment" and not config.get("ancestral_sequence", {}).get("skip_extraction", False):

    include: "workflow/rules/ancestral_generation.smk"

    print("Including alignment-based ancestor extraction rules")


# Skipping ancestor extraction - using pre-computed ancestral sequences

# Conditionally load ancestral reconstruction
if should_include_ancestral_reconstruction():

    include: "workflow/rules/ancestral_reconstruction.smk"

    print("Using ancestral reconstruction workflow")


include: "workflow/rules/variant_derivation.smk"
include: "workflow/rules/variant_simulation.smk"
include: "workflow/rules/summary_report.smk"


##### Conditionally load annotation methods #####
if should_include_vep():

    include: "workflow/rules/vep_annotation.smk"

    print("Including VEP annotation rules")


if should_include_snpeff():

    include: "workflow/rules/snpeff_annotation.smk"

    print("Including SNPEff annotation rules")


include: "workflow/rules/phast_annotation.smk"
include: "workflow/rules/combine_annotations.smk"
include: "workflow/rules/train_test_model.smk"
include: "workflow/rules/score_variants.smk"


##### Workflow-specific target rules #####


rule preprocessing_only:
    """Run only preprocessing steps"""
    input:
        (
            expand("resources/genome/Sus_scrofa_ref_chr{chr}.fa", chr=config["chromosomes"]["karyotype"])
            if should_include_preprocessing()
            else []
        ),
        "results/preprocessing_validation/validation_report.txt" if should_include_preprocessing() else [],


rule annotation_only:
    """Run only annotation steps"""
    input:
        expand(
            "results/annotation/{annotator}/{type}/chr{chr}_{annotator}.tsv",
            annotator=config["annotation"]["variant_annotator"],
            chr=config["chromosomes"]["karyotype"],
            type=["simulated", "derived"],
        ),


rule conservation_only:
    """Run only conservation annotation"""
    input:
        expand("results/annotation/constraint/constraint_chr{chr}.bed", chr=config["chromosomes"]["karyotype"]),


rule model_only:
    """Train model without scoring"""
    input:
        "results/model/All/full.mod.pickle",
        "results/model/All/full.scaler.pickle",
        "results/model/All/full.mod.weights.csv",


rule test_dag:
    """CI target: validates DAG for all core modules (1-8)."""
    input:
        # Module 1: Extract ancestor
        expand(
            "results/ancestral_seq/{ancestor}/chr{chr}.fa",
            ancestor=config["mark_ancestor"]["name_ancestor"],
            chr=config["chromosomes"]["karyotype"],
        ),
        # Module 2: Derive variants
        expand("results/derived_variants/singletons/chr{chr}.vcf", chr=config["chromosomes"]["karyotype"]),
        # Module 3: Simulate variants
        expand("results/simulated_variants/trimmed_snps/chr{chr}.vcf", chr=config["chromosomes"]["karyotype"]),
        # Module 4: Summary report
        "results/visualisation/stats_report.html",
        # Module 5: Annotate variants (VEP or SNPEff depending on config)
        expand(
            "results/annotation/vep/{type}/chr{chr}_vep.tsv",
            type=["simulated", "derived"],
            chr=config["chromosomes"]["karyotype"],
        )
        if should_include_vep()
        else [],
        expand(
            "results/annotation/snpeff/{type}/chr{chr}_snpeff.tsv",
            type=["simulated", "derived"],
            chr=config["chromosomes"]["karyotype"],
        )
        if should_include_snpeff()
        else [],
        # Module 6: Combine variant annotations with conservation constraint
        expand(
            "results/annotation/constraint/constraint_chr{chr}.bed",
            chr=config["chromosomes"]["karyotype"],
        ),
        expand(
            "results/dataset/{type}/chr{chr}_annotated.tsv",
            type=["simulated", "derived"],
            chr=config["chromosomes"]["karyotype"],
        ),
        # Module 7: Train test model
        "results/model/All/full.mod.pickle",
        "results/model/All/full.scaler.pickle",
        "results/model/All/full.mod.weights.csv",
        # Module 8: Score variants (CADD scores)
        expand("results/cadd_scores/chr{chr}.tsv.gz", chr=config["chromosomes"]["score"]),
        expand("results/cadd_scores/chr{chr}.tsv.gz.tbi", chr=config["chromosomes"]["score"]),
        "results/cadd_scores/scoring_summary.txt",


##### target rules #####
rule all:
    input:
        # Module 1: Extract ancestor (required for Module 5)
        expand(
            "results/ancestral_seq/{ancestor}/chr{chr}.fa",
            ancestor=config["mark_ancestor"]["name_ancestor"],
            chr=config["chromosomes"]["karyotype"],
        ),
        # Module 2: Derive variants (required for Module 5)
        expand("results/derived_variants/singletons/chr{chr}.vcf", chr=config["chromosomes"]["karyotype"]),
        # Module 3: Simulate variants (required for Module 5)
        expand("results/simulated_variants/trimmed_snps/chr{chr}.vcf", chr=config["chromosomes"]["karyotype"]),
        # Module 4: Summary report (standalone)
        "results/visualisation/stats_report.html",
        # Module 5: Annotate variants (required for Module 6)
        expand(
            "results/annotation/vep/{type}/chr{chr}_vep.tsv",
            chr=config["chromosomes"]["karyotype"],
            type=["simulated", "derived"],
        )
        if should_include_vep()
        else [],
        expand(
            "results/annotation/snpeff/{type}/chr{chr}_snpeff.tsv",
            chr=config["chromosomes"]["karyotype"],
            type=["simulated", "derived"],
        )
        if should_include_snpeff()
        else [],
        # Module 6: Combine annotations (required for Module 7)
        expand("results/annotation/constraint/constraint_chr{chr}.bed", chr=config["chromosomes"]["karyotype"]),
        expand(
            "results/dataset/{type}/chr{chr}_annotated.tsv",
            type=["simulated", "derived"],
            chr=config["chromosomes"]["karyotype"],
        ),
        # Module 7: Train test model (required for Module 8)
        "results/model/All/full.mod.pickle",
        "results/model/All/full.scaler.pickle",
        "results/model/All/full.mod.weights.csv",
        # Module 8: Score variants (final outputs)
        expand("results/cadd_scores/chr{chr}.tsv.gz", chr=config["chromosomes"]["score"]),
        expand("results/cadd_scores/chr{chr}.tsv.gz.tbi", chr=config["chromosomes"]["score"]),
        "results/cadd_scores/scoring_summary.txt",
