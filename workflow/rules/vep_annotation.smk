# VEP Annotation Module
# Contains vep_cache, run_vep, process_vep rules

wildcard_constraints:
    part="[a-zA-Z0-9-]+",

rule vep_cache:
    params:
        version_species=config["annotation"]["vep"]["cache"]["install_params"],
    log:
        "results/logs/annotate_vars/vep_cache.log",
    output:
        directory(config["annotation"]["vep"]["cache"]["directory"]),
    conda:
        get_conda_env("annotation")
    shell:
        "vep_install -a cf -n {params.version_species} -c {output} --CONVERT 2> {log}"

rule run_vep:
    input:
        script="../scripts/vep_annotation/vep.sh",
        vcf="{folder}/{file}.vcf.gz",
        cache=(
            rules.vep_cache.output
            if config["annotation"]["vep"]["cache"]["should_install"] == "True"
            else []
        ),
    params:
        cache_dir=config["annotation"]["vep"]["cache"]["directory"],
        species_name=config["species_name"],
    log:
        "{folder}/{file}_vep.log",
    conda:
        get_conda_env("annotation")
    threads: 2
    output:
        temp("{folder}/{file}_vep_output.tsv"),
    shell:
        """
         mkdir -p $(dirname {output})
         chmod +x {input.script} 2>> {log} && {input.script} {input.vcf} {output} {params.cache_dir} {params.species_name} {threads} 2>> {log} && [[ -s {output} ]] 2>> {log}
         """

rule process_vep:
    input:
        script="../scripts/vep_annotation/VEP_process.py",
        vcf="{folder}/chr{chr}.vcf.gz",
        index="{folder}/chr{chr}.vcf.gz.tbi",
        vep="{folder}/chr{chr}_vep_output.tsv",
        genome=config["generate_variants"]["reference_genome_wildcard"],
        grantham=[
            config["annotation"]["grantham_matrix"],
        ],
    params:
        output_type=lambda wildcards: (
            "derived" if "derived_variants" in wildcards.folder else "simulated"
        ),
    log:
        "{folder}/chr{chr}_process_vep.log",
    conda:
        get_conda_env("common")
    output:
        vep_tsv="{folder}/chr{chr}_vep.tsv",
    shell:
        """
         mkdir -p $(dirname {output}) results/annotation/vep/{params.output_type}
         python3 {input.script} -v {input.vep} -s {input.vcf} -r {input.genome} -g {input.grantham} -o {output.vep_tsv} 2>> {log} && cp {output.vep_tsv} results/annotation/vep/{params.output_type}/chr{wildcards.chr}_vep.tsv 2>> {log}
         """