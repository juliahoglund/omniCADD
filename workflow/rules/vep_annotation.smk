# VEP Annotation Module
# Contains vep_cache, run_vep, process_vep rules


wildcard_constraints:
    part="[a-zA-Z0-9-]+",


rule vep_cache:
    output:
        directory(config["annotation"]["vep"]["cache"]["directory"]),
    log:
        "results/logs/vep_cache.log",
    conda:
        get_conda_env("annotation")
    threads: get_resource("vep_cache", "threads")
    resources:
        mem_mb=get_resource("vep_cache", "mem_mb"),
        runtime=get_resource("vep_cache", "runtime"),
        time=get_resource("vep_cache", "time"),
        partition=get_resource("vep_cache", "partition"),
    params:
        version_species=config["annotation"]["vep"]["cache"]["install_params"],
    shell:
        "vep_install -a cf -n {params.version_species} -c {output} --CONVERT 2> {log}"


rule run_vep:
    input:
        script="workflow/scripts/vep_annotation/vep.sh",
        vcf="{folder}/{file}.vcf.gz",
        cache=(rules.vep_cache.output if config["annotation"]["vep"]["cache"]["should_install"] == "True" else []),
    output:
        temp("{folder}/{file}_vep_output.tsv"),
    log:
        "results/logs/run_vep/{folder}/{file}.log",
    conda:
        get_conda_env("annotation")
    threads: get_resource("run_vep", "threads")
    resources:
        mem_mb=get_resource("run_vep", "mem_mb"),
        runtime=get_resource("run_vep", "runtime"),
        time=get_resource("run_vep", "time"),
        partition=get_resource("run_vep", "partition"),
    params:
        cache_dir=config["annotation"]["vep"]["cache"]["directory"],
        species_name=config["species_name"],
    shell:
        """
         mkdir -p $(dirname {output})
         chmod +x {input.script} 2>> {log} && {input.script} {input.vcf} {output} {params.cache_dir} {params.species_name} {threads} 2>> {log} && [[ -s {output} ]] 2>> {log}
         """


rule process_vep:
    input:
        script="workflow/scripts/vep_annotation/VEP_process.py",
        vcf=get_vep_source_vcf,
        index=get_vep_source_index,
        vep=get_vep_source_vep_output,
        genome=config["generate_variants"]["reference_genome_wildcard"],
        grantham=config["annotation"]["grantham_matrix"],
    output:
        vep_tsv="results/annotation/vep/{type}/chr{chr}_vep.tsv",
    log:
        "results/logs/process_vep/{type}/chr{chr}.log",
    conda:
        get_conda_env("common")
    threads: get_resource("process_vep", "threads")
    resources:
        mem_mb=get_resource("process_vep", "mem_mb"),
        runtime=get_resource("process_vep", "runtime"),
        time=get_resource("process_vep", "time"),
        partition=get_resource("process_vep", "partition"),
    shell:
        """
         mkdir -p $(dirname {output.vep_tsv})
         python3 {input.script} -v {input.vep} -s {input.vcf} -r {input.genome} -g {input.grantham} -o {output.vep_tsv} 2>> {log}
         """
