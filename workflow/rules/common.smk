"""
 Basic helper rules for the other parts of the workflow.
 Configuration loading helper functions are also found here.

 :Author: Job van Schipstal
 :Date: 21-9-2023
"""

"""
Global wildcard constraints, ease matching of wildcards in rules.
Chr is constrained to only be numbers or letters.
Name, label and file may not contain /, they may not be sub-folders.
"""
wildcard_constraints:   
     chr="[a-zA-Z0-9]+",
     file="[^/]+",
     label="[^/]+",
     name="[^/]+"

"""
Bgzip_tabix combines bgzip and the tabix rule to reduce overhead.
This prioritises the combined rule over just tabix,
desired since they produce the same output.
Bgzip_validation_variants has the highest priority since 
it also handles moving the variants into the results folder.
"""


# ruleorder: get_validation_variants > bgzip_tabix > tabix


"""  
 Unzip MAF files since the tool needs a seekable file,
 it is more effective to not have to decompress multiple times.
 Marked as temp so automatically deleted when longer needed, leaving only the compressed original.
"""
rule unzip_maf:
     input:
          "{folder}/{file}.maf.gz"
     output:
          "{folder}/{file}.maf"  # TODO Mark TEMP after testing
     shell:
          "gzip -dc {input} > {output}"

"""
Counts variants in a VCF file, by counting all lines not starting with a comment or whitespace
"""
rule count_variants:
     input:
          "{folder}/{file}.vcf"
     output:
          "{folder}/{file}.vcf.count"
     shell:
          "grep -c '^[^#\S]' {input} > {output}"

"""
Sums counts of all per-chromosome VCF files in a folder.
Since this is only used for training we sum for the training chromosomes.
"""
rule add_counts:
     input:
          expand("{{folder}}/chr{chr}.vcf.count",
               chr=config["chromosomes"]["karyotype"])
     output:
          report("{folder}/total.count", category="Logs")
     shell:
         '''
         cat {input} | awk "{{s+=\$1}} END {{print s}}" > {output}
         '''

"""
Compress VCF file using bgzip and index using tabix
"""
rule bgzip_tabix:
     input:
          "{folder}/{file}.vcf"
     conda:
          "../envs/mainpython.yml"
     output:
          vcf=temp("{folder}/{file}.vcf.gz"),
          index=temp("{folder}/{file}.vcf.gz.tbi")
     shell:
          "bgzip -c {input} > {output.vcf} && "
          "tabix -p vcf -f {output.vcf}"

"""
Index VCF using tabix
"""
rule tabix:
     input:
          "{folder}/{file}.vcf.gz"
     conda:
          "../envs/mainpython.yml"
     output:
          temp("{folder}/{file}.vcf.gz.tbi")
     shell:
          "tabix -p vcf -f {input}"


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
