# -*- snakemake -*-

"""
 Module that applies a logistic regression model to the dataset generated.
 Trains a model which is then validated using a hold-out test set
 and databases of known variants.
 Additionally, the model is used to score all variants,
 which are then used to generate the whole-genome CADD scores.

 :Author: Job van Schipstal
 :Date: 23-10-2023

 Scripts are based upon the work of Christian Gross,
 but written for scikit-learn instead of turi create.

 :Extension and modification: Julia Höglund
 :Date: 30-4-2024
"""

wildcard_constraints:
    cols="[^/]+"

# what are the folds?
def get_folds(excluding = None) -> list:
    """
    Get list of numbers, one for each fold that is to be taken as input.
    :param excluding: optional int(-like), fold to exclude (def None)
    :return: List of numbers, usefull for snakemake.expand
    """
    folds = list(range(config["model"]["n_folds"]))
    if excluding:
        excluding = int(excluding)
        if excluding in folds:
            folds.remove(excluding)
    return folds

"""
Loads the different dataset chunks and merges them.
The dataset is then split into n_folds which are each written to disk
"""
# change to fit file names later
rule fold_data:
    input:
        derived = expand("results/dataset/derived/chr{chr}.npz",
                  chr=config["chromosomes"]["train"]),
        derived_m = expand("results/dataset/derived/chr{chr}.npz.meta.csv.gz",
                    chr=config["chromosomes"]["train"]),
        derived_c = expand("results/dataset/derived/chr{chr}.npz.columns.csv",
                    chr=config["chromosomes"]["train"]),
        simulated = expand("results/dataset/simulated/chr{chr}.npz",
                    chr=config["chromosomes"]["train"]),
        simulated_m = expand("results/dataset/simulated/chr{chr}.npz.meta.csv.gz",
                      chr=config["chromosomes"]["train"]),
        simulated_c = expand("results/dataset/simulated/chr{chr}.npz.columns.csv",
                      chr=config["chromosomes"]["train"]),

        script=workflow.source_path(MODEL_P + "fold_data.py"),
        lib=workflow.source_path(MODEL_P + "data_helper.py")
    conda: 
        "../envs/model.yml"
    priority: 20
    threads: 4
    #resources:
    #    mem_mb = int(config["dataset_memory_mb"] * 2)
    output:
        test = expand("results/dataset/fold_{fold}.npz",
               fold = get_folds()),
        test_m = expand("results/dataset/fold_{fold}.npz.meta.csv.gz",
                 fold = get_folds()),
        test_c = expand("results/dataset/fold_{fold}.npz.columns.csv",
                 fold = get_folds())
    shell:
        "python3 {input.script} "
        " -m {input.lib} "
        " -n {threads} "
        " -i {input.derived} {input.simulated} "
        " -o {output.test}"

# change to fit file names later
rule train_model:
    input:
        test = "results/dataset/fold_{fold}.npz",
        test_m = "results/dataset/fold_{fold}.npz.meta.csv.gz",
        test_c ="results/dataset/fold_{fold}.npz.columns.csv",
          
        train = lambda wildcards: expand("results/dataset/fold_{fold}.npz",
                                          fold=get_folds(wildcards.fold)),
        train_m = lambda wildcards: expand("results/dataset/fold_{fold}.npz.meta.csv.gz",
                                            fold=get_folds(wildcards.fold)),
        train_c = lambda wildcards: expand("results/dataset/fold_{fold}.npz.columns.csv",
                                            fold=get_folds(wildcards.fold)),

        sel_cols = lambda wildcards: [] if wildcards.cols == "All" else \
                     config["model"]["column_subsets"][wildcards.cols],
        script = workflow.source_path(SCRIPTS_7 + "train_model.py"),
        lib = workflow.source_path(MODEL_P + "data_helper.py")
    params:
        c = config["model"]["test_params"]["c"],
        max_iter = config["model"]["test_params"]["max_iter"],
        file_pattern = "results/model/{cols}/fold_{fold}_[C]C_[ITER]iter.mod",
        sel_cols = lambda wildcards: "All" if wildcards.cols == "All" else \
            config["model"]["column_subsets"][wildcards.cols]
    conda:
        "../envs/mainpython.yml"
    priority: 20
    resources:
        mem_mb = config["dataset_memory_mb"]
    threads: len(config["model"]["test_params"]["c"]) * \
             len(config["model"]["test_params"]["max_iter"])
    output:
        model = expand("results/model/{{cols}}/fold_{{fold}}_{c}C_{iter}iter.mod.pickle",
                       c=config["model"]["test_params"]["c"],
                       iter=config["model"]["test_params"]["max_iter"]),
        stats = expand("results/model/{{cols}}/fold_{{fold}}_{c}C_{iter}iter.mod.stats.txt",
                       c=config["model"]["test_params"]["c"],
                       iter=config["model"]["test_params"]["max_iter"]),
        weights = expand("results/model/{{cols}}/fold_{{fold}}_{c}C_{iter}iter.mod.weights.csv",
                         c=config["model"]["test_params"]["c"],
                         iter=config["model"]["test_params"]["max_iter"]),
        probs = expand("results/model/{{cols}}/fold_{{fold}}_{c}C_{iter}iter.mod.pred.csv.gz",
                       c=config["model"]["test_params"]["c"],
                       iter=config["model"]["test_params"]["max_iter"]),
        scaler = "results/model/{cols}/fold_{fold}.scaler.pickle"
    shell:
         """python3 {input.script} \
         -m {input.lib} \
         --train {input.train} \
         --test {input.test} \
         --columns {params.sel_cols} \
         -c {params.c} \
         -i {params.max_iter} \
         --file-pattern {params.file_pattern} \
         -n {threads} \
         --save-weights \
         --save-scaler {output.scaler}"""