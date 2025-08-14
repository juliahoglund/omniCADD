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

def get_folds(excluding = None) -> list:
    """
    Get list of numbers, one for each fold that is to be taken as input.
    :param excluding: optional int(-like), fold to exclude (def None)
    :return: List of numbers, useful for snakemake.expand
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
        script=workflow.source_path(SCRIPTS_7 + "fold_data.py"),
        lib=workflow.source_path(SCRIPTS_7 + "data_helper.py")
    conda: 
        get_conda_env("model")
    priority: 20
    threads: 4
    resources:
        mem_mb = lambda wildcards, attempt: min(128000, int(config["memory"]["dataset_mb"] * attempt)),  # Aggressive scaling for large datasets
        runtime = lambda wildcards, attempt: min(960, 120 * attempt)  # Up to 16 hours for large datasets
    output:
        test = expand("results/dataset/fold_{fold}.npz",
               fold = get_folds()),
        test_m = expand("results/dataset/fold_{fold}.npz.meta.csv.gz",
                 fold = get_folds()),
        test_c = expand("results/dataset/fold_{fold}.npz.columns.csv",
                 fold = get_folds())
    benchmark:
        "logs/benchmarks/fold_data.tsv"
    shell:
        ensure_dirs(*[f"results/dataset/fold_{fold}.npz" for fold in get_folds()]) +
        "python3 {input.script} "
        " -m {input.lib} "
        " -n {threads} "
        " -i {input.derived} {input.simulated} "
        " -o {output.test}"

# Add function to handle fold dependencies
def get_train_folds(fold):
    """Get all folds except the test fold"""
    all_folds = list(range(1, config["model"]["n_folds"] + 1))
    return [f for f in all_folds if f != int(fold)]

rule train_model:
    input:
        test = "results/dataset/fold_{fold}.npz",
        test_m = "results/dataset/fold_{fold}.npz.meta.csv.gz",
        test_c ="results/dataset/fold_{fold}.npz.columns.csv",
        train = lambda wildcards: expand("results/dataset/fold_{fold}.npz",
                                          fold=get_train_folds(wildcards.fold)),
        train_m = lambda wildcards: expand("results/dataset/fold_{fold}.npz.meta.csv.gz",
                                            fold=get_train_folds(wildcards.fold)),
        train_c = lambda wildcards: expand("results/dataset/fold_{fold}.npz.columns.csv",
                                            fold=get_train_folds(wildcards.fold)),
        script = workflow.source_path(SCRIPTS_7 + "train_model.py"),
        lib = workflow.source_path(SCRIPTS_7 + "data_helper.py")
    params:
        c = config["model"]["test_params"]["c"],
        max_iter = config["model"]["test_params"]["max_iter"],
        file_pattern = "results/model/{cols}/fold_{fold}_[C]C_[ITER]iter.mod",
        sel_cols = lambda wildcards: "All" if wildcards.cols == "All" else \
            config["model"]["column_subsets"][wildcards.cols]
    conda:
        get_conda_env("model")
    priority: 20
    resources:
        mem_mb = lambda wildcards, attempt: min(64000, config["memory"]["dataset_mb"] * attempt),
        runtime = lambda wildcards, attempt: min(720, 60 * attempt)
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
    benchmark:
        "logs/benchmarks/train_model_{cols}_fold_{fold}.tsv"
    log:
        "results/logs/model/train_{cols}_fold_{fold}.log"
    shell:
         '''
         ''' + ensure_dirs("results/model/{wildcards.cols}", "results/logs/model", "logs/benchmarks") + '''
         
         python3 {input.script} \
         -m {input.lib} \
         --train {input.train} \
         --test {input.test} \
         --columns {params.sel_cols} \
         -c {params.c} \
         -i {params.max_iter} \
         --file-pattern {params.file_pattern} \
         -n {threads} \
         --save-weights \
         --save-scaler {output.scaler} > {log} 2>&1
         '''

rule final_model:
    input:
          train=expand("results/dataset/fold_{fold}.npz",
                       fold=get_folds()),
          train_m=expand("results/dataset/fold_{fold}.npz.meta.csv.gz",
                         fold=get_folds()),
          train_c=expand("results/dataset/fold_{fold}.npz.columns.csv",
                         fold=get_folds()),
          script=workflow.source_path(SCRIPTS_7 + "train_model.py"),
          lib=workflow.source_path(SCRIPTS_7 + "data_helper.py")
    params:
        c=config["model"]["final_params"]["c"],
        max_iter=config["model"]["final_params"]["max_iter"],
        file_pattern="results/model/{cols}/full.mod",
        sel_cols= lambda wildcards: "All" if wildcards.cols == "All" else \
            config["model"]["column_subsets"][wildcards.cols]
    conda:
         get_conda_env("model")
    priority: 20
    resources:
        mem_mb=lambda wildcards, attempt: min(96000, config["memory"]["dataset_mb"] * attempt),  # Final model needs lots of memory
        runtime=lambda wildcards, attempt: min(1440, 180 * attempt)  # Up to 24 hours
    output:
          model="results/model/{cols}/full.mod.pickle",
          scaler="results/model/{cols}/full.scaler.pickle",
          weights="results/model/{cols}/full.mod.weights.csv"
    benchmark:
        "logs/benchmarks/final_model_{cols}.tsv"
    log:
        "results/logs/model/final_model_{cols}.log"
    shell:
        '''
        ''' + ensure_dirs("results/model/{wildcards.cols}", "results/logs/model", "logs/benchmarks") + '''
        
        python3 {input.script} \\
         -m {input.lib} \\
         --train {input.train} \\
         --columns {params.sel_cols} \\
         -c {params.c} \\
         -i {params.max_iter} \\
         --file-pattern {params.file_pattern} \\
         --save-weights \\
         --save-scaler {output.scaler} > {log} 2>&1
         '''

rule evaluate_models:
    input:
        models = expand("results/model/{{cols}}/fold_{fold}_{c}C_{iter}iter.mod.pickle",
                       fold=get_folds(),
                       c=config["model"]["test_params"]["c"],
                       iter=config["model"]["test_params"]["max_iter"]),
        stats = expand("results/model/{{cols}}/fold_{fold}_{c}C_{iter}iter.mod.stats.txt",
                      fold=get_folds(),
                      c=config["model"]["test_params"]["c"],
                      iter=config["model"]["test_params"]["max_iter"]),
        script = workflow.source_path(SCRIPTS_7 + "evaluate_models.py")
    conda:
        get_conda_env("model")
    resources:
        mem_mb=4000,
        runtime=60
    output:
        summary="results/model/{cols}/model_evaluation_summary.tsv",
        best_params="results/model/{cols}/best_parameters.json"
    benchmark:
        "logs/benchmarks/evaluate_models_{cols}.tsv"
    log:
        "results/logs/model/evaluate_models_{cols}.log"
    shell:
        ensure_dirs("results/model/{wildcards.cols}", "results/logs/model", "logs/benchmarks") +
        "python3 {input.script} "
        "--stats {input.stats} "
        "--summary {output.summary} "
        "--best-params {output.best_params} > {log} 2>&1"



