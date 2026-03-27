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
:Date: 2026-03-17
"""


wildcard_constraints:
    cols="[^/]+",


# Loads the different dataset chunks and merges them.
# The dataset is then split into n_folds which are each written to disk
rule fold_data:
    input:
        derived=expand("results/dataset/derived/chr{chr}.npz", chr=config["chromosomes"]["train"]),
        derived_m=expand(
            "results/dataset/derived/chr{chr}.npz.meta.csv.gz",
            chr=config["chromosomes"]["train"],
        ),
        derived_c=expand(
            "results/dataset/derived/chr{chr}.npz.columns.csv",
            chr=config["chromosomes"]["train"],
        ),
        simulated=expand(
            "results/dataset/simulated/chr{chr}.npz",
            chr=config["chromosomes"]["train"],
        ),
        simulated_m=expand(
            "results/dataset/simulated/chr{chr}.npz.meta.csv.gz",
            chr=config["chromosomes"]["train"],
        ),
        simulated_c=expand(
            "results/dataset/simulated/chr{chr}.npz.columns.csv",
            chr=config["chromosomes"]["train"],
        ),
        script="workflow/scripts/train_test_model/fold_data.py",
        lib="workflow/scripts/utilities/data_helper.py",
    output:
        test=expand("results/dataset/fold_{fold}.npz", fold=get_folds()),
        test_m=expand("results/dataset/fold_{fold}.npz.meta.csv.gz", fold=get_folds()),
        test_c=expand("results/dataset/fold_{fold}.npz.columns.csv", fold=get_folds()),
    log:
        "results/logs/train_test_model/fold_data.log",
    benchmark:
        "logs/benchmarks/fold_data.tsv"
    priority: 20
    conda:
        get_conda_env("model")
    threads: get_resource("fold_data", "threads")
    resources:
        mem_mb=lambda wildcards, attempt: min(
            get_resource("fold_data", "mem_mb") * 2, get_resource("fold_data", "mem_mb") * attempt
        ),
        time=get_resource("fold_data", "time"),
        partition=get_resource("fold_data", "partition"),
    shell:
        """
        mkdir -p results/dataset
        python3 {input.script} -m {input.lib} -n {threads} -i {input.derived} {input.simulated} -o {output.test} 2> {log}
        """


rule train_model:
    input:
        test="results/dataset/fold_{fold}.npz",
        test_m="results/dataset/fold_{fold}.npz.meta.csv.gz",
        test_c="results/dataset/fold_{fold}.npz.columns.csv",
        train=lambda wildcards: expand("results/dataset/fold_{fold}.npz", fold=get_train_folds(wildcards.fold)),
        train_m=lambda wildcards: expand(
            "results/dataset/fold_{fold}.npz.meta.csv.gz",
            fold=get_train_folds(wildcards.fold),
        ),
        train_c=lambda wildcards: expand(
            "results/dataset/fold_{fold}.npz.columns.csv",
            fold=get_train_folds(wildcards.fold),
        ),
        script="workflow/scripts/train_test_model/train_model.py",
        lib="workflow/scripts/utilities/data_helper.py",
    output:
        model=expand(
            "results/model/{{cols}}/fold_{{fold}}_{c}C_{iter}iter.mod.pickle",
            c=config["model"]["test_params"]["c"],
            iter=config["model"]["test_params"]["max_iter"],
        ),
        stats=expand(
            "results/model/{{cols}}/fold_{{fold}}_{c}C_{iter}iter.mod.stats.txt",
            c=config["model"]["test_params"]["c"],
            iter=config["model"]["test_params"]["max_iter"],
        ),
        weights=expand(
            "results/model/{{cols}}/fold_{{fold}}_{c}C_{iter}iter.mod.weights.csv",
            c=config["model"]["test_params"]["c"],
            iter=config["model"]["test_params"]["max_iter"],
        ),
        probs=expand(
            "results/model/{{cols}}/fold_{{fold}}_{c}C_{iter}iter.mod.pred.csv.gz",
            c=config["model"]["test_params"]["c"],
            iter=config["model"]["test_params"]["max_iter"],
        ),
        scaler="results/model/{cols}/fold_{fold}.scaler.pickle",
    log:
        "results/logs/model/train_{cols}_fold_{fold}.log",
    benchmark:
        "logs/benchmarks/train_model_{cols}_fold_{fold}.tsv"
    priority: 20
    conda:
        get_conda_env("model")
    threads: get_resource("train_model", "threads")
    resources:
        mem_mb=lambda wildcards, attempt: min(
            get_resource("train_model", "mem_mb") * 2, get_resource("train_model", "mem_mb") * attempt
        ),
        time=get_resource("train_model", "time"),
        partition=get_resource("train_model", "partition"),
    params:
        c=config["model"]["test_params"]["c"],
        max_iter=config["model"]["test_params"]["max_iter"],
        file_pattern="results/model/{cols}/fold_{fold}_[C]C_[ITER]iter.mod",
        sel_cols=lambda wildcards: ("All" if wildcards.cols == "All" else config["model"]["column_subsets"][wildcards.cols]),
    shell:
        """
         mkdir -p results/model/{wildcards.cols} results/logs/model logs/benchmarks
         python3 {input.script} -m {input.lib} --train {input.train} --test {input.test} --columns {params.sel_cols} -c {params.c} -i {params.max_iter} --file-pattern {params.file_pattern} -n {threads} --save-weights --save-scaler {output.scaler} > {log} 2>&1
         """


rule final_model:
    input:
        train=expand("results/dataset/fold_{fold}.npz", fold=get_folds()),
        train_m=expand("results/dataset/fold_{fold}.npz.meta.csv.gz", fold=get_folds()),
        train_c=expand("results/dataset/fold_{fold}.npz.columns.csv", fold=get_folds()),
        script="workflow/scripts/train_test_model/train_model.py",
        lib="workflow/scripts/utilities/data_helper.py",
    output:
        model="results/model/{cols}/full.mod.pickle",
        scaler="results/model/{cols}/full.scaler.pickle",
        weights="results/model/{cols}/full.mod.weights.csv",
    log:
        "results/logs/model/final_model_{cols}.log",
    benchmark:
        "logs/benchmarks/final_model_{cols}.tsv"
    priority: 20
    conda:
        get_conda_env("model")
    threads: get_resource("final_model", "threads")
    resources:
        mem_mb=lambda wildcards, attempt: min(
            get_resource("final_model", "mem_mb") * 2, get_resource("final_model", "mem_mb") * attempt
        ),
        time=get_resource("final_model", "time"),
        partition=get_resource("final_model", "partition"),
    params:
        c=config["model"]["final_params"]["c"],
        max_iter=config["model"]["final_params"]["max_iter"],
        file_pattern="results/model/{cols}/full.mod",
        sel_cols=lambda wildcards: ("All" if wildcards.cols == "All" else config["model"]["column_subsets"][wildcards.cols]),
    shell:
        """
        mkdir -p results/model/{wildcards.cols} results/logs/model logs/benchmarks
        python3 {input.script} -m {input.lib} --train {input.train} --columns {params.sel_cols} -c {params.c} -i {params.max_iter} --file-pattern {params.file_pattern} --save-weights --save-scaler {output.scaler} > {log} 2>&1
        """


rule evaluate_models:
    input:
        models=expand(
            "results/model/{{cols}}/fold_{fold}_{c}C_{iter}iter.mod.pickle",
            fold=get_folds(),
            c=config["model"]["test_params"]["c"],
            iter=config["model"]["test_params"]["max_iter"],
        ),
        stats=expand(
            "results/model/{{cols}}/fold_{fold}_{c}C_{iter}iter.mod.stats.txt",
            fold=get_folds(),
            c=config["model"]["test_params"]["c"],
            iter=config["model"]["test_params"]["max_iter"],
        ),
        script="workflow/scripts/train_test_model/evaluate_models.py",
    output:
        summary="results/model/{cols}/model_evaluation_summary.tsv",
        best_params="results/model/{cols}/best_parameters.json",
    log:
        "results/logs/model/evaluate_models_{cols}.log",
    benchmark:
        "logs/benchmarks/evaluate_models_{cols}.tsv"
    conda:
        get_conda_env("model")
    threads: get_resource("evaluate_models", "threads")
    resources:
        mem_mb=get_resource("evaluate_models", "mem_mb"),
        time=get_resource("evaluate_models", "time"),
        partition=get_resource("evaluate_models", "partition"),
    shell:
        """
        mkdir -p results/model/{wildcards.cols} results/logs/model logs/benchmarks
        python3 {input.script} --stats {input.stats} --summary {output.summary} --best-params {output.best_params} > {log} 2>&1
        """
