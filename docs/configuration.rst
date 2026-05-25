Configuration
=============

Complete guide to configuring the omniCADD pipeline.

Configuration Files
-------------------

The pipeline uses several configuration files:

.. code-block:: text

   config/
   ├── config.yaml                    # Main configuration
   ├── annot_combinations_config.tsv  # Annotation combinations
   ├── annot_processing_config.tsv    # Processing parameters
   └── slurm/                         # Cluster configurations

Main Configuration File
-----------------------

The ``config/config.yaml`` file controls all aspects of the pipeline.

Species Configuration
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: yaml

   # Species identifiers
   species_name: "sus_scrofa"    # Must match alignment and tree
   comparison: "bos_taurus"       # Species for ancestral node
   
   # Display names
   species_display: "Pig"
   comparison_display: "Cattle"

**Important**: Species names must exactly match those in:

* Multi-species alignment (MAF files)
* Phylogenetic tree (Newick file)
* VEP cache (if using Ensembl VEP)

Chromosome Configuration
~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: yaml

   chromosomes:
     # All autosomal chromosomes
     autosomes: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18]
     
     # Complete karyotype (including sex chromosomes)
     karyotype: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, "X", "Y"]
     
     # Chromosomes for model training
     train: [1, 2, 3, 4, 5]
     
     # Chromosomes to score
     score: [1, 2, 3, 4, 5]

**Configuration Tips**:

* **Testing**: Use 1-2 chromosomes: ``train: [1]``, ``score: [1]``
* **Production**: Use all chromosomes for training and scoring

Ancestral Sequence Configuration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Two modes are available for ancestral sequence reconstruction:

Mode 1: Alignment (Recommended)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Extract ancestor from multi-species alignments (e.g., Ensembl EPO):

.. code-block:: yaml

   ancestor:
     mode: "alignment"
     alignment_path: "resources/alignment/"
     species_tree: "resources/tree_43_mammals.nwk"
     ancestor_species: "sus_scrofa"
     comparison_species: "bos_taurus"

**Requirements**:
* EPO or comparable multi-species alignment
* Phylogenetic tree matching alignment
* At least 3 species (reference + 2 outgroups)


Mode 2: Reconstruct (Under Development)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Reconstruct ancestor from pairwise alignment with outgroup:

.. code-block:: yaml

   ancestor:
     mode: "reconstruct"
     outgroup_genome: "resources/outgroup/outgroup_genome.fa.gz"
     alignment_tool: "minimap2"  # or "lastz"
     alignment_params: "-x asm20"  # minimap2 params

**Use case**: Species not in public multi-species alignments

Variant Annotation
^^^^^^^^^^^^^^^^^^

.. code-block:: yaml

   annotation:
     # Variant effect predictor
     variant_annotator: "vep"  # Options: "vep", "snpeff", "both"
     
     # VEP settings
     vep_version: 111             # Ensembl VEP cache version
     vep_auto_install: true       # Auto-download VEP cache
     vep_cache_dir: "~/.vep"      # VEP cache location
     vep_species: "sus_scrofa"    # VEP species name
     vep_assembly: "Sscrofa11.1"  # Assembly version
     
     # Additional VEP options
     vep_plugins: []              # VEP plugins to enable
     vep_custom_annotations: []   # Custom annotation files

**VEP Options**:
* ``vep``: Ensembl Variant Effect Predictor (recommended, when available)
* ``snpeff``: SNPEff annotation (faster, less detailed)
* ``both``: Run both VEP and SNPEff (slower but most comprehensive)

Conservation Annotation
^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: yaml

   annotation:
     conservation:
       # PHAST conservation
       phast: true
       phast_params:
         mod_file: "REV"  # Substitution model
         
       # GERP conservation
       gerp: true
       gerp_params:
         neutral_rate: "auto"  # or specific value

Model Configuration
~~~~~~~~~~~~~~~~~~~

Training Parameters
^^^^^^^^^^^^^^^^^^^

.. code-block:: yaml

   model:
     algorithm: "logistic_regression"
     k_folds: 5  # Cross-validation folds
     random_seed: 42  # For reproducibility
     
     # Grid search parameters
     grid_search:
       C_values: [0.001, 0.01, 0.1, 1, 10]  # Regularization strength
       max_iter: [1000, 2000, 5000]          # Maximum iterations
       
     # Feature selection
     feature_selection:
       enabled: false
       method: "recursive"  # or "lasso"
       n_features: 50

**Grid Search**:
* Pipeline tests all combinations of ``C_values`` and ``max_iter``
* Total models trained: ``len(C_values) × len(max_iter) × k_folds``
* Example: 5 × 3 × 5 = 75 models

**Optimization Tips**:

* **Testing**: Reduce grid search space
  
  .. code-block:: yaml
  
     C_values: [0.1, 1]
     max_iter: [1000]

* **Production**: Use full grid search
* **Large datasets**: May need higher ``max_iter``

Variant Configuration
~~~~~~~~~~~~~~~~~~~~~

Derived Variants
^^^^^^^^^^^^^^^^

.. code-block:: yaml

   derived_variants:
     af_threshold: 0.90           # Minimum allele frequency
     quality_filter: "QUAL>30"    # VCF quality threshold
     depth_filter: "DP>10"        # Depth filter
     singleton_distance: 100      # Min bp between variants

**AF Threshold**:
* **0.90** (default): Fixed/nearly-fixed variants

Simulated Variants
^^^^^^^^^^^^^^^^^^

.. code-block:: yaml

   simulated_variants:
     overestimation_factor: 5    # Mutation rate multiplier
     simulate_indels: false      # Include INDELs (slow, scoring not implemented)
     indel_fraction: 0.1         # If enabled, fraction of INDELs
     
     # Simulation parameters
     mutation_model: "HKY85"     # Substitution model
     transition_transversion: 2.0  # Ti/Tv ratio

**Overestimation Factor**:
* Increases mutation rate to generate more simulated variants
* Higher value = more training data (benign class)
* **Default**: 5x (recommended)
* **Range**: 3-10x

Container Configuration
~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: yaml

   containers:
     use_docker: false           # Use Docker
     use_singularity: true       # Use Singularity/Apptainer
     docker_prefix: "juliahoglund/"  # Docker Hub username
     singularity_cache: "~/.singularity/cache"

**Container Usage**:
* **Local**: Docker recommended
* **HPC**: Singularity/Apptainer required
* **Both disabled**: Use conda environments only

Resource Configuration
~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: yaml

   resources:
     # Default resources for rules
     default_threads: 4
     default_mem_mb: 16000
     
     # Rule-specific overrides
     annotate_vep:
       threads: 8
       mem_mb: 32000
       runtime: "12:00:00"
     
     train_model:
       threads: 16
       mem_mb: 64000
       runtime: "24:00:00"

Annotation Combination Configuration
-------------------------------------

File: ``config/annot_combinations_config.tsv``

Test different annotation combinations:

.. code-block:: text

   Combination         VEP    PHAST  GERP
   full               TRUE   TRUE   TRUE
   vep_only           TRUE   FALSE  FALSE
   conservation_only  FALSE  TRUE   TRUE
   vep_phast          TRUE   TRUE   FALSE
   vep_gerp           TRUE   FALSE  TRUE

**Purpose**: Evaluate feature importance by comparing model performance across annotation subsets.

Annotation Processing Configuration
------------------------------------

File: ``config/annot_processing_config.tsv``

.. code-block:: text

   Parameter              Value
   chunk_size             10000
   min_alignment_depth    3
   max_gap_fraction       0.5
   imputation_method      mean

**Parameters**:
* ``chunk_size``: Alignment positions per chunk (smaller = less memory)
* ``min_alignment_depth``: Minimum species in alignment
* ``max_gap_fraction``: Maximum gap proportion in alignment column
* ``imputation_method``: How to handle missing values (mean, median, mode)

Environment-Specific Configuration
-----------------------------------

Testing Configuration
~~~~~~~~~~~~~~~~~~~~~

For quick testing and development:

.. code-block:: yaml

   # Minimal chromosomes
   chromosomes:
     train: [1]
     score: [1]
   
   # Reduced annotations
   annotation:
     conservation:
       gerp: false  # Disable GERP for speed
   
   # Smaller grid search
   model:
     k_folds: 3
     grid_search:
       C_values: [0.1, 1]
       max_iter: [1000]
   
   # No INDELs
   simulated_variants:
     simulate_indels: false

Production Configuration
~~~~~~~~~~~~~~~~~~~~~~~~

For full genome analysis:

.. code-block:: yaml

   # All chromosomes
   chromosomes:
     train: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18]
     score: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, "X", "Y"]
   
   # Full annotations
   annotation:
     conservation:
       phast: true
       gerp: true
   
   # Complete grid search
   model:
     k_folds: 5
     grid_search:
       C_values: [0.001, 0.01, 0.1, 1, 10]
       max_iter: [1000, 2000, 5000]

Configuration Validation
-------------------------

Check Configuration Syntax
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   # Validate YAML syntax
   python -c "import yaml; yaml.safe_load(open('config/config.yaml'))"
   
   # If valid, no output. If invalid, shows error

Dry Run Validation
~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   # Test configuration without execution
   snakemake -n
   
   # Verbose output
   snakemake -n -p

Check Resource Estimates
~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   # Estimate resources needed
   snakemake --summary
   
   # Show rule execution order
   snakemake --dag | dot -Tpdf > dag.pdf
